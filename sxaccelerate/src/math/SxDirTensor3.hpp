// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------

template<class T>
SxDirTensor3<T>::SxDirTensor3 ()
{
   n1 = n2 = n3 = -1;
}
   
template<class T>
SxDirTensor3<T>::SxDirTensor3 (int n1_, int n2_, int n3_)
{
   n1 = n1_;  n2 = n2_;  n3 = n3_;

   tensor.resize (n1);

   for (int i1 = 0; i1 < n1; i1++)  {
      tensor(i1) = SxVector<T> (n2,n3);
//      tensor(i1).reformat (n2,n3);
   }
}

template<class T>
SxDirTensor3<T>::SxDirTensor3 (const SxDirTensor3<T> &t)
{
   n1 = t.n1, n2 = t.n2, n3 = t.n3;

   SX_CHECK (t.tensor.getSize() == n1, t.tensor.getSize(), n1);

   tensor.resize (n1);

   for (int i1 = 0; i1 < n1; i1++)  {
      tensor(i1) = SxVector<T> (n2,n3);
//      tensor(i1).reformat (n2,n3);
      tensor(i1).copy ( t(i1) );
   }
}

template<class T>
void SxDirTensor3<T>::reformat (int n1_, int n2_, int n3_)
{
   n1 = n1_;  n2 = n2_;  n3 = n3_;

   tensor.resize (n1);

   for (int i1 = 0; i1 < n1; i1++)  {
      tensor(i1) = SxVector<T> (n2,n3);
//      tensor(i1).reformat (n2,n3);
   }
}

template<class T>
void SxDirTensor3<T>::set (const SxVector<T> &mat, int i1)
{
   SX_CHECK (i1 >= 0 && i1 < n1, i1, n1);

   tensor(i1) = mat; // copy!
}

template<class T>
void SxDirTensor3<T>::set (SxVector<T> &&mat, int i1)
{
   SX_CHECK (i1 >= 0 && i1 < n1, i1, n1);

   tensor(i1) = std::move (mat);
}

template<class T>
void SxDirTensor3<T>::set (const T &elem, int i1, int i2, int i3)
{
   SX_CHECK (i1 >= 0 && i1 < n1, i1, n1);
   SX_CHECK (i2 >= 0 && i2 < n2, i2, n2);
   SX_CHECK (i3 >= 0 && i3 < n3, i3, n3);

   tensor(i1)(i2,i3) = elem;
}

template<class T>
void SxDirTensor3<T>::setAllElem (const T &in)
{
   SX_CHECK (tensor.getSize() == n1, tensor.getSize(), n1);

   for (int i1 = 0; i1 < n1; i1++)  {
      tensor(i1).set (in);
   }
}

template<class T>
void SxDirTensor3<T>::setDiagonals (const T &in)
{
   SX_CHECK (tensor.getSize() == n1, tensor.getSize(), n1);
   SX_CHECK (n2 == n3, n2, n3);  // otherwise not quadratic

   for (int i1 = 0; i1 < n1; i1++)  {
      tensor(i1).reformat (n2, n2);
      tensor(i1).set (0.);
      tensor(i1).diag ().set (in);
   }
}

template<class T>
SxVector<T> &SxDirTensor3<T>::operator() (int i1)
{
   SX_CHECK (i1 >= 0 && i1 < n1, i1, n1);
   return tensor(i1);
}

template<class T>
const SxVector<T> &SxDirTensor3<T>::operator() (int i1) const
{
   SX_CHECK (i1 >= 0 && i1 < n1, i1, n1);
   return tensor(i1);
}

template<class T>
T &SxDirTensor3<T>::operator() (int i1, int i2, int i3)
{
   SX_CHECK (i1 >= 0 && i1 < n1, i1, n1);
   SX_CHECK (i2 >= 0 && i2 < n2, i2, n2);
   SX_CHECK (i3 >= 0 && i3 < n3, i3, n3);
   return tensor(i1)(i2,i3);
}

template<class T>
const T
&SxDirTensor3<T>::operator() (int i1, int i2, int i3) const
{
   SX_CHECK (i1 >= 0 && i1 < n1, i1, n1);
   SX_CHECK (i2 >= 0 && i2 < n2, i2, n2);
   SX_CHECK (i3 >= 0 && i3 < n3, i3, n3);
   return tensor(i1)(i2,i3);
}

template<class T>
SxDirTensor3<T> SxDirTensor3<T>::operator= (const SxDirTensor3<T> &in)
{
   if (this != &in) tensor = in.tensor; // copy
   return *this;
}

template<class T>
SxDirTensor3<T> SxDirTensor3<T>::operator= (SxDirTensor3<T> &&in)
{
   if (this != &in) tensor = std::move (in.tensor); // copy
   // destroy in
   in.tensor.resize (0);
   in.n1 = in.n2 = in.n3 = 0;
   return *this;
}

template<class T>
void SxDirTensor3<T>::print ()
{
   print (n1);
}

template<class T>
void SxDirTensor3<T>::print (int z)
{
   SX_CHECK (z > 0 && z <= n1, z, n1);

   cout << SX_SEPARATOR;
   if (z < n1)  {
      sxprintf ("SxDirTensor3  ---  %d of %d matrices printed:\n\n", z, n1);
   }

   for (int i1 = 0; i1 < z; i1++)  {
      sxprintf ("i1 = %d:\n\n", i1);
      cout << tensor(i1) << endl;
   }

   cout << SX_SEPARATOR;
}

template<class T>
void SxDirTensor3<T>::print (const SxString &s)
{
   print (s, n1);
}

template<class T>
void SxDirTensor3<T>::print (const SxString &s, int z)
{
   SX_CHECK (z > 0 && z <= n1, z, n1);

   cout << SX_SEPARATOR;

   sxprintf ("SxDirTensor3:  %s  ---  %d of %d matrices printed:\n\n",
             s.ascii(), z, n1);

   for (int i1 = 0; i1 < z; i1++)  {
      sxprintf ("%s --> i1 = %d:\n\n", s.ascii(), i1);
      cout << tensor(i1) << endl;
      if (i1 != z-1)  {
         cout << "- - - - - - - - - - - - - - - - - - - - - - -" << endl;
      }  else  {
         cout << SX_SEPARATOR;
      }
   }
}

template<class T>
void SxDirTensor3<T>::writeDirTensor3 (const SxString &file) const
{
   try  {
      SxBinIO io (file, SxBinIO::BINARY_WRITE_ONLY);
      writeDirTensor3 (io);
      io.setMode (SxBinIO::WRITE_DATA);
      writeDirTensor3 (io);
      io.close ();
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }
}

template<class T>
void SxDirTensor3<T>::writeDirTensor3 (const SxBinIO &io) const
{
   int iElem, nElem = n1 * n2 * n3;
   int i1, i2, i3;
   SxVector<T> elems(nElem);

   // --- encode data: tensor --> one very long vector
   iElem = 0;
   for (i1 = 0; i1 < n1; i1++)  {
      for (i2 = 0; i2 < n2; i2++)  {
         for (i3 = 0; i3 < n3; i3++)  {
            elems(iElem) = tensor(i1)(i2,i3);
            iElem++;
         }
      }
   }

   SX_CHECK (iElem == nElem, iElem, nElem);

   try  {
      // --- write header, create dimensions
      io.addDimension ("n1", n1);
      io.addDimension ("n2", n2);
      io.addDimension ("n3", n3);
      io.addDimension ("nElem", nElem);

      // --- write data
      io.writeVec ("elems", elems, "nElem");
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }
}

template<class T>
void SxDirTensor3<T>::readDirTensor3 (const SxString &file)
{
   try  {
      SxBinIO io (file, SxBinIO::BINARY_READ_ONLY);
      readDirTensor3 (io);
      io.close ();
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }
}

template<class T>
void SxDirTensor3<T>::readDirTensor3 (const SxBinIO &io)
{
   int iElem, nElem = -1;
   int i1, i2, i3, n1In = -1, n2In = -1, n3In = -1;
   SxVector<T> elems;

   try  {
      n1In  = io.getDimension ("n1");
      n2In  = io.getDimension ("n2");
      n3In  = io.getDimension ("n3");
      nElem = io.getDimension ("nElem");

      elems.resize (nElem);

      io.readVec ("elems", &elems, nElem);
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

   if (n1In * n2In * n3In != nElem)  {
      sxprintf ("| Error:  file \"%s\" corrupt, n1 x n2 x n3 = %d, but "
                "nElem = %d\n",
                io.filename.ascii(), n1In * n2In * n3In, nElem);
      SX_QUIT;
   }

   reformat (n1In, n2In, n3In);

   // --- decode data
   iElem = 0;
   for (i1 = 0; i1 < n1; i1++)  {
      for (i2 = 0; i2 < n2; i2++)  {
         for (i3 = 0; i3 < n3; i3++)  {
            tensor(i1)(i2,i3) = elems(iElem);
            iElem++;
         }
      }
   }

   SX_CHECK (iElem == nElem, iElem, nElem);
}

template<class T>
void SxDirTensor3<T>::readTensor3 (const SxString &file)
{
   try  {
      SxBinIO io (file, SxBinIO::BINARY_READ_ONLY);
      readTensor3 (io);
      io.close ();
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }
}

template<class T>
void SxDirTensor3<T>::readTensor3 (const SxBinIO &io)
{
   int iElem, nElem = -1;
   int i1, i2, i3, n1In = -1, n2In = -1, n3In = -1;
   SxVector<T> elems;

   try  {
      n1In  = io.getDimension ("nk");
      n2In  = io.getDimension ("nm");
      n3In  = io.getDimension ("nn");
      nElem = io.getDimension ("nElem");

      elems.resize (nElem);

      io.readVec ("kmnElems", &elems, nElem);
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

   if (n1In * n2In * n3In != nElem)  {
      sxprintf ("| Error: file \"%s\" is corrupt, since nk x nm x nn = %d, "
                "but nElem = %d.\n",
                io.filename.ascii(), n1In * n2In * n3In, nElem);
      SX_QUIT;
   }

   reformat (n1In, n2In, n3In);

   // --- decode data
   iElem = 0;
   for (i1 = 0; i1 < n1; i1++)  {
      for (i2 = 0; i2 < n2; i2++)  {
         for (i3 = 0; i3 < n3; i3++)  {
            tensor(i1)(i2,i3) = elems(iElem);
            iElem++;
         }
      }
   }

   SX_CHECK (iElem == nElem, iElem, nElem);
}

template<class T>
void SxDirTensor3<T>::readUkmn (const SxString &file)
{
   int iElem, nElem = -1;
   int ik, im, in, nkIn = -1, dimIn = -1;
   SxVector<T> kmnElems;

   try  {
      SxBinIO io (file, SxBinIO::BINARY_READ_ONLY);

      nkIn  = io.getDimension ("nk");
      dimIn = io.getDimension ("dim");
      nElem = io.getDimension ("nElem");

      kmnElems.resize (nElem);

      io.readVec ("kmnElems", &kmnElems, nElem);

      if (nkIn * dimIn * dimIn != nElem)  {
         sxprintf ("| Error:  file \"%s\" corrupt, nk x dim x dim = %d, but "
                   "nElem = %d\n",
                   io.filename.ascii(), nkIn * dimIn * dimIn, nElem);
         SX_QUIT;
      }

      io.close ();
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

   reformat (nkIn, dimIn, dimIn);

   // --- decode data
   iElem = 0;
   for (ik = 0; ik < nkIn; ik++)  {
      for (im = 0; im < dimIn; im++)  {
         for (in = 0; in < dimIn; in++)  {
            tensor(ik)(im,in) = kmnElems(iElem);
            iElem++;
         }
      }
   }

   SX_CHECK (iElem == nElem, iElem, nElem);
}
