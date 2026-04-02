// ---------------------------------------------------------------------------
//
//      The ab-initio based multiscale library
//
//                  S / P H I / n X
//
//      Copyright:  Max-Planck-Institut fuer Eisenforschung GmbH
//                  40237 Duesseldorf, Germany
//
//      Contact:    https://sxrepo.mpie.de
//      Authors:    see sphinx/AUTHORS
//      License:    see sphinx/LICENSE
//
// ---------------------------------------------------------------------------
#include <SxFourierInterpol.h>
#include <SxConstants.h>

SxFourierInterpol::SxFourierInterpol (SxBinIO &io)
{
   SX_CHECK (io.mode == SxBinIO::BINARY_READ_ONLY);

   try {
      if (! io.contains ("cell"))  {
         cout << endl << "'cell' is missing in netCDF file '";
         cout << io.filename << "'." << endl;
         cout << "This doesn't seem to be a netCDF mesh file." << endl;
         SX_QUIT;
      }
      cell.read (io);
      SxVector3<int> mesh;
      io.read ("dim", &mesh);

      int meshSize = mesh.product ();
      SxVector<double> realSpace (meshSize);
      io.read ("mesh", &realSpace, meshSize);

      // --- Fourier transform to reciprocal space
      SxFFT3d fft(SxFFT::Reverse, mesh(0), mesh(1), mesh(2), cell.volume);
      SxVector<SX_FFT_COMPLEX> in (realSpace), out(meshSize);
      fft.fftReverse (meshSize, in.elements, out.elements);
      meshInG = SxVector<PrecCoeffG> (out);

      // --- set up G-vectors
      gVec.reformat (meshSize,3);
      SxVector3<int> relG;
      int dim;
      for (int i = 0; i < meshSize; i++)  {
         relG = fft.mesh.getMeshVec (i, SxMesh3D::Origin);
         for (dim = 0; dim < 3; dim++)  {
            gVec(i,dim) = relG(dim);
         }
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

}

void SxFourierInterpol::print (ostream &out)
{
   int dim;
   for (int i = 0; i < meshInG.getSize (); i++)  {
      for (dim = 0; dim < 3; dim++)
         out << gVec(i,dim) << ' ';
      out << meshInG(i) << endl;
   }
}

SxVector<double> SxFourierInterpol::interpolate (const SxAtomicStructure &rVec)
{
   SX_CHECK (meshInG.getSize () > 0, meshInG.getSize ());
   int ng = (int)meshInG.getSize ();
   int nr = rVec.getNAtoms ();
   SX_CHECK (nr > 0, nr);

   // alternative 1:
   // e^iGr = e^i(gx rx) * e^(gy ry) * e^i(gz rz)
   // use that gx, gy, gz is integer and setup
   // exp_x (gx,r), exp_y(gy,r), exp_z(gz,r)
   // and then calculate
   // expikr (ir, ig) = exp_x(gx,r) * exp_y(gy,r) * exp_z(gz,r)

   gVec = gVec.transpose ();
   int maxX = 0, maxY = 0, maxZ = 0;
   for (int ig = 0; ig < ng; ++ig)  {
      if (abs(gVec(0,ig)) > maxX) maxX = abs(gVec(0,ig));
      if (abs(gVec(1,ig)) > maxY) maxY = abs(gVec(1,ig));
      if (abs(gVec(2,ig)) > maxZ) maxZ = abs(gVec(2,ig));
   }
   SxCell recCell = cell.getReciprocalCell ();
   SxVector<SxComplex16> phiX(2*maxX+1), phiY(2*maxY+1), phiZ(2*maxZ+1);

   SxVector<double> res(nr);
   double sum, gr;
   for (int ir = 0; ir < nr; ++ir)  {
      gr = rVec(ir) ^ recCell.basis(0);
      for (int i = 0; i <= 2*maxX; ++i)
         sincos(gr * (i-maxX), &phiX(i).im, &phiX(i).re);
      gr = rVec(ir) ^ recCell.basis(1);
      for (int i = 0; i <= 2*maxY; ++i)
         sincos(gr * (i-maxY), &phiY(i).im, &phiY(i).re);
      gr = rVec(ir) ^ recCell.basis(2);
      for (int i = 0; i <= 2*maxZ; ++i)
         sincos(gr * (i-maxZ), &phiZ(i).im, &phiZ(i).re);
      sum = 0.;
      for (int ig = 0; ig < ng; ++ig)
         sum += (  phiX(gVec(0,ig) + maxX) * phiY(gVec(1,ig) + maxY)
                 * phiZ(gVec(2,ig) + maxZ) * meshInG(ig)).re;
      res(ir) = sum / sqrt(cell.volume);
   }

   gVec = gVec.transpose ();
   return res;

   /*
   // --- alternative 2: use matrix-matrix multiplication

   SxVector<double> gVecCart, gr; // ng x 3
   SxVector<PrecCoeffG> expikr(nr, ng);
   gVecCart = SxVector<double> (gVec) ^ SxVector<double> (TWO_PI * cell.inv);
   // (cout << "G ^ r...").flush ();
   gr = (gVecCart ^ rVec.coords).transpose (); // TODO: avoid transpose

   // (cout << "done\n").flush ();
   SxVector<double>::Iterator grIt = gr.begin ();
   SxVector<PrecCoeffG>::Iterator expikrIt = expikr.begin ();
   // (cout << "exp(iGr)...").flush ();
   for (int i = 0; i < gr.getSize (); i++, grIt++, expikrIt++)  {
      (*expikrIt).re = cos(*grIt);
      (*expikrIt).im = sin(*grIt);
   }
   SX_CHECK (grIt == gr.end ());
   SX_CHECK (expikrIt == expikr.end ());
   // (cout << "done\n").flush ();

   SX_CHECK (cell.volume > 0., cell.volume);

   // (cout << "sum c(G) * exp(iGr)...").flush ();
   SxVector<double> result = (expikr ^ meshInG) / sqrt(cell.volume);
   // (cout << "done\n").flush ();

   return result;
   */
}

void SxFourierInterpol::condense (double threshold)
{
   int meshSize = (int)meshInG.getSize ();
   int newMeshSize = 0;

   SxVector<int> map(meshSize);

   SxVector<PrecCoeffG>::Iterator it = meshInG.begin ();
   SxVector<int>::Iterator mapIt = map.begin ();

   for (int i = 0; i < meshSize; i++, ++it)
      if ((*it).absSqr () >= threshold)  {
         *mapIt++ = i;
         newMeshSize++;
      }

   SxVector<PrecCoeffG> newMeshInG(newMeshSize);
   SxVector<int> newGVec(newMeshSize,3);

   it = newMeshInG.begin ();
   mapIt = map.begin ();
   for (int i = 0; i < newMeshSize; i++, mapIt++, it++)  {
      *it = meshInG(*mapIt);
      newGVec(i,0) = gVec(*mapIt,0);
      newGVec(i,1) = gVec(*mapIt,1);
      newGVec(i,2) = gVec(*mapIt,2);
   }

   meshInG = newMeshInG;
   gVec = newGVec;

}

