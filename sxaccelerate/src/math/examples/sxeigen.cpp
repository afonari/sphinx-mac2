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

#include <SxUtil.h>
#include <SxVector.h>
#include <SxSymMatrix.h>
#include <SxMathLib.h>
#include <SxConstants.h>  /* for complex I */
#include <SxEigensystem.h>
#include <stdio.h>
#include <iostream>

#define TYPE_M1 SxComplex8
#define TYPE_M2 SxComplex16
#define TYPE_M3 float
#define TYPE_M4 double

/** 
  \example eigen.cpp
  \author  Sixten Boeck
  */
int main ()
{
   initSPHInXMath ();

   int n = 4;
   SxVector<TYPE_M1>    m1(n,n);
   SxSymMatrix<TYPE_M2> m2(n);
   SxVector<TYPE_M3>    m3(n,n);
   SxSymMatrix<TYPE_M4> m4(n);

   // --- set up matrix m1 and m3
   int r, c;
   for (r=0; r < n; r++)  {
      for (c=r; c < n; c++)  {
         m1(r,c) = TYPE_M2(sqrt((double)(r*c+1)),
                           cbrt((double)(r*c+1)));
         m1(c,r) = m1(r,c).conj();
         m3(r,c) = ((double)(r+c+1));
         m3(c,r) = m3(r,c);
         if (r==c)  m1(r,c).im = 0.;
      }
   }
   // --- initialize matrix m2/4. m2/4 = m1/3, but trigonally packed
   //     (use only upper-right triangle)
   for (r=0; r < n; r++)
      for (c=r; c < n; c++)  {
         m2(r,c) = TYPE_M2(sqrt((double)(r*c+1)),
                           cbrt((double)(r*c+1)));
         m4(r,c) = sqrt((double)(r*c+1));
         if (r==c)  m2(r,c).im = 0.;
      }
               
   // --- print matrices
   cout << "m1\n" << m1 << endl;
   cout << "m2\n" << m2 << endl;
   cout << "m3\n" << m3 << endl;
   cout << "m4\n" << m4 << endl;

   cout << "m1 ^ m1 =\n" << (m1 ^ m1) << endl; 
   cout << "m1 ^ m2 =\n" << (SxVector<TYPE_M2>(m1) ^ m2.expand()) << endl;
   cout << "m3 ^ m1 =\n" << (SxVector<TYPE_M1>(m3) ^ m1) << endl;
   cout << "m2 ^ m4 =\n" << (m2 ^ SxVector<TYPE_M2>(m4.expand())) << endl;

   // --- compute eigensystems
   SxEigensystem<TYPE_M1> eig1(m1);
   SxEigensystem<TYPE_M3> eig3(m3);
   SxSymEigensystem<TYPE_M2> eig2(m2);
   SxSymEigensystem<TYPE_M4> eig4(m4);

   cout << "e1: "; eig1.print(true);
   cout << "e2: "; eig2.print(true);
   cout << "e3: "; eig3.print(true);
   cout << "e4: "; eig4.print(true);

   // --- check eigensystem
   printf ("Check eigensystem\n");
   for (r=0; r < n; r++)  {
      cout << "||m1^V("<<r<<") - v("<<r<<")*V("<<r<<")|| = "
        << (   (m1^eig1.vecs.colRef(r)) 
            -  eig1.vals(r)*eig1.vecs.colRef(r)).absSqr().sum()
        << endl;
   }
   for (r=0; r < n; r++)  {
      cout << "||m2^V("<<r<<") - v("<<r<<")*V("<<r<<")|| = "
        << (   (m2.expand()^eig2.vecs.colRef(r)) 
            -  eig2.vals(r)*eig2.vecs.colRef(r)).absSqr().sum()
        << endl;
   }
   for (r=0; r < n; r++)  {
      cout << "||m3^V("<<r<<") - v("<<r<<")*V("<<r<<")|| = "
        << (   (m3^eig3.vecs.colRef(r)) 
            -  eig3.vals(r)*eig3.vecs.colRef(r)).absSqr().sum()
        << endl;
   }
   for (r=0; r < n; r++)  {
      cout << "||m4^V("<<r<<") - v("<<r<<")*V("<<r<<")|| = "
        << (   (m4.expand()^eig4.vecs.colRef(r)) 
            -  eig4.vals(r)*eig4.vecs.colRef(r)).absSqr().sum()
        << endl;
   }
   cout << endl;

   // --- check norm of eigenvectors
   printf ("Check norm\n");
   for (r=0; r < n; r++)  {
      cout << "||V("<<r<<")|| = "
        << dot((eig1.vecs.colRef(r)), eig1.vecs.colRef(r)) << endl;
   }
   for (r=0; r < n; r++)  {
      cout << "||V("<<r<<")|| = "
        << dot((eig3.vecs.colRef(r)), eig3.vecs.colRef(r)) << endl;
   }


   // --- Cholesky decomposition, M = L * L^  (M symmetric)
   m1 = SxVector<TYPE_M1> (
           SxList<TYPE_M1> () << 10 <<  6       << 5
                              <<  6 << 22       << 1.+2.*I
                              <<  5 <<  1.-2.*I << 3
        ).transpose ();
   m1.reshape (3,3);
   cout << "m1 = \n" << m1;
   SxVector<TYPE_M1> L = m1.choleskyDecomposition();
   cout << "L = \n" << L;
   cout << "L ^ L^* - M = " 
        << ((L ^ L.adjoint()) - m1).absSqr().sum() << endl;


   return 0;
}

