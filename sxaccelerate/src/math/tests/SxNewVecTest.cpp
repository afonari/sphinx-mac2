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

#include <SxVector.h>
#include <SxEigensystem.h>
#include <SxCLI.h>
#include <SxSymMatrix.h>
#include <SxDiagMat.h>

int main (int argc, char** argv)
{
   SxCLI cli (argc, argv);
   int M = cli.option ("-M","matrix size","matrix size for matmult test")
           .toInt (500);
   cli.finalize ();

   // enable caching of memory allocations
   SX_ALLOC_CACHE;

   // simple vector
   SxVector<double> a(10);
   cout << "sizeof(a) = " << sizeof(a) << endl;
   SX_LOOP(i)
      a(i) = i;
   cout << "a=" << a << endl;

   SxVector<double> b(10);
   SX_LOOP(i)
      b(i) = 2 * i + 1;
   cout << "b=" << b << endl;

   // vector-vector arithmetic
   cout << "a+b=" << (a+b) << endl;
   cout << "a-b=" << (a-b) << endl;
   cout << "a*b=" << (a*b) << endl;
   cout << "a/b=" << (a/b) << endl;

   // a 10x10 matrix
   SxVector<double> mat;
   mat.reformat (10, 10);
   SX_LOOP(i) mat(i) = i+1;

   // a 5x10 matrix
   SxVector<double> mat2;
   mat2.reformat (5, 10);
   SX_LOOP(i) mat2(i) = i+1;

   // submatrix = reference to contiguous subset of rows/columns
   //             with gaps between (contiguous) columns
   SxVecRef<double, SubMatrix> sub = mat.blockRef (0,0,5,10);
   SX_LOOP(i) sub(i) = -sub(i);

   // vector-vector arithmetic with one operand being SubMatrix
   cout << mat << endl;
   cout << "mat+sub=" << (mat2 + sub) << endl;
   cout << "sub+mat=" << (sub+mat2) << endl;
   cout << "mat-sub=" << (mat2 - sub) << endl;
   cout << "sub-mat=" << (sub-mat2) << endl;
   cout << "mat*sub=" << (mat2 * sub) << endl;
   cout << "sub*mat=" << (sub*mat2) << endl;
   cout << "mat/sub=" << (mat2 / sub) << endl;
   cout << "sub/mat=" << (sub/mat2) << endl;

   // vector-vector arithmetic with both operands being SubMatrix
   SxVecRef<double, SubMatrix> sub2 = mat.blockRef (5,0,5,10);
   cout << sub2 << endl;
   cout << "sub+sub=" << (sub + sub2) << endl;
   cout << "sub-sub=" << (sub - sub2) << endl;
   cout << "sub*sub=" << (sub * sub2) << endl;
   cout << "sub/sub=" << (sub / sub2) << endl;

   // Strided = reference to a strided subset
   SxVecRef<double, Strided> diag = mat.diag ();

   // vector-vector arithmetic with strided & contiguous operands
   cout << "diag=" << diag << endl;
   cout << "a+diag" << (a + diag) << endl;
   cout << "diag+a" << (diag + a) << endl;
   cout << "a-diag" << (a - diag) << endl;
   cout << "diag-a" << (diag - a) << endl;
   cout << "a*diag" << (a * diag) << endl;
   cout << "diag*a" << (diag * a) << endl;
   cout << "a/diag" << (a / diag) << endl;
   cout << "diag/b" << (diag / b) << endl;

   // references to columns (contiguous)
   SX_LOOP(c)  {
      cout << "mat column " << (c+1) << ": " << mat.colRef (c) << endl;
      cout << "sub column " << (c+1) << ": " << sub.colRef (c) << endl;
   }

   // references to rows (strided)
   SX_LOOP(r)
      cout << "mat row " << (r+1) << ": " << mat.rowRef (r) << endl;
   SX_LOOP(r)
      cout << "sub row " << (r+1) << ": " << sub.rowRef (r) << endl;

   // --- swap rows 3 (index 2) and 4 (index 3)
   {
      // copy
      SxVector<double> oldRow2 = mat.rowRef (2);
      // copy 3->2 (in place)
      mat.rowRef(2) = mat.rowRef(3);
      // copy ->3 (in place)
      mat.rowRef(3) = oldRow2;
   }
   cout << mat << endl;

   // GeneralMatrix = reference with stride for rows, and possibly gaps
   //                 between columns
   SxVecRef<double, GeneralMatrix> genMat = mat.getRef<GeneralMatrix> (1, 5, 10, 2);
   cout << "General matrix: every second row" << endl;
   cout << "gen=" << genMat << endl;

   // vector arithmetic: operands are SubMatrix and GeneralMatrix
   cout << "gen + sub = " << (genMat + sub) << endl;
   cout << "gen - sub = " << (genMat - sub) << endl;
   cout << "gen * sub = " << (genMat * sub) << endl;
   cout << "gen / sub = " << (genMat / sub) << endl;

   cout << "sub + gen = " << (sub + genMat) << endl;
   cout << "sub - gen = " << (sub - genMat) << endl;
   cout << "sub * gen = " << (sub * genMat) << endl;
   cout << "sub / gen = " << (sub / genMat) << endl;

   // in-place arithmetic: contiguous vector
   SxVector<double> c = b;
   cout << "b=" << b << endl << "c=" << c << endl;
   c += b;
   cout << "b=" << b << endl << "c=" << c << endl;
   c = b;
   c -= b;
   cout << "b=" << b << endl << "c=" << c << endl;
   c = b;
   c *= b;
   cout << "b=" << b << endl << "c=" << c << endl;
   c = b;
   c /= b;
   cout << "b=" << b << endl << "c=" << c << endl;

   // in-place arithmetic: strided vector
   SxVecRef<double,Strided> d = mat.rowRef (2);
   cout << "d=" << d << endl;
   c = b;
   c += d;
   cout << "c=" << c << endl;
   c = b;
   c -= d;
   cout << "c=" << c << endl;
   c = b;
   c *= d;
   cout << "c=" << c << endl;
   c = b;
   c /= d;
   cout << "c=" << c << endl;

   // --- scalar-vector operations
   cout << (a + 5.) << endl;
   cout << (a - 5.) << endl;
   cout << (a * 5.) << endl;
   cout << (a / 5.) << endl;

   cout << (5. + a) << endl;
   cout << (5. - a) << endl;
   cout << (5. * a) << endl;
   cout << (5. / b) << endl;

   cout << ((a + a) + 5.) << endl;


   // in-place arithmetic with scalar operand
   c = b; c+= 5.; cout << c << endl;
   c = b; c-= 5.; cout << c << endl;
   c = b; c*= 5.; cout << c << endl;
   c = b; c/= 5.; cout << c << endl;

   // in-place arithmetic with scalar operand
   // and strided vector
   cout << "rowRef and scalar" << endl;
   cout << mat.rowRef(2) << endl;
   cout << (mat.rowRef (2) + 5.) << endl;
   cout << (mat.rowRef (2) - 5.) << endl;
   cout << (mat.rowRef (2) * 5.) << endl;
   cout << (mat.rowRef (2) / 5.) << endl;
   cout << (5. + mat.rowRef (2)) << endl;
   cout << (5. - mat.rowRef (2)) << endl;
   cout << (5. * mat.rowRef (2)) << endl;
   cout << (5. / mat.rowRef (2)) << endl;

   SxVector<double> mat3 = mat;

   // in-place arithmetic: type mixing (here: int and double)
   cout << mat3 << endl;
   mat3.rowRef (2) += 20;
   cout << mat3 << endl;

   // --- integer type
   SxVector<int> iVec(10);
   SX_LOOP(i) iVec(i) = i + 1;
   cout << (iVec / 2) << endl;

   cout << (a + iVec / 2) << endl;

   // computations involving temporaries may switch to in-place
   // algorithm, reusing the memory of the about-to-die vector
   cout << "Compute on temp:" << endl;
   cout << ((a + a + 2. * a) + (3. * a)) << endl;

   // --- functions
   cout << "a.sqr=" << a.sqr () << endl;
   cout << "a.cub=" << a.cub () << endl;
   cout << "a.abs=" << a.abs () << endl;
   cout << "a.absSqr=" << a.absSqr () << endl;
   cout << "exp(a)=" << exp(a) << endl;
   cout << "erf(a)=" << erf(a) << endl;
   cout << "f(a)=" << (a - a.sqr () + 3.) << endl;
   cout << "f(a)=" << (a - a.sqr () + 3.).abs () << endl;

   SxVecRef<double,Strided> strided = mat.rowRef (2);
   cout << "strided=" << strided << endl;
   cout << "strided.sqr=" << strided.sqr () << endl;
   cout << "strided.cub=" << strided.cub () << endl;
   cout << "strided.abs=" << strided.abs () << endl;
   cout << "strided.absSqr=" << strided.absSqr () << endl;
   cout << "exp(strided)=" << exp(strided) << endl;

   // reduction (sum, product, normSqr) for Simd
   {
      cout << "--- Simd reductions" << endl;
      double vSum = 0.;
      SxVector<double> &v = mat;
      SX_LOOP(i) vSum += v(i);
      cout << "vSum = " << vSum << endl;
      cout << "v.sum () = " << v.sum () << endl;
      double vProd = 1.;
      SX_LOOP(i) vProd *= v(i);
      cout << "vProd = " << vProd << endl;
      cout << "v.product () = " << v.product () << endl;
      double v2Sum = 0.;
      SX_LOOP(i) v2Sum += v(i) * v(i);
      cout << "v2Sum = " << v2Sum << endl;
      cout << "v.normSqr () = " << v.normSqr () << endl;
      cout << "v.absSqr ().sum () = " << v.absSqr ().sum () << endl;
   }
   // reduction (sum, product, normSqr) for SimdCol
   {
      cout << "--- SimdCol reductions" << endl;
      double vSum = 0.;
      SxVecRef<double,SubMatrix> &v = sub;
      SX_LOOP(i) vSum += v(i);
      cout << "vSum = " << vSum << endl;
      cout << "v.sum () = " << v.sum () << endl;
      double vProd = 1.;
      SX_LOOP(i) vProd *= v(i);
      cout << "vProd = " << vProd << endl;
      cout << "v.product () = " << v.product () << endl;
      double v2Sum = 0.;
      SX_LOOP(i) v2Sum += v(i) * v(i);
      cout << "v2Sum = " << v2Sum << endl;
      cout << "v.normSqr () = " << v.normSqr () << endl;
      cout << "v.absSqr ().sum () = " << v.absSqr ().sum () << endl;
   }
   // reduction (sum, product, normSqr) for UseIterator
   {
      cout << "--- UseIterator reductions" << endl;
      double vSum = 0.;
      SxVecRef<double,Strided> &v = strided;
      SX_LOOP(i) vSum += v(i);
      cout << "vSum = " << vSum << endl;
      cout << "v.sum () = " << v.sum () << endl;
      double vProd = 1.;
      SX_LOOP(i) vProd *= v(i);
      cout << "vProd = " << vProd << endl;
      cout << "v.product () = " << v.product () << endl;
      double v2Sum = 0.;
      SX_LOOP(i) v2Sum += v(i) * v(i);
      cout << "v2Sum = " << v2Sum << endl;
      cout << "v.normSqr () = " << v.normSqr () << endl;
      cout << "v.absSqr ().sum () = " << v.absSqr ().sum () << endl;
   }
   cout << "---" << endl;

   // --- matrix access with (i,j)
   SX_LOOP2(j,i)
      cout << "mat2(" << i << "," << j << ")=" << mat2(i,j) << endl;
   SX_LOOP2(j,i)
      cout << "sub(" << i << "," << j << ")=" << sub(i,j) << endl;
   SX_LOOP2(j,i)
      cout << "genMat(" << i << "," << j << ")=" << genMat(i,j) << endl;

   cout << "---" << endl;

   // minimum and maximum values
   cout << a.minval () << ' ' << a.maxval () << endl;
   cout << mat.minval () << ' ' << mat.maxval () << endl;

   // sorting
   cout << a.getSorted ((a * -1.).getSortIdx ()) << endl;
   cout << a.integrate (0.01) << endl;

   // test some functions specific for complex numbers
   SxVector<SxComplex16> complexVec(10);
   SX_LOOP(i) complexVec(i) = SxComplex16::phase (i);
   complexVec.reshape (2,5);
   cout << complexVec << endl;
   cout << complexVec.real () << endl;
   cout << complexVec.imag () << endl;
   cout << complexVec.conj () << endl;
   // this is a strided reference!
   cout << complexVec.rowRef (1).real () << endl;
   // this is a strided reference!
   cout << complexVec.rowRef (1).imag () << endl;
   // this is (presently) a new vector
   cout << complexVec.rowRef (1).conj () << endl;

   // test transpose for complex vectors (rare use)
   cout << complexVec.transpose () << endl;

   // --- test transpose for large matrices (uses cache-oblivious algorithm)
   SxVector<double> mat4;
   mat4.reformat (5000,5000);
   SX_LOOP2(i,j)
      mat4(i,j) = i + 1e-4 * j;

   SxVector<double> mat4T = mat4.transpose ();
   for (int i = 0; i < 10; ++i)
      for (int j = 0; j < 10; ++j)
         cout << mat4(i,j) << " " << mat4T(i,j) << endl;

   {
      double res = 0.;
      SX_LOOP2(i,j)
         res += sqr (mat4(i,j) - mat4T(j,i));
      cout << "transposeCheck: " << res << endl;
   }

   // --- test transpose for non-contiguous matrices
   cout << "Submatrix transpose:" << endl;
   cout << sub << endl;
   cout << sub.transpose () << endl;

   cout << "General matrix transpose:" << endl;
   cout << genMat << endl;
   cout << genMat.transpose () << endl;

   cout << "Strided matrix transpose:" << endl;
   cout << complexVec.real () << endl;
   cout << complexVec.real ().transpose () << endl;

   // test adjoint
   cout << complexVec.transpose () << endl;
   cout << complexVec.adjoint () << endl;
   cout << (complexVec.adjoint () - complexVec.transpose ().conj ()).normSqr () << endl;

   cout << "---" << endl;

   initSPHInXMath ();
   {
      SxVector<double> x = a;
      x.normalize ();
      cout << x << endl << x.norm () << endl;
      x.randomize ();
      cout << x << endl << x.norm () << endl;
   }

   // SxIdx interface
   cout << a(SxIdx(1,4)) << endl;
   cout << mat(SxIdx(1,20)) << endl;

   // axpy interface
   {
      SxVector<double> x = a;
      x.plus_assign_ax (0.74, b);
      cout << x << endl;
      cout << (x - (a + 0.74 * b)).norm () << endl;
   }

   cout << "--- scalar product" << endl;
   cout << dot(a,b) << endl;
   cout << (a * b).sum () << endl;


   cout << "--- matrix multiply" << endl;

   // matrix-vector
   cout << (mat ^ a) << endl;
   {
      SxVector<double> x(10);
      SX_LOOP(i)
         x(i) = (mat.rowRef (i) * a).sum ();
      cout << x << endl;
   }

   // matrix-matrix
   cout << (mat ^ mat) << endl;
   {
      SxVector<double> x;
      x.reformat (10, 10);
      SX_LOOP2(i,j)
         x(i,j) = (mat.rowRef (i) * mat.colRef (j)).sum ();
      cout << (x - (mat ^ mat)).norm () << endl;
   }

   // matrix-matrix (first operand is submatrix)
   cout << (sub ^ mat) << endl;
   {
      SxVector<double> x;
      x.reformat (5, 10);
      SX_LOOP2(i,j)
         x(i,j) = (sub.rowRef (i) * mat.colRef (j)).sum ();
      cout << (x - (sub ^ mat)).norm () << endl;
   }

   // matrix-matrix (2nd operand is submatrix)
   SxVector<double> mat10x5 = sub.transpose ();
   cout << (mat10x5 ^ sub) << endl;
   {
      SxVector<double> x;
      x.reformat (10, 10);
      SX_LOOP2(i,j)
         x(i,j) = (sub.colRef (i) * sub.colRef (j)).sum ();
      cout << (x - (sub.transpose () ^ sub)).norm () << endl;
   }

   // matrix-matrix (2nd operand is general matrix)
   cout << (mat10x5 ^ genMat) << endl;
   {
      SxVector<double> x;
      x.reformat (10, 10);
      SX_LOOP2(i,j)
         x(i,j) = (mat10x5.rowRef (i) * genMat.colRef (j)).sum ();
      cout << (x - (mat10x5 ^ genMat)).norm () << endl;
   }

   // matrix-matrix (operands are general matrix/SubMatrix)
   {
      SxVector<double> mat5;
      mat5.reformat (M,M);
      SX_LOOP2(i,j)
         mat5(i,j) = i + 1e-4 * j;
      SxVecRef<double,GeneralMatrix> genMat5
         = mat5.getRef<GeneralMatrix>(1, M/2, M, 2);
      SxVecRef<double,SubMatrix> blockMat5
         = mat5.blockRef(0,0,M,M/2);
      SxVector<double> x;
      x.reformat (blockMat5.getNRows (), genMat5.getNCols ());
      SX_LOOP2(i,j)
         x(i,j) = dot(blockMat5.rowRef (i), genMat5.colRef (j));
      cout << (x - (blockMat5 ^ genMat5)).norm () / x.norm () << endl;
      x.reformat (M/2, M/2);
      SX_LOOP2(i,j)
         x(i,j) = dot(genMat5.rowRef (i), blockMat5.colRef (j));
      cout << (x - (genMat5 ^ blockMat5)).norm () / x.norm () << endl;
   }
   /*
   cout << (genMat ^ mat) << endl;
   {
      SxVector<double> x;
      x.reformat (5, 10);
      SX_LOOP2(i,j)
         x(i,j) = (genMat.rowRef (i) * mat.colRef (j)).sum ();
      cout << (x - (genMat ^ mat)).norm () << endl;
   }
   */

   cout << "--- overlap matrices" << endl;
   cout << (mat.transpose () ^ mat) << endl;
   cout << mat.overlap (mat) << endl;
   cout << (mat.overlap (mat) - (mat.transpose () ^ mat)).norm () << endl;
   cout << (complexVec.overlap (complexVec) - (complexVec.adjoint () ^ complexVec)).norm () << endl;
   cout << (sub.overlap (mat2) - (sub.transpose () ^ mat2)).norm () << endl;

   cout << "--- in place rotation" << endl;
   {
      SxVector<double> mat6 = mat - 1.;
      mat6.rotate (mat);
      cout << (mat6 - ((mat - 1.) ^ mat)) << endl;
   }

   cout << "--- inverse" << endl;
   {
      SxVector<double> matI = mat;
      matI.diag () += 1.;
      SxVector<double> one;
      one.reformat (matI.getNRows (), matI.getNCols ());
      one.set (0.);
      one.diag () = 1.;
      cout << matI.inverse () << endl;
      cout << ((matI.inverse () ^ matI) - one).norm () << endl;
      cout << ((matI + 0.).inverse () ^ (matI)).diag () << endl;
   }

   cout << "--- solve" << endl;
   {
      SxVector<double> x = mat.solve (a);
      SxVector<double> R = (mat ^ x) - a;
      cout << "R= " << R << endl;
      cout << (R.reshape (1,10) ^ mat).norm () << endl;
      // temporary versions
      cout << (x - mat.solve (a + 0.)).norm () << endl;
      cout << (x - (mat + 0.).solve (a)).norm () << endl;
      cout << (x - (mat + 0.).solve (a + 0.)).norm () << endl;
   }

   cout << "--- pseudo-inverse" << endl;
   {
      SxVector<double> matI = mat.pseudoInverse ();
      SxVector<double> x = matI ^ a;
      SxVector<double> R = (mat ^ x) - a;
      cout << "R= " << R << endl;
      cout << (R.reshape (1,10) ^ mat).norm () << endl;

      SxVector<double> mat7 = mat2;
      for (int i = 0; i < mat7.getNRows (); i++)
         mat7(i,i)+=1.;
      cout << symEigenvalues (mat7 ^ mat7.transpose ()) << endl;
      matI = mat7.pseudoInverse ();
      cout << (mat7 ^ matI) << endl; // should be unit matrix
      cout << symEigenvalues (matI ^ mat7) << endl; // should be 5x zero and 5x one
   }

   cout << "--- cholesky" << endl;
   {
      SxVector<double> sym = mat.overlap (mat);
      sym.diag () += 1.;
      SxVector<double> L = sym.choleskyDecomposition ();
      cout << ((L ^ L.transpose ()) - sym).norm () / sym.norm () << endl;
      SxVector<double> R = sym.choleskyDecomposition (UpperRight);
      cout << ((R.transpose () ^ R) - sym).norm () / sym.norm () << endl;
   }

   cout << "--- eigensystem " << endl;
   {
      SxVector<double> sym = mat.overlap (mat);
      SxVector<SxComplex16> symC = sym;
      cout << sym << endl;
      SxSymEigensystem<double> eig(sym);
      cout << eig.vals << endl;
      SxVector<double> D;
      D.reformat (sym.getNCols (), sym.getNCols ());
      D.set (0.);
      D.diag () = eig.vals;
      cout << eig.vecs.overlap (eig.vecs) << endl;
      cout << ((eig.vecs ^ D ^ eig.vecs.transpose ()) - sym).norm () / sym.norm () << endl;
      sym.diag () += 1.;
      // temporary version
      cout << symEigenvalues (sym.inverse ()) << endl;
      // demonstrate that temporary version overwrites original memory
      cout << SxSymEigensystem<double> (std::move (sym)).vals << endl;
      cout << sym << endl;

      // Hermitean matrices
      cout << (symC - symC.adjoint ()).norm () / symC.norm () << endl;
      SxSymEigensystem<SxComplex16> eigC(symC);
      cout << eigC.vals << endl;
   }

   {
      cout << "Symmetric packed" << endl;
      SxVector<double> sym = mat.overlap (mat);
      sym.diag () += 0.1; // make sure sym is invertible

      SxSymMatrix<double> pSym(mat.getNRows ());
      for (ssize_t j = 0; j < pSym.getNCols (); j++)  {
         for (ssize_t i = 0; i <= j; i++)  {
            pSym(i,j) = sym(i,j);
         }
      }
      cout << "Expand: " << endl;
      cout << pSym << endl;
      cout << pSym.expand () << endl;
      cout << (pSym.expand () - sym).normSqr () << endl;
      cout << "Inverse: " << endl;
      SxSymMatrix<double> pSymInv = pSym.inverse ();
      cout << (pSym.expand () ^ pSymInv.expand ()) << endl;

      cout << (pSymInv.expand () - sym.inverse ()).normSqr () << endl;
      cout << "matmul test: " << ((pSym ^ mat) - (pSym.expand () ^ mat)).normSqr () << endl;

      pSym.resize (5000);
      for (ssize_t j = 0; j < pSym.getNCols (); j++)  {
         for (ssize_t i = 0; i <= j; i++)  {
            pSym(i,j) = i + 1e-4 * j;
         }
      }
      cout << "Expand (5000): " << endl;
      SxVector<double> expanded = pSym.expand ();
      cout << (expanded - expanded.transpose ()).normSqr () << endl;
      for (ssize_t j = 0; j < pSym.getNCols (); j++)  {
         for (ssize_t i = 0; i <= j; i++)  {
            SX_CHECK (pSym(i,j) == expanded(i,j), i, j);
         }
      }
   }

   {
      cout << "Diagonal matrix" << endl;
      SxDiagMat<double> D(a);
      SxVector<double> D2(10,10);
      D2.set (0.); D2.diag () = a;
      cout << (D ^ b) << '=' << (a * b) << endl;
      cout << (D ^ b) << '=' << (a * b) << endl;
      cout << "D^mat=" << (D ^ mat) << endl;
      cout << "check vs expanded diag: " << ((D ^ mat) - (D2 ^ mat)).norm () << endl;
      cout << "mat^D=" << (mat ^ D) << endl;
      cout << "check vs expanded diag: " << ((mat ^ D) - (mat ^ D2)).norm () << endl;
      cout << "strided: " << (D ^ mat.rowRef (2)) << endl;
      cout << "subMat(1..2,*) ^ D: " << (mat.blockRef (1,0,2,10) ^ D) << endl;
      cout << "D ^ subMat(*,1..2): " << (D ^ mat.blockRef (0,1,10,2)) << endl;
   }


   destroySFHIngXMath ();
}
