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
#ifndef _SX_BLASLIB_H_
#define _SX_BLASLIB_H_

#include <SxComplex.h>
#include <SxUtil.h>
#include <SxConfig.h>
#include <SxMath.h>

enum UPLO   { UpperRight, LowerLeft };
enum EIGCMD { All, ValuesOnly, VectorsOnly, OptSize };

//------------------------------------------------------------------------------
// norm of vectors
//------------------------------------------------------------------------------
// --- The next 5 lines were generated from snippets/SxBlasLib.h snippet norm2
SX_EXPORT_MATH float  norm2 (const float *vec, ssize_t n, int incx = 1);
SX_EXPORT_MATH double norm2 (const double *vec, ssize_t n, int incx = 1);
SX_EXPORT_MATH float  norm2 (const SxComplex8 *vec, ssize_t n, int incx = 1);
SX_EXPORT_MATH double norm2 (const SxComplex16 *vec, ssize_t n, int incx = 1);
// --- norm2

//------------------------------------------------------------------------------
// scale vectors
//------------------------------------------------------------------------------
SX_EXPORT_MATH void scale (float *vec,  const float alpha, int n);
SX_EXPORT_MATH void scale (double *vec, const double alpha, int n);
SX_EXPORT_MATH void scale (SxComplex8 *vec,  const SxComplex8 &alpha, int n);
SX_EXPORT_MATH void scale (SxComplex16 *vec, const SxComplex16 &alpha, int n);


//------------------------------------------------------------------------------
// Y += a*X
//------------------------------------------------------------------------------
// --- The next 24 lines were generated from ./snippets/SxBlasLib.h snippet axpy
SX_EXPORT_MATH void axpy (float *yOut, int incy, const float &alpha, const float *xIn, int incx, ssize_t n);
inline void axpy (float *yOut, const float &alpha, const float *xIn, ssize_t n)
{
   axpy (yOut, 1, alpha, xIn, 1, n);
}

SX_EXPORT_MATH void axpy (double *yOut, int incy, const double &alpha, const double *xIn, int incx, ssize_t n);
inline void axpy (double *yOut, const double &alpha, const double *xIn, ssize_t n)
{
   axpy (yOut, 1, alpha, xIn, 1, n);
}

SX_EXPORT_MATH void axpy (SxComplex8 *yOut, int incy, const SxComplex8 &alpha, const SxComplex8 *xIn, int incx, ssize_t n);
inline void axpy (SxComplex8 *yOut, const SxComplex8 &alpha, const SxComplex8 *xIn, ssize_t n)
{
   axpy (yOut, 1, alpha, xIn, 1, n);
}

SX_EXPORT_MATH void axpy (SxComplex16 *yOut, int incy, const SxComplex16 &alpha, const SxComplex16 *xIn, int incx, ssize_t n);
inline void axpy (SxComplex16 *yOut, const SxComplex16 &alpha, const SxComplex16 *xIn, ssize_t n)
{
   axpy (yOut, 1, alpha, xIn, 1, n);
}
// --- axpy


//------------------------------------------------------------------------------
// scalar product of vectors
//------------------------------------------------------------------------------
SX_EXPORT_MATH float  scalarProduct (const float *aVec, const float *bVec, 
                                     int n);
SX_EXPORT_MATH double scalarProduct (const double *aVec, const double *bVec, 
                                     int n);
SX_EXPORT_MATH SxComplex8 scalarProduct (const SxComplex8 *aVec, 
                                         const SxComplex8 *bVec, int n);
SX_EXPORT_MATH SxComplex16 scalarProduct (const SxComplex16 *aVec, 
                                          const SxComplex16 *bVec, int n);

//------------------------------------------------------------------------------
// general matrix-matrix multiplication
// and symmetric packed matrix-vector multiplication
//------------------------------------------------------------------------------
// --- The next 92 lines were generated from snippets/SxBlasLib.h snippet matmult
SX_EXPORT_MATH void matmult (float *resMat, const float beta, const float *aMat,
                             const float *bMat, int aMatRows, int aMatCols,
                             int bMatCols, int lda, int ldb, int ldc);
// wrapper for compact matrices
inline void matmult (float *resMat, const float *aMat,
                     const float *bMat, int aMatRows, int aMatCols,
                     int bMatCols)
{
   matmult(resMat, 0., aMat, bMat,
           aMatRows, aMatCols, bMatCols,
           aMatRows, aMatCols, aMatRows);
}
/** \brief multiply packed symmetric matrix with a vector
  @param n number of elements in vector
  @param resVec target vector. Target vector is always stride 1.
  @param aPacked packed symmetric/Hermitean n x n matrix
  @param uplo   packing style
  @param xVec  vector to be multiplied with matrix
  @param xStride element stride of vector x
*/
SX_EXPORT_MATH void symmatmult (int n, float *resVec, const float *aPacked,
                                enum UPLO uplo, const float *xVec, int xStride);

SX_EXPORT_MATH void matmult (double *resMat, const double beta, const double *aMat,
                             const double *bMat, int aMatRows, int aMatCols,
                             int bMatCols, int lda, int ldb, int ldc);
// wrapper for compact matrices
inline void matmult (double *resMat, const double *aMat,
                     const double *bMat, int aMatRows, int aMatCols,
                     int bMatCols)
{
   matmult(resMat, 0., aMat, bMat,
           aMatRows, aMatCols, bMatCols,
           aMatRows, aMatCols, aMatRows);
}
/** \brief multiply packed symmetric matrix with a vector
  @param n number of elements in vector
  @param resVec target vector. Target vector is always stride 1.
  @param aPacked packed symmetric/Hermitean n x n matrix
  @param uplo   packing style
  @param xVec  vector to be multiplied with matrix
  @param xStride element stride of vector x
*/
SX_EXPORT_MATH void symmatmult (int n, double *resVec, const double *aPacked,
                                enum UPLO uplo, const double *xVec, int xStride);

SX_EXPORT_MATH void matmult (SxComplex8 *resMat, const SxComplex8 &beta, const SxComplex8 *aMat,
                             const SxComplex8 *bMat, int aMatRows, int aMatCols,
                             int bMatCols, int lda, int ldb, int ldc);
// wrapper for compact matrices
inline void matmult (SxComplex8 *resMat, const SxComplex8 *aMat,
                     const SxComplex8 *bMat, int aMatRows, int aMatCols,
                     int bMatCols)
{
   matmult(resMat, 0., aMat, bMat,
           aMatRows, aMatCols, bMatCols,
           aMatRows, aMatCols, aMatRows);
}
/** \brief multiply packed symmetric matrix with a vector
  @param n number of elements in vector
  @param resVec target vector. Target vector is always stride 1.
  @param aPacked packed symmetric/Hermitean n x n matrix
  @param uplo   packing style
  @param xVec  vector to be multiplied with matrix
  @param xStride element stride of vector x
*/
SX_EXPORT_MATH void symmatmult (int n, SxComplex8 *resVec, const SxComplex8 *aPacked,
                                enum UPLO uplo, const SxComplex8 *xVec, int xStride);

SX_EXPORT_MATH void matmult (SxComplex16 *resMat, const SxComplex16 &beta, const SxComplex16 *aMat,
                             const SxComplex16 *bMat, int aMatRows, int aMatCols,
                             int bMatCols, int lda, int ldb, int ldc);
// wrapper for compact matrices
inline void matmult (SxComplex16 *resMat, const SxComplex16 *aMat,
                     const SxComplex16 *bMat, int aMatRows, int aMatCols,
                     int bMatCols)
{
   matmult(resMat, 0., aMat, bMat,
           aMatRows, aMatCols, bMatCols,
           aMatRows, aMatCols, aMatRows);
}
/** \brief multiply packed symmetric matrix with a vector
  @param n number of elements in vector
  @param resVec target vector. Target vector is always stride 1.
  @param aPacked packed symmetric/Hermitean n x n matrix
  @param uplo   packing style
  @param xVec  vector to be multiplied with matrix
  @param xStride element stride of vector x
*/
SX_EXPORT_MATH void symmatmult (int n, SxComplex16 *resVec, const SxComplex16 *aPacked,
                                enum UPLO uplo, const SxComplex16 *xVec, int xStride);
// --- matmult

//------------------------------------------------------------------------------
// overlap matrices
//------------------------------------------------------------------------------
// computes top(A)^T* top(B)
// top(X) restricts matrix to sumSize many rows (from top)
SX_EXPORT_MATH void matovlp (float *resMat, 
                             const float *aMat, const float *bMat, 
                             int aMatRows, int aMatCols, int bMatRows,
                             int bMatCols, int sumSize);
SX_EXPORT_MATH void matovlp (double *resMat, 
                             const double *aMat, const double *bMat, 
                             int aMatRows, int aMatCols, int bMatRows,
                             int bMatCols, int sumSize);
SX_EXPORT_MATH void matovlp (SxComplex8 *resMat, 
                             const SxComplex8 *aMat, const SxComplex8 *bMat, 
                             int aMatRows, int aMatCols, int bMatRows,
                             int bMatCols, int sumSize);
SX_EXPORT_MATH void matovlp (SxComplex16 *resMat, 
                             const SxComplex16 *aMat, const SxComplex16 *bMat, 
                             int aMatRows, int aMatCols, int bMatRows,
                             int bMatCols, int sumSize);
//------------------------------------------------------------------------------
// Matrix decompositions
//------------------------------------------------------------------------------
SX_EXPORT_MATH void cholesky (float *resMat, enum UPLO,   // modifies inMat!!!
                              float *inMat, int n);
SX_EXPORT_MATH void cholesky (double *resMat, enum UPLO, 
                              double *inMat, int n);
SX_EXPORT_MATH void cholesky (SxComplex8 *resMat, enum UPLO, 
                              SxComplex8 *inMat, int n);
SX_EXPORT_MATH void cholesky (SxComplex16 *resMat, enum UPLO, 
                              SxComplex16 *inMat, int n);

SX_EXPORT_MATH void singularValueDecomp (float *mat, int nRows, int nCols,
                                         float *vals,
                                         float *left,
                                         float *right, // V^H
                                         bool zeroSpace);
SX_EXPORT_MATH void singularValueDecomp (double *mat, int nRows, int nCols,
                                         double *vals,
                                         double *left,
                                         double *right, // V^H
                                         bool zeroSpace);
SX_EXPORT_MATH void singularValueDecomp (SxComplex8 *mat, int nRows, int nCols,
                                         float *vals,
                                         SxComplex8 *left,
                                         SxComplex8 *right, // V^H
                                         bool zeroSpace);
SX_EXPORT_MATH void singularValueDecomp (SxComplex16 *mat, int nRows, int nCols,
                                         double *vals,
                                         SxComplex16 *left,
                                         SxComplex16 *right, // V^H
                                         bool zeroSpace);

//------------------------------------------------------------------------------
// matrix inversion
//------------------------------------------------------------------------------
SX_EXPORT_MATH void matInverse (float       *mat, int nRows, int nCols);
SX_EXPORT_MATH void matInverse (double      *mat, int nRows, int nCols);
SX_EXPORT_MATH void matInverse (SxComplex8  *mat, int nRows, int nCols);
SX_EXPORT_MATH void matInverse (SxComplex16 *mat, int nRows, int nCols);
SX_EXPORT_MATH void matInverseTri (float       *mat, int nRows, enum UPLO);
SX_EXPORT_MATH void matInverseTri (double      *mat, int nRows, enum UPLO);
SX_EXPORT_MATH void matInverseTri (SxComplex8  *mat, int nRows, enum UPLO);
SX_EXPORT_MATH void matInverseTri (SxComplex16 *mat, int nRows, enum UPLO);

//------------------------------------------------------------------------------
// Linear equation solver (least sqare based)
//------------------------------------------------------------------------------
SX_EXPORT_MATH void solveLinEq (float *mat, int nRows, int nCols, 
                                float *b,   int bCols);
SX_EXPORT_MATH void solveLinEq (double *mat, int nRows, int nCols, 
                                double *b,   int bCols);
SX_EXPORT_MATH void solveLinEq (SxComplex8 *mat, int nRows, int nCols, 
                                SxComplex8 *b,   int bCols);
SX_EXPORT_MATH void solveLinEq (SxComplex16 *mat, int nRows, int nCols, 
                                SxComplex16 *b,   int bCols);

//------------------------------------------------------------------------------
// eigensolver
//------------------------------------------------------------------------------
// modifies inMat!!!
// --- The next 25 lines were generated from snippets/SxBlasLib.h snippet eigensolver
SX_EXPORT_MATH int matEigensolver (SxComplex8 *eigVals, float *eigVecs,
                                   float *inMat, int n, EIGCMD cmd=All,
                                   int size=0);
SX_EXPORT_MATH void matSymEigensolver (float  *eigVals, float *eigVecs,
                                       float *inMat, int n,
                                       EIGCMD cmd=All);
SX_EXPORT_MATH int matEigensolver (SxComplex16 *eigVals, double *eigVecs,
                                   double *inMat, int n, EIGCMD cmd=All,
                                   int size=0);
SX_EXPORT_MATH void matSymEigensolver (double *eigVals, double *eigVecs,
                                       double *inMat, int n,
                                       EIGCMD cmd=All);
SX_EXPORT_MATH int matEigensolver (SxComplex8 *eigVals, SxComplex8 *eigVecs,
                                   SxComplex8 *inMat, int n, EIGCMD cmd=All,
                                   int size=0);
SX_EXPORT_MATH void matSymEigensolver (float  *eigVals, SxComplex8 *eigVecs,
                                       SxComplex8 *inMat, int n,
                                       EIGCMD cmd=All);
SX_EXPORT_MATH int matEigensolver (SxComplex16 *eigVals, SxComplex16 *eigVecs,
                                   SxComplex16 *inMat, int n, EIGCMD cmd=All,
                                   int size=0);
SX_EXPORT_MATH void matSymEigensolver (double *eigVals, SxComplex16 *eigVecs,
                                       SxComplex16 *inMat, int n,
                                       EIGCMD cmd=All);
// --- eigensolver
// modified inMat!!!
SX_EXPORT_MATH void matEigensolverTri (float *eigVals, float *eigVecs,
                                       float *inMat,
                                       int n, enum UPLO, EIGCMD cmd=All);
SX_EXPORT_MATH void matEigensolverTri (double *eigVals, double *eigVecs,
                                       double *inMat,
                                       int n, enum UPLO, EIGCMD cmd=All);
SX_EXPORT_MATH void matEigensolverTri (float *eigVals, SxComplex8 *eigVecs,
                                       SxComplex8 *inMat,
                                       int n, enum UPLO, EIGCMD cmd=All);
SX_EXPORT_MATH void matEigensolverTri (double *eigVals, SxComplex16 *eigVecs,
                                       SxComplex16 *inMat,
                                       int n, enum UPLO, EIGCMD cmd=All);


//------------------------------------------------------------------------------
// Some missing pototypes
//------------------------------------------------------------------------------
extern "C" {
//   zgetri (int, void *, int, int *, void *, int, int &);
}


#endif /* _SX_BLASLIB_H_ */
