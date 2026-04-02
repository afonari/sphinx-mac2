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
// ref1 - Comp. Phys. Comm (128), 1-45 (2000)

//#include <stdio.h>
//#include <string.h>
//#include <stdlib.h>
//#include <iostream>
#include <math.h>
#include <SxBlasLib.h>
#include <SxError.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif

//#include <SxRandom.h>



// --- BLAS, LAPACK
#if defined(USE_VECLIB)
//#   include <veclib.h>   // --- don't use HP header files!
//#   include <lapack.h>   //     more details in SxMLIB.h
#     include <SxMLIB.h>
#elif defined (ESSE_ESSL)
    // --- before including ESSL we have to do a messy definition
    //     otherwise it doesn't compile. i don't know a better way.
    //     see also: SxFFT.h
#   define _ESV_COMPLEX_
#   include <complex>
#   include <essl.h>
#elif defined (USE_INTEL_MKL)
// mkl uses 'long long', which is not ISO C++. Disable errors here.
//#pragma GCC diagnostic push
//#pragma GCC diagnostic warning "-Wlong-long"
    extern "C"  {
#      include <mkl.h>
    }
//#pragma GCC diagnostic pop
#elif defined (USE_ACML)
    extern "C"  {
#      include <acml.h>
    }
#elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
#  define blasint int
       extern "C"  {
#         include <f2c.h>
#         include <cblas.h>
#         include <clapack.h>
       }
#elif defined (USE_ATLAS)
#  if defined (USE_ACCELERATE_FRAMEWORK)
#     include <Accelerate/Accelerate.h>
      typedef __CLPK_integer integer;
      typedef __CLPK_logical logical;
      typedef __CLPK_real real;
      typedef __CLPK_doublereal doublereal;
      typedef __CLPK_ftnlen ftnlen;
      typedef __CLPK_complex complex;
      typedef __CLPK_doublecomplex doublecomplex;
#  else
       extern "C"  { 
#         include <f2c.h>
#         include <cblas.h>
#         include <clapack.h>
       }
#  endif /* USE_ACCELERATE_FRAMEWORK */
#elif defined (USE_NETLIB)
// define lapacke_ types
#  define lapack_complex_float  SxComplex8
#  define lapack_complex_double SxComplex16
   extern "C"  {
#         include <cblas.h>
#         include <lapacke.h>
   }
#else
#   error "No numeric library specified"
// make sure that we do not drown in error messages
#define SX_IGNORE_THE_REST_OF_THE_FILE
#endif

#ifndef SX_IGNORE_THE_REST_OF_THE_FILE
//------------------------------------------------------------------------------
// BLAS/LAPACK error handling   
//------------------------------------------------------------------------------
#if   defined (USE_VECLIB)
   // error handling not yet supported
#elif defined (USE_ESSL)
   // error handling not yet supported
#elif defined (USE_INTEL_MKL)
   // error handling not yet supported
#elif defined (USE_ACML)
   // error handling not yet supported
#elif defined (USE_GOTO)
   // error handling not yet supported
#else
#  include <stdarg.h>

#ifdef MACOSX
#  if ( DIST_VERSION_L >= 1070L )
      extern "C" void cblas_xerbla (int p, char *rout, char *form, ...)
#  else
      extern "C" void cblas_xerbla (int p, const char *rout, const char *form, ...)
#  endif /* DIST_VERSION_L */
#else
   extern "C" void cblas_xerbla (int p, const char *rout, const char *form, ...)
#endif /* MACOSX */
   { 
      //sxprintf ("\nA BLAS/LAPACK error has occured!\n");
      std::cout << "\nA BLAS/LAPACK error has occured!\n";
      // --- original code from ATLAS/interfaces/blas/C/src/cblas_xerbla.c
      va_list argptr;
      va_start(argptr, form);
#     ifdef GCCWIN
         //if (p)  sxprintf("Parameter %d to routine %s was incorrect\n", p, rout);
         if (p) std::cout << "Parameter " << p << " to routine "
                          << rout << " was incorrect\n";
         vprintf(form, argptr);
#     else
         if (p)
            fprintf(stderr, "Parameter %d to routine %s was incorrect\n", 
                    p, rout);
         vfprintf(stderr, form, argptr);
#     endif /* CYGWIN */
      va_end(argptr);
      // --- end of original ATLAS code
      SX_EXIT;
   }

#endif /* USE_VECLIB */


//------------------------------------------------------------------------------
// norm of vectors
//------------------------------------------------------------------------------

// --- The next 96 lines were generated from snippets/SxBlasLib.cpp snippet norm2
float norm2 (const float *vec, ssize_t n, int incx)
{
   if (n < 1) return 0;
   int nn = int(n);
   if (nn < n)  {
      // BLAS indices limited to integer range: split vector into two
      ssize_t n2 = n / 2;
      double res1 = norm2(vec, n2, incx);
      double res2 = norm2(vec + n2 * incx, n - n2, incx);
      return (float)sqrt(res1*res1 + res2*res2);
   }
#  if   defined (USE_VECLIB)
      return snrm2 (&nn, const_cast<float*>(vec), &incx);
#  elif defined (USE_ESSL)
      return snrm2 (nn, vec, incx);
#  elif defined (USE_INTEL_MKL)
      return snrm2 (&nn, const_cast<float*>(vec), &incx);
#  elif defined (USE_ACML)
      return snrm2 (nn, const_cast<float*>(vec), incx);
#  else /* also: USE_GOTO */
      return cblas_snrm2 (nn, vec, incx);
#  endif
}

double norm2 (const double *vec, ssize_t n, int incx)
{
   if (n < 1) return 0;
   int nn = int(n);
   if (nn < n)  {
      // BLAS indices limited to integer range: split vector into two
      ssize_t n2 = n / 2;
      double res1 = norm2(vec, n2, incx);
      double res2 = norm2(vec + n2 * incx, n - n2, incx);
      return (double)sqrt(res1*res1 + res2*res2);
   }
#  if   defined (USE_VECLIB)
      return dnrm2 (&nn, const_cast<double*>(vec), &incx);
#  elif defined (USE_ESSL)
      return dnrm2 (nn, vec, incx);
#  elif defined (USE_INTEL_MKL)
      return dnrm2 (&nn, const_cast<double*>(vec), &incx);
#  elif defined (USE_ACML)
      return dnrm2 (nn, const_cast<double*>(vec), incx);
#  else /* also: USE_GOTO */
      return cblas_dnrm2 (nn, vec, incx);
#  endif
}

float norm2 (const SxComplex8 *vec, ssize_t n, int incx)
{
   if (n < 1) return 0;
   int nn = int(n);
   if (nn < n)  {
      // BLAS indices limited to integer range: split vector into two
      ssize_t n2 = n / 2;
      double res1 = norm2(vec, n2, incx);
      double res2 = norm2(vec + n2 * incx, n - n2, incx);
      return (float)sqrt(res1*res1 + res2*res2);
   }
#  if   defined (USE_VECLIB)
      return scnrm2 (&nn, (complex8_t*)const_cast<SxComplex8*>(vec), &incx);
#  elif defined (USE_ESSL)
      return scnrm2 (nn, (const complex<float>*)vec, incx);
#  elif defined (USE_INTEL_MKL)
      return scnrm2 (&nn, (MKL_Complex8*)const_cast<SxComplex8*>(vec), &incx);
#  elif defined (USE_ACML)
      return scnrm2 (nn, (complex*)const_cast<SxComplex8*>(vec), incx);
#  else /* also: USE_GOTO */
      return cblas_scnrm2 (nn, vec, incx);
#  endif
}

double norm2 (const SxComplex16 *vec, ssize_t n, int incx)
{
   if (n < 1) return 0;
   int nn = int(n);
   if (nn < n)  {
      // BLAS indices limited to integer range: split vector into two
      ssize_t n2 = n / 2;
      double res1 = norm2(vec, n2, incx);
      double res2 = norm2(vec + n2 * incx, n - n2, incx);
      return (double)sqrt(res1*res1 + res2*res2);
   }
#  if   defined (USE_VECLIB)
      return dznrm2 (&nn, (complex16_t*)const_cast<SxComplex16*>(vec), &incx);
#  elif defined (USE_ESSL)
      return dznrm2 (nn, (const complex<double>*)vec, incx);
#  elif defined (USE_INTEL_MKL)
      return dznrm2 (&nn, (MKL_Complex16*)const_cast<SxComplex16*>(vec), &incx);
#  elif defined (USE_ACML)
      return dznrm2 (nn, (doublecomplex*)const_cast<SxComplex16*>(vec), incx);
#  else /* also: USE_GOTO */
      return cblas_dznrm2 (nn, vec, incx);
#  endif
}
// --- norm2


//------------------------------------------------------------------------------
// scale vectors
//------------------------------------------------------------------------------
void scale (float *vec,  const float alpha, int n)
{
   int incx = 1;
#  if   defined (USE_VECLIB)
      sscal (&n, (float *)&alpha, (float *)vec, &incx);
#  elif defined (USE_ESSL)
      sscal (n, alpha, vec, incx);
#  elif defined (USE_INTEL_MKL)
      sscal (&n, const_cast<float *>(&alpha), (float *)vec, &incx);
#  elif defined (USE_ACML)
      sscal (n, alpha, vec, incx);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_sscal (n, alpha, vec, incx);
#  else
      cblas_sscal (n, alpha, vec, incx);
#  endif      
}


void scale (double *vec, const double alpha, int n)
{
   int incx = 1;
#  if   defined (USE_VECLIB)
      dscal (&n, (double *)&alpha, (double *)vec, &incx);
#  elif defined (USE_ESSL)
      dscal (n, alpha, vec, incx);
#  elif defined (USE_INTEL_MKL)
      dscal (&n, const_cast<double *>(&alpha), (double *)vec, &incx);
#  elif defined (USE_ACML)
      dscal (n, alpha, vec, incx);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_dscal (n, alpha, vec, incx);
#  else
      cblas_dscal (n, alpha, vec, incx);
#  endif      
}

void scale (SxComplex8 *vec,  const SxComplex8 &alpha, int n)
{
   int incx = 1;
#  if   defined (USE_VECLIB)
      cscal (&n, (complex8_t *)&alpha, (complex8_t *)vec, &incx);
#  elif defined (USE_ESSL)
      const complex<float> alphaTmp (alpha.re, alpha.im);
      cscal (n, alphaTmp, (complex<float> *)vec, incx);
#  elif defined (USE_INTEL_MKL)
      cscal (&n, (MKL_Complex8 *)const_cast<SxComplex8*>(&alpha),
             (MKL_Complex8 *)vec, &incx);
#  elif defined (USE_ACML)
      cscal (n, (complex *)const_cast<SxComplex8 *>(&alpha), 
             (complex *)vec, incx);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_cscal (n, (float*)&alpha, (float*)vec, incx);
#  else
      cblas_cscal (n, &alpha, vec, incx);
#  endif      
}


void scale (SxComplex16 *vec, const SxComplex16 &alpha, int n)
{
   int incx = 1;
#  if   defined (USE_VECLIB)
      zscal (&n, (complex16_t *)&alpha, (complex16_t *)vec, &incx);
#  elif defined (USE_ESSL)
      const complex<double> alphaTmp (alpha.re, alpha.im);
      zscal (n, alphaTmp, (complex<double> *)vec, incx);
#  elif defined (USE_INTEL_MKL)
      zscal (&n, (MKL_Complex16 *)const_cast<SxComplex16*>(&alpha),
             (MKL_Complex16 *)vec, &incx);
#  elif defined (USE_ACML)
      zscal (n, (doublecomplex *)const_cast<SxComplex16 *>(&alpha), 
                (doublecomplex *)vec, incx);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_zscal (n, (double*)&alpha, (double*)vec, incx);
#  else
      cblas_zscal (n, &alpha, vec, incx);
#  endif      
}


//------------------------------------------------------------------------------
// Y += a*X 
//------------------------------------------------------------------------------
// --- The next 58 lines were generated from snippets/SxBlasLib.cpp snippet axpy_real
void axpy (float *yOut, int incy, const float &alpha, const float *xIn, int incx, ssize_t n)
{
   if (n < 1) return;
   int nn = int(n);
   if (nn < n)  {
      // BLAS indices limited to integer range: split vector into two
      ssize_t n2 = n / 2;
      axpy (yOut, incy, alpha, xIn, incx, n2);
      axpy (yOut + n2 * incy, alpha, xIn + n2 * incx, n - n2);
   }
#  if   defined (USE_VECLIB)
      saxpy (&nn, (float *)&alpha,
             (float *)xIn, &incx, (float *)yOut, &incy);
#  elif defined (USE_ESSL)
      saxpy (nn, a,
            (float *)xIn, incx, (float *)yOut, incy);
#  elif defined (USE_INTEL_MKL)
      saxpy (&nn, const_cast<float *>(&alpha),
            const_cast<float *>(xIn), &incx, yOut, &incy);
#  elif defined (USE_ACML)
      saxpy (nn, alpha, (float *)const_cast<float *>(xIn),
             incx, (float *)yOut, incy);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_saxpy (nn, alpha, const_cast<float*>(xIn), incx, yOut, incy);
#  else
      cblas_saxpy (nn, alpha, xIn, incx, yOut, incy);
#  endif
}

void axpy (double *yOut, int incy, const double &alpha, const double *xIn, int incx, ssize_t n)
{
   if (n < 1) return;
   int nn = int(n);
   if (nn < n)  {
      // BLAS indices limited to integer range: split vector into two
      ssize_t n2 = n / 2;
      axpy (yOut, incy, alpha, xIn, incx, n2);
      axpy (yOut + n2 * incy, alpha, xIn + n2 * incx, n - n2);
   }
#  if   defined (USE_VECLIB)
      daxpy (&nn, (double *)&alpha,
             (double *)xIn, &incx, (double *)yOut, &incy);
#  elif defined (USE_ESSL)
      daxpy (nn, a,
            (double *)xIn, incx, (double *)yOut, incy);
#  elif defined (USE_INTEL_MKL)
      daxpy (&nn, const_cast<double *>(&alpha),
            const_cast<double *>(xIn), &incx, yOut, &incy);
#  elif defined (USE_ACML)
      daxpy (nn, alpha, (double *)const_cast<double *>(xIn),
             incx, (double *)yOut, incy);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_daxpy (nn, alpha, const_cast<double*>(xIn), incx, yOut, incy);
#  else
      cblas_daxpy (nn, alpha, xIn, incx, yOut, incy);
#  endif
}
// --- axpy_real
// --- The next 66 lines were generated from snippets/SxBlasLib.cpp snippet axpy_complex
void axpy (SxComplex8 *yOut, int incy,
           const SxComplex8 &alpha, const SxComplex8 *xIn, int incx, ssize_t n)
{
   if (n < 1) return;
   int nn = int(n);
   if (nn < n)  {
      // BLAS indices limited to integer range: split vector into two
      ssize_t n2 = n / 2;
      axpy (yOut, incy, alpha, xIn, incx, n2);
      axpy (yOut + n2 * incy, alpha, xIn + n2 * incx, n - n2);
   }
#  if   defined (USE_VECLIB)
      caxpy (&nn, (complex8_t *)&alpha,
             (complex8_t *)xIn, &incx, (complex8_t *)yOut, &incy);
#  elif defined (USE_ESSL)
      const complex<float> alphaTmp (alpha.re, alpha.im);
      caxpy (nn, alphaTmp,
            (complex<float> *)xIn, incx, (complex<float> *)yOut, incy);
#  elif defined (USE_INTEL_MKL)
      caxpy (&nn, (MKL_Complex8 *)const_cast<SxComplex8 *>(&alpha),
            (MKL_Complex8 *)const_cast<SxComplex8 *>(xIn), &incx,
            (MKL_Complex8 *)yOut, &incy);
#  elif defined (USE_ACML)
      caxpy (nn, (complex *)const_cast<SxComplex8 *>(&alpha),
                (complex *)const_cast<SxComplex8 *>(xIn), incx,
                (complex *)yOut, incy);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_caxpy (nn, (float*)&alpha, (float*)xIn, incx, (float*)yOut, incy);
#  else
      cblas_caxpy (nn, &alpha, xIn, incx, yOut, incy);
#  endif
}

void axpy (SxComplex16 *yOut, int incy,
           const SxComplex16 &alpha, const SxComplex16 *xIn, int incx, ssize_t n)
{
   if (n < 1) return;
   int nn = int(n);
   if (nn < n)  {
      // BLAS indices limited to integer range: split vector into two
      ssize_t n2 = n / 2;
      axpy (yOut, incy, alpha, xIn, incx, n2);
      axpy (yOut + n2 * incy, alpha, xIn + n2 * incx, n - n2);
   }
#  if   defined (USE_VECLIB)
      zaxpy (&nn, (complex16_t *)&alpha,
             (complex16_t *)xIn, &incx, (complex16_t *)yOut, &incy);
#  elif defined (USE_ESSL)
      const complex<double> alphaTmp (alpha.re, alpha.im);
      zaxpy (nn, alphaTmp,
            (complex<double> *)xIn, incx, (complex<double> *)yOut, incy);
#  elif defined (USE_INTEL_MKL)
      zaxpy (&nn, (MKL_Complex16 *)const_cast<SxComplex16 *>(&alpha),
            (MKL_Complex16 *)const_cast<SxComplex16 *>(xIn), &incx,
            (MKL_Complex16 *)yOut, &incy);
#  elif defined (USE_ACML)
      zaxpy (nn, (doublecomplex *)const_cast<SxComplex16 *>(&alpha),
                (doublecomplex *)const_cast<SxComplex16 *>(xIn), incx,
                (doublecomplex *)yOut, incy);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_zaxpy (nn, (double*)&alpha, (double*)xIn, incx, (double*)yOut, incy);
#  else
      cblas_zaxpy (nn, &alpha, xIn, incx, yOut, incy);
#  endif
}
// --- axpy_complex

//------------------------------------------------------------------------------
// scalar product
//------------------------------------------------------------------------------

float scalarProduct (const float *aVec, const float *bVec, int n)
{
   const float *aPtr = aVec;
   const float *bPtr = bVec;
   float res = 0.;
   for (int i=0; i < n; i++, aPtr++, bPtr++)
      res += *aPtr * *bPtr;
   return res;
}



double scalarProduct (const double *aVec, const double *bVec, int n)
{
   const double *aPtr = aVec;
   const double *bPtr = bVec;
   double res = 0;
   for (int i=0; i < n; i++, aPtr++, bPtr++)  
      res += *aPtr * *bPtr;
   return res;
}


SxComplex8 scalarProduct (const SxComplex8 *aVec, const SxComplex8 *bVec, int n)
{
   const SxComplex8 *aPtr = aVec;
   const SxComplex8 *bPtr = bVec;
   SxComplex8 res = (SxComplex8)0.;
   for (int i=0; i < n; i++, aPtr++, bPtr++)
      res += aPtr->conj() * *bPtr;
   return res;
}


SxComplex16 scalarProduct (const SxComplex16 *aVec, const SxComplex16 *bVec, int n)
{
   const SxComplex16 *aPtr = aVec;
   const SxComplex16 *bPtr = bVec;
   SxComplex16 res = (SxComplex8)0.;
   for (int i=0; i < n; i++, aPtr++, bPtr++)
      res += aPtr->conj() * *bPtr;
   return res;
}


//------------------------------------------------------------------------------
// general matrix-matrix multiplication
//------------------------------------------------------------------------------
void matmult (float *resMat, float beta, const float *aMat, const float *bMat,
              int aMatRows, int aMatCols, int bMatCols, int lda, int ldb, int ldc)
{
   float alpha = 1.0;
#  if   defined (USE_VECLIB)
      char noTrans = 'N';
      sgemm (&noTrans, &noTrans,
             &aMatRows, &bMatCols, &aMatCols,
             &alpha, (float *)aMat, &lda,
             (float *)bMat, &ldb, &beta, (float *)resMat, &ldc, 0, 0);
#  elif defined (USE_ESSL)
      const char noTrans = 'N';
      sgemm (&noTrans, &noTrans,
             aMatRows, bMatCols, aMatCols,
             alpha, aMat, lda,
             bMat, ldb, beta, resMat, ldc);
#  elif defined (USE_INTEL_MKL)
      char noTrans = 'N';
      sgemm (&noTrans, &noTrans, &aMatRows, &bMatCols, &aMatCols,
             &alpha, const_cast<float *>(aMat), &lda,
             const_cast<float *>(bMat), &ldb, &beta, resMat, &ldc);
#  elif defined (USE_ACML)
      char noTrans = 'N';
      sgemm (noTrans, noTrans, aMatRows, bMatCols, aMatCols,
             alpha, const_cast<float *>(aMat), lda,
             const_cast<float *>(bMat), ldb, beta, resMat, ldc);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_sgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                   aMatRows, bMatCols, aMatCols,
                   alpha, const_cast<float*>(aMat), lda,
                   const_cast<float*>(bMat), ldb, beta, resMat, ldc);
#  else
      cblas_sgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                   aMatRows, bMatCols, aMatCols,
                   alpha, aMat, lda,
                   bMat, ldb, beta, resMat, ldc);
#  endif
}


void matmult (double *resMat, double beta, const double *aMat, const double *bMat,
              int aMatRows, int aMatCols, int bMatCols, int lda, int ldb, int ldc)
{
   double alpha = 1.0;
#  if   defined (USE_VECLIB)
      char noTrans = 'N';
      dgemm (&noTrans, &noTrans,
             &aMatRows, &bMatCols, &aMatCols,
             &alpha, (double *)aMat, &lda,
             (double *)bMat, &ldb, &beta, (double *)resMat, &ldc, 0, 0);
#  elif defined (USE_ESSL)
      const char noTrans = 'N';
      dgemm (&noTrans, &noTrans,
             aMatRows, bMatCols, aMatCols,
             alpha, aMat, lda,
             bMat, ldb, beta, resMat, ldc);
#  elif defined (USE_INTEL_MKL)
      char noTrans = 'N';
      dgemm (&noTrans, &noTrans,
             &aMatRows, &bMatCols, &aMatCols,
             &alpha, const_cast<double *>(aMat), &lda,
             const_cast<double *>(bMat), &ldb, &beta, resMat, &ldc);
#  elif defined (USE_ACML)
      char noTrans = 'N';
      dgemm (noTrans, noTrans,
             aMatRows, bMatCols, aMatCols,
             alpha, const_cast<double *>(aMat), lda,
             const_cast<double *>(bMat), ldb, beta, resMat, ldc);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                         aMatRows, bMatCols, aMatCols,
                         alpha, const_cast<double*>(aMat), lda,
                         const_cast<double*>(bMat), ldb, beta, resMat, ldc);
#  else
      cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                   aMatRows, bMatCols, aMatCols,
                   alpha, aMat, lda,
                   bMat, ldb, beta, resMat, ldc);
#  endif
}


void matmult (SxComplex8 *resMat, const SxComplex8 &beta,
              const SxComplex8 *aMat, const SxComplex8 *bMat,
              int aMatRows, int aMatCols, int bMatCols,
              int lda, int ldb, int ldc)
{
   SxComplex8 alpha (1.0, 0.0);
#  if   defined (USE_VECLIB)
      char noTrans = 'N';
      cgemm (&noTrans, &noTrans,
             &aMatRows, &bMatCols, &aMatCols,
             (complex8_t *)&alpha, (complex8_t *)aMat, &lda,
             (complex8_t *)bMat,   &ldb, (const complex8_t *)&beta,
             (complex8_t *)resMat, &ldc, 0, 0);
#  elif defined (USE_ESSL)
      const char noTrans = 'N';
      complex<float> alphaT (alpha.re, alpha.im);
      complex<float> betaT  (beta.re,  beta.im);
      cgemm (&noTrans, &noTrans,
             aMatRows, bMatCols, aMatCols,
             alphaT, (complex<float> *)aMat, lda,
             (complex<float> *)bMat,   ldb, betaT,
             (complex<float> *)resMat, ldc);
#  elif defined (USE_INTEL_MKL)
      char noTrans = 'N';
      cgemm (&noTrans, &noTrans, &aMatRows, &bMatCols, &aMatCols,
             (const MKL_Complex8 *)&alpha,
             (const MKL_Complex8 *)aMat, &lda,
             (const MKL_Complex8 *)bMat, &ldb,
             (const MKL_Complex8 *)&beta,
             (MKL_Complex8 *)resMat, &ldc);
#  elif defined (USE_ACML)
      char noTrans = 'N';
      cgemm (noTrans, noTrans,
             aMatRows, bMatCols, aMatCols,
             (complex *)const_cast<SxComplex8 *>(&alpha),
             (complex *)const_cast<SxComplex8 *>(aMat), lda,
             (complex *)const_cast<SxComplex8 *>(bMat), ldb,
             (complex *)const_cast<SxComplex8 *>(&beta),
             (complex *)resMat, ldc);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      cblas_cgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                   aMatRows, bMatCols, aMatCols,
                   (const float*)&alpha, (const float*)aMat, lda,
                   (const float*)bMat, ldb, (const float*)&beta,
                   (float*)resMat, ldc);
#  else
      cblas_cgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                   aMatRows, bMatCols, aMatCols,
                   &alpha, aMat, lda,
                   bMat, ldb, &beta, resMat, ldc);
#  endif
}


void matmult (SxComplex16 *resMat, const SxComplex16 &beta,
              const SxComplex16 *aMat, const SxComplex16 *bMat,
              int aMatRows, int aMatCols, int bMatCols,
              int lda, int ldb, int ldc)
{
   SxComplex16 alpha (1.0, 0.0);
#  if   defined (USE_VECLIB)
      char noTrans = 'N';
      zgemm (&noTrans, &noTrans,
             &aMatRows, &bMatCols, &aMatCols,
             (complex16_t *)&alpha, (complex16_t *)aMat, &lda,
             (complex16_t *)bMat,   &ldb, (complex16_t *)&beta,
             (complex16_t *)resMat, &ldc, 0, 0);
#  elif defined (USE_ESSL)
      const char noTrans = 'N';
      complex<double> alphaT (alpha.re, alpha.im);
      complex<double> betaT  (beta.re,  beta.im);
      zgemm (&noTrans, &noTrans,
             aMatRows, bMatCols, aMatCols,
             alphaT, (complex<double> *)aMat, lda,
             (complex<double> *)bMat,   ldb, betaT,
             (complex<double> *)resMat, ldc);
#  elif defined (USE_INTEL_MKL)
      char noTrans = 'N';
      zgemm (&noTrans, &noTrans, &aMatRows, &bMatCols, &aMatCols,
             (const MKL_Complex16 *)&alpha,
             (const MKL_Complex16 *)aMat, &lda,
             (const MKL_Complex16 *)bMat, &ldb,
             (const MKL_Complex16 *)&beta, (MKL_Complex16 *)resMat, &ldc);
#  elif defined (USE_ACML)
      char noTrans = 'N';
      zgemm (noTrans, noTrans,
             aMatRows, bMatCols, aMatCols,
             (doublecomplex *)const_cast<SxComplex16 *>(&alpha),
             (doublecomplex *)const_cast<SxComplex16 *>(aMat), lda,
             (doublecomplex *)const_cast<SxComplex16 *>(bMat), ldb,
             (doublecomplex *)const_cast<SxComplex16 *>(&beta),
             (doublecomplex *)resMat, ldc);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      if (bMatCols == 1)  {
         cblas_zgemv (CblasColMajor, CblasNoTrans, aMatRows, aMatCols,
                      (double*)&alpha, (double*)aMat, lda,
                      (double*)bMat, 1, (double*)&beta, (double*)resMat, 1);
      } else if (aMatRows == 1)  {
         cblas_zgemv (CblasColMajor, CblasTrans, aMatCols, bMatCols,
                      (double*)&alpha, (double*)bMat, ldb,
                      (double*)aMat, 1, (double*)&beta, (double*)resMat, 1);
      } else {
         cblas_zgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                      aMatRows, bMatCols, aMatCols,
                      (double*)&alpha, (double*)aMat, lda,
                      (double*)bMat, ldb, (double*)&beta, (double*)resMat,
                      ldc);
      }
#  else
      if (bMatCols == 1)  {
         cblas_zgemv (CblasColMajor, CblasNoTrans, aMatRows, aMatCols,
                      &alpha, aMat, lda,
                      bMat, 1, &beta, resMat, 1);
      } else if (aMatRows == 1)  {
         cblas_zgemv (CblasColMajor, CblasTrans, aMatCols, bMatCols,
                      &alpha, bMat, ldb,
                      aMat, 1, &beta, resMat, 1);
      } else {
#if defined(USE_OPENMP) && defined(USE_ATLAS)
         if (aMatRows > 1024)  {
#           pragma omp parallel
            {
               int nThreads = omp_get_num_threads ();
               int nb = aMatRows / nThreads;
               // distribute remaining elements over all threads
               if (nb * nThreads < aMatRows) nb++;
               int offset = omp_get_thread_num () * nb;
               // reduce nb for last thread
               if (offset + nb > aMatRows) nb = aMatRows - offset;
               cblas_zgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                            nb, bMatCols, aMatCols,
                            &alpha, aMat + offset, lda,
                            bMat, ldb, &beta, resMat + offset, ldc);
            }
         } else // no openMP parallelism ...
#endif
         cblas_zgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                      aMatRows, bMatCols, aMatCols,
                      &alpha, aMat, lda,
                      bMat, ldb, &beta, resMat, ldc);
      }
#  endif
}

//------------------------------------------------------------------------------
// overlap matrices
//------------------------------------------------------------------------------
void matovlp (float *resMat, 
              const float *aMat, const float *bMat, 
              int aMatRows, int aMatCols, int bMatRows, int bMatCols,
              int sumSize)
{
   float alpha = 1.0, beta = 0.0;
#  if   defined (USE_VECLIB)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      sgemm (&conjTrans, &noTrans, 
             &aMatCols, &bMatCols, &sumSize,
             &alpha, (float *)aMat, &aMatRows,
             (float *)bMat,   &bMatRows, &beta, 
             (float *)resMat, &aMatCols, 0, 0);
#  elif defined (USE_ESSL)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      sgemm (&conjTrans, &noTrans, 
             aMatCols, bMatCols, sumSize,
             alpha, aMat, aMatRows,
             bMat, bMatRows, beta, resMat, aMatCols);
#  elif defined (USE_INTEL_MKL)
//      SX_EXIT; // not tested
      // change by khr
      char noTrans = 'N', conjTrans = 'C';
      sgemm (&conjTrans, &noTrans, 
             &aMatCols, &bMatCols, &sumSize,
             &alpha, const_cast<float *>(aMat), &aMatRows,
             const_cast<float *>(bMat), &bMatRows, &beta, resMat, &aMatCols);
#  elif defined (USE_ACML)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      sgemm (conjTrans, noTrans, 
             aMatCols, bMatCols, sumSize,
             alpha, const_cast<float *>(aMat), aMatRows,
             const_cast<float *>(bMat), bMatRows, 
             beta, resMat, aMatCols);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      if (bMatCols == 1)  {
            cblas_sgemv (CblasColMajor, CblasConjTrans, sumSize, aMatCols,
                         alpha, const_cast<float*>(aMat), aMatRows,
                         const_cast<float*>(bMat), 1, beta, resMat, 1);
         } else if (aMatCols == 1)  {
            cblas_sgemv (CblasColMajor, CblasConjTrans, sumSize, bMatCols,
                         alpha, const_cast<float*>(bMat), bMatRows,
                         const_cast<float*>(aMat), 1, beta, resMat, 1);
         } else {
            cblas_sgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                         aMatCols, bMatCols, sumSize,
                         alpha, const_cast<float *>(aMat), aMatRows,
                         const_cast<float*>(bMat), bMatRows, beta, resMat, aMatCols);
      }
#  else  
      if (bMatCols == 1)  {
         cblas_sgemv (CblasColMajor, CblasConjTrans, sumSize, aMatCols,
                      alpha, aMat, aMatRows,
                      bMat, 1, beta, resMat, 1);
      } else if (aMatCols == 1)  {
         cblas_sgemv (CblasColMajor, CblasConjTrans, sumSize, bMatCols,
                      alpha, bMat, bMatRows,
                      aMat, 1, beta, resMat, 1);
      } else {
         cblas_sgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                      aMatCols, bMatCols, sumSize, 
                      alpha, aMat, aMatRows, 
                      bMat, bMatRows, beta, resMat, aMatCols);
      }
#  endif      
}

void matovlp (double *resMat, 
              const double *aMat, const double *bMat, 
              int aMatRows, int aMatCols, int bMatRows, int bMatCols,
              int sumSize)
{
   double alpha = 1.0, beta = 0.0;
#  if   defined (USE_VECLIB)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      dgemm (&conjTrans, &noTrans, 
             &aMatCols, &bMatCols, &sumSize,
             &alpha, (double *)aMat, &aMatRows,
             (double *)bMat,   &bMatRows, &beta, 
             (double *)resMat, &aMatCols, 0, 0);
#  elif defined (USE_ESSL)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      dgemm (&conjTrans, &noTrans, 
             aMatCols, bMatCols, sumSize,
             alpha, aMat, aMatRows,
             bMat, bMatRows, beta, resMat, aMatCols);
#  elif defined (USE_INTEL_MKL)
//      SX_EXIT; // not tested
      // change by khr
      char noTrans = 'N', conjTrans = 'C';
      dgemm (&conjTrans, &noTrans, 
             &aMatCols, &bMatCols, &sumSize,
             &alpha, const_cast<double *>(aMat), &aMatRows,
             const_cast<double *>(bMat), &bMatRows, &beta, resMat, &aMatCols);
#  elif defined (USE_ACML)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      dgemm (conjTrans, noTrans, 
             aMatCols, bMatCols, sumSize,
             alpha, const_cast<double *>(aMat), aMatRows,
             const_cast<double *>(bMat), bMatRows, 
             beta, resMat, aMatCols);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      if (bMatCols == 1)  {
         cblas_dgemv (CblasColMajor, CblasConjTrans, sumSize, aMatCols,
                      alpha, const_cast<double*>(aMat), aMatRows,
                      const_cast<double*>(bMat), 1, beta, resMat, 1);
      } else if (aMatCols == 1)  {
         cblas_dgemv (CblasColMajor, CblasConjTrans, sumSize, bMatCols,
                      alpha, const_cast<double*>(bMat), bMatRows,
                      const_cast<double*>(aMat), 1, beta, resMat, 1);
      } else {
         cblas_dgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                      aMatCols, bMatCols, sumSize,
                      alpha, const_cast<double*>(aMat), aMatRows,
                      const_cast<double*>(bMat), bMatRows, beta, resMat, aMatCols);
      }
#  else  
      if (bMatCols == 1)  {
         cblas_dgemv (CblasColMajor, CblasConjTrans, sumSize, aMatCols,
                      alpha, aMat, aMatRows,
                      bMat, 1, beta, resMat, 1);
      } else if (aMatCols == 1)  {
         cblas_dgemv (CblasColMajor, CblasConjTrans, sumSize, bMatCols,
                      alpha, bMat, bMatRows,
                      aMat, 1, beta, resMat, 1);
      } else {
         cblas_dgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                      aMatCols, bMatCols, sumSize, 
                      alpha, aMat, aMatRows, 
                      bMat, bMatRows, beta, resMat, aMatCols);
      }
#  endif      
}

void matovlp (SxComplex8 *resMat, 
              const SxComplex8 *aMat, const SxComplex8 *bMat, 
              int aMatRows, int aMatCols, int bMatRows, int bMatCols,
              int sumSize)
{
   SxComplex8 alpha (1.0, 0.0), beta (0.0, 0.0);
#  if   defined (USE_VECLIB)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      cgemm (&conjTrans, &noTrans, 
             &aMatCols, &bMatCols, &sumSize,
             (complex8_t *)&alpha, (complex8_t *)aMat, &aMatRows,
             (complex8_t *)bMat,   &bMatRows, (complex8_t *)&beta, 
             (complex8_t *)resMat, &aMatCols, 0, 0);
#  elif defined (USE_ESSL)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      complex<float> alphaT (alpha.re, alpha.im);
      complex<float> betaT  (beta.re,  beta.im);
      cgemm (&conjTrans, &noTrans, 
             aMatCols, bMatCols, sumSize,
             alphaT, (complex<float> *)aMat, aMatRows,
             (complex<float> *)bMat,   bMatRows, betaT, 
             (complex<float> *)resMat, aMatCols);
#  elif defined (USE_INTEL_MKL)
      char noTrans = 'N', conjTrans = 'C';
      cgemm (&conjTrans, &noTrans, &aMatCols, &bMatCols, &sumSize,
             (MKL_Complex8 *)&alpha,
             (MKL_Complex8 *)const_cast<SxComplex8 *>(aMat), &aMatRows,
             (MKL_Complex8 *)const_cast<SxComplex8 *>(bMat), &bMatRows,
             (MKL_Complex8 *)&beta, (MKL_Complex8 *)resMat, &aMatCols);
#  elif defined (USE_ACML)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      cgemm (conjTrans, noTrans, 
             aMatCols, bMatCols, sumSize,
             (complex *)const_cast<SxComplex8 *>(&alpha), 
             (complex *)const_cast<SxComplex8 *>(aMat), aMatRows,
             (complex *)const_cast<SxComplex8 *>(bMat), bMatRows, 
             (complex *)const_cast<SxComplex8 *>(&beta), 
             (complex *)resMat, aMatCols);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      if (bMatCols == 1)  {
         cblas_cgemv (CblasColMajor, CblasConjTrans, sumSize, aMatCols,
                      (float*)&alpha, (float*)aMat, aMatRows,
                      (float*)bMat, 1, (float*)&beta, (float*)resMat, 1);
      } else if (aMatCols == 1)  {
         cblas_cgemv (CblasColMajor, CblasConjTrans, sumSize, bMatCols,
                      (float*)&alpha, (float*)bMat, bMatRows,
                      (float*)aMat, 1, (float*)&beta, (float*)resMat, 1);
         // now resMat contains (B.adjoint () ^ A)
         // perform conjugate (note: res is a vector)
         for (ssize_t i = 0; i < bMatCols; ++i)
            resMat[i].im = -resMat[i].im;
      } else {
         cblas_cgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                      aMatCols, bMatCols, sumSize,
                      (float*)&alpha, (float*)aMat, aMatRows,
                      (float*)bMat, bMatRows, (float*)&beta, (float*)resMat,
                      aMatCols);
      }
#  else  
      if (bMatCols == 1)  {
         cblas_cgemv (CblasColMajor, CblasConjTrans, sumSize, aMatCols,
                      &alpha, aMat, aMatRows,
                      bMat, 1, &beta, resMat, 1);
      } else if (aMatCols == 1)  {
         cblas_cgemv (CblasColMajor, CblasConjTrans, sumSize, bMatCols,
                      &alpha, bMat, bMatRows,
                      aMat, 1, &beta, resMat, 1);
         // now resMat contains (B.adjoint () ^ A)
         // perform conjugate (note: res is a vector)
         for (ssize_t i = 0; i < bMatCols; ++i)
            resMat[i].im = -resMat[i].im;
      } else {
         cblas_cgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                      aMatCols, bMatCols, sumSize, 
                      &alpha, aMat, aMatRows, 
                      bMat, bMatRows, &beta, resMat, aMatCols);
      }
#  endif      
}

/*
#include <SxTimer.h>
enum BlasLibTimer { Matovlp };
REGISTER_TIMERS (BlasLibTimer)
{
   regTimer (Matovlp, "matovlp");
}
*/

void matovlp (SxComplex16 *resMat, 
              const SxComplex16 *aMat, const SxComplex16 *bMat, 
              int aMatRows, int aMatCols, int bMatRows, int bMatCols,
              int sumSize)
{
   SxComplex16 alpha (1.0, 0.0), beta (0.0, 0.0);
#  if   defined (USE_VECLIB)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      zgemm (&conjTrans, &noTrans, 
             &aMatCols, &bMatCols, &sumSize,
             (complex16_t *)&alpha, (complex16_t *)aMat, &aMatRows,
             (complex16_t *)bMat,   &bMatRows, (complex16_t *)&beta, 
             (complex16_t *)resMat, &aMatCols, 0, 0);
#  elif defined (USE_ESSL)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      complex<double> alphaT (alpha.re, alpha.im);
      complex<double> betaT  (beta.re,  beta.im);
      zgemm (&conjTrans, &noTrans, 
             aMatCols, bMatCols, sumSize,
             alphaT, (complex<double> *)aMat, aMatRows,
             (complex<double> *)bMat,   bMatRows, betaT, 
             (complex<double> *)resMat, aMatCols);
#  elif defined (USE_INTEL_MKL)
      char noTrans = 'N', conjTrans = 'C';
      if (aMat == bMat && aMatRows == bMatRows && aMatCols == bMatCols)  {
         // special case: A.overlap (A) => use Hermitean routine ...
         char uplo = 'U';
         zherk (&uplo, &conjTrans, &aMatCols, &sumSize,
                &alpha.re,
                (MKL_Complex16 *)const_cast<SxComplex16 *>(aMat), &aMatRows,
                &beta.re, (MKL_Complex16 *)resMat, &aMatCols);
         // ... and copy upper half to lower half
         for (int i = 0; i < aMatCols; i++)
            for (int j = 0; j < i; j++)
               resMat[i + aMatCols * j] = resMat[j + aMatCols * i].conj ();
      } else {
         zgemm (&conjTrans, &noTrans, &aMatCols, &bMatCols, &sumSize,
                (MKL_Complex16 *)&alpha,
                (MKL_Complex16 *)const_cast<SxComplex16 *>(aMat), &aMatRows,
                (MKL_Complex16 *)const_cast<SxComplex16 *>(bMat), &bMatRows,
                (MKL_Complex16 *)&beta, (MKL_Complex16 *)resMat, &aMatCols);
      }
#  elif defined (USE_ACML)
      SX_EXIT; // not tested
      char noTrans = 'N', conjTrans = 'C';
      zgemm (conjTrans, noTrans, 
             aMatCols, bMatCols, sumSize,
             (doublecomplex *)const_cast<SxComplex16 *>(&alpha), 
             (doublecomplex *)const_cast<SxComplex16 *>(aMat), aMatRows,
             (doublecomplex *)const_cast<SxComplex16 *>(bMat), bMatRows, 
             (doublecomplex *)const_cast<SxComplex16 *>(&beta), 
             (doublecomplex *)resMat, aMatCols);
#  elif defined (USE_GOTO) // khr: GotoBLAS, experimental!
      if (bMatCols == 1)  {
            cblas_zgemv (CblasColMajor, CblasConjTrans, sumSize, aMatCols,
                         (double*)&alpha, (double*)aMat, aMatRows,
                         (double*)bMat, 1, (double*)&beta, (double*)resMat, 1);
         } else if (aMatCols == 1)  {
            cblas_zgemv (CblasColMajor, CblasConjTrans, sumSize, bMatCols,
                         (double*)&alpha, (double*)bMat, bMatRows,
                         (double*)aMat, 1, (double*)&beta, (double*)resMat, 1);
            // now resMat contains (B.adjoint () ^ A)
            // perform conjugate (note: res is a vector)
            for (ssize_t i = 0; i < bMatCols; ++i)
               resMat[i].im = -resMat[i].im;
         } else {
            cblas_zgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                         aMatCols, bMatCols, sumSize,
                         (double*)&alpha, (double*)aMat, aMatRows,
                         (double*)bMat, bMatRows, (double*)&beta,
                         (double*)resMat, aMatCols);
      }
#  else  
      if (bMatCols == 1)  {
         cblas_zgemv (CblasColMajor, CblasConjTrans, sumSize, aMatCols,
                      &alpha, aMat, aMatRows,
                      bMat, 1, &beta, resMat, 1);
      } else if (aMatCols == 1)  {
         cblas_zgemv (CblasColMajor, CblasConjTrans, sumSize, bMatCols,
                      &alpha, bMat, bMatRows,
                      aMat, 1, &beta, resMat, 1);
         // now resMat contains (B.adjoint () ^ A)
         // perform conjugate (note: res is a vector)
         for (ssize_t i = 0; i < bMatCols; ++i)
            resMat[i].im = -resMat[i].im;
      } else if (aMat == bMat && aMatRows == bMatRows && aMatCols == bMatCols)  {
         // special case: A.overlap (A) => use Hermitean routine ...
         cblas_zherk (CblasColMajor, CblasUpper, CblasConjTrans, aMatCols, sumSize,
                      1.0, aMat, aMatRows,
                      0.0, resMat, aMatCols);
         // ... and copy upper half to lower half
         for (int i = 0; i < aMatCols; i++)
            for (int j = 0; j < i; j++)
               resMat[i + aMatCols * j] = resMat[j + aMatCols * i].conj ();
      } else {
#ifndef USE_OPENMP
         cblas_zgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                      aMatCols, bMatCols, sumSize, 
                      &alpha, aMat, aMatRows, 
                      bMat, bMatRows, &beta, resMat, aMatCols);
#else
         //CLOCK (Matovlp);
         int nb = 128;
         int np = sumSize / nb;

         // --- compute contribution from rest elements
         if (sumSize > nb * np)  {
            cblas_zgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                         aMatCols, bMatCols, sumSize - nb * np, 
                         &alpha, aMat + nb * np, aMatRows, 
                         bMat + nb * np, bMatRows, &beta, resMat, aMatCols);
         } else {
            for (int i = 0; i < aMatCols * bMatCols; ++i)
               resMat[i].re = resMat[i].im = 0.;

         }
#        pragma omp parallel
         {
            SxComplex16 *part = NULL;
            // --- compute partial results (in parallel)
#           pragma omp for
            for (int ip = 0; ip < np; ip++) {
               if (!part)  {
                  part = new SxComplex16[aMatCols * bMatCols];
                  cblas_zgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                               aMatCols, bMatCols, nb, 
                               &alpha, aMat + ip * nb, aMatRows, 
                               bMat + ip * nb, bMatRows, &beta, part, aMatCols);
               } else {
                  cblas_zgemm (CblasColMajor, CblasConjTrans, CblasNoTrans,
                               aMatCols, bMatCols, nb, 
                               &alpha, aMat + ip * nb, aMatRows, 
                               bMat + ip * nb, bMatRows, &alpha, part, aMatCols);
               }
            }
            // --- sum partial results
            if (part)  {
#              pragma omp critical
               cblas_zaxpy (aMatCols * bMatCols, &alpha, part, 1, resMat, 1);
               delete [] part;
            }
         }
#     endif
      }
#  endif      
}

//------------------------------------------------------------------------------
// Matrix decompositions
//------------------------------------------------------------------------------
void cholesky (float *resMat, enum UPLO uplo, float * /*inMat*/, int n)
{
   char uploChar = (uplo == UpperRight) ? 'U' : 'L';
#  if   defined (USE_VECLIB)
      int err = 0;
      spotrf (&uploChar, &n, (float *)resMat, &n, &err, 0);
#  elif defined (USE_ESSL)
      int err = 0;
      spotrf (&uploChar, n, resMat, n, err);
#  elif defined (USE_INTEL_MKL)
      int err = 0;      
      spotrf (&uploChar, &n, resMat, &n, &err);
#  elif defined (USE_ACML)
      int err = 0;      
      spotrf (uploChar, n, resMat, n, &err);
#  elif defined (USE_NETLIB)
      int err = LAPACKE_spotrf (LAPACK_COL_MAJOR, uploChar, n, resMat, n);
#  else 
      integer rank  = (integer)n, err = 0;
      spotrf_ (&uploChar, &rank, (real *)resMat, &rank, &err);
#   endif      
   if ( err )  { // TODO: throw exception
      std::cout << "cholesky err=" << err << std::endl;
      SX_EXIT;
   }
}


void cholesky (double *resMat, enum UPLO uplo, double * /*inMat*/, int n)
{
   char uploChar = (uplo == UpperRight) ? 'U' : 'L';
#  if   defined (USE_VECLIB)
      int err = 0;
      dpotrf (&uploChar, &n, (double *)resMat, &n, &err, 0);
#  elif defined (USE_ESSL)
      int err = 0;
      dpotrf (&uploChar, n, resMat, n, err);
#  elif defined (USE_INTEL_MKL)
      int err = 0;      
      dpotrf (&uploChar, &n, resMat, &n, &err);
#  elif defined (USE_ACML)
      int err = 0;      
      dpotrf (uploChar, n, resMat, n, &err);
#  elif defined (USE_NETLIB)
      int err = LAPACKE_dpotrf (LAPACK_COL_MAJOR, uploChar, n, resMat, n);
#  else 
      integer rank  = (integer)n, err = 0;
      dpotrf_ (&uploChar, &rank, (doublereal *)resMat, &rank, &err);
#  endif      
      if ( err )  { // TODO: throw exception
         std::cout << "cholesky err=" << err << std::endl;
         SX_EXIT;
      }
}


void cholesky (SxComplex8 *resMat, enum UPLO uplo, 
               SxComplex8 * /*inMat*/, int n)
{
   char uploChar = (uplo == UpperRight) ? 'U' : 'L';
#  if   defined (USE_VECLIB)
      int err = 0;
      cpotrf (&uploChar, &n, (complex8_t *)resMat, &n, &err, 0);
#  elif defined (USE_ESSL)
      int err = 0;
      cpotrf (&uploChar, n, (complex<float> *)resMat, n, err);
#  elif defined (USE_INTEL_MKL)
      int err = 0;      
      cpotrf (&uploChar, &n, (MKL_Complex8 *)resMat, &n, &err);
#  elif defined (USE_ACML)
      int err = 0;      
      cpotrf (uploChar, n, (complex *)resMat, n, &err);
#  elif defined (USE_NETLIB)
      int err = LAPACKE_cpotrf (LAPACK_COL_MAJOR, uploChar, n, resMat, n);
#  else 
      integer rank  = (integer)n, err = 0;
      cpotrf_ (&uploChar, &rank, (complex *)resMat, &rank, &err);
#  endif      
   if ( err )  { // TODO: throw exception
      std::cout << "cholesky err=" << err << std::endl;
      SX_EXIT;
   }
}


void cholesky (SxComplex16 *resMat, enum UPLO uplo, 
               SxComplex16 * /*inMat*/, int n)
{
   char uploChar = (uplo == UpperRight) ? 'U' : 'L';
#  if   defined (USE_VECLIB)
      int err = 0;
      zpotrf (&uploChar, &n, (complex16_t *)resMat, &n, &err, 0);
#  elif defined (USE_ESSL)
      int err = 0;
      zpotrf (&uploChar, n, (complex<double> *)resMat, n, err);
#  elif defined (USE_INTEL_MKL)
      int err = 0;      
      zpotrf(&uploChar, &n, (MKL_Complex16 *)resMat, &n, &err);
#  elif defined (USE_ACML)
      int err = 0;      
      zpotrf(uploChar, n, (doublecomplex *)resMat, n, &err);
#  elif defined (USE_NETLIB)
      int err = LAPACKE_zpotrf (LAPACK_COL_MAJOR, uploChar, n, resMat, n);
#  else 
      integer rank  = (integer)n, err = 0;
      zpotrf_ (&uploChar, &rank, (doublecomplex *)resMat, &rank, &err);
#  endif   
   if ( err )  { // TODO: throw exception
      std::cout << "cholesky err=" << err << std::endl;
      SX_EXIT;
   }
}

// --- The next 208 lines were generated from math/snippets/SxBlasLib.cpp snippet SVDREAL
void singularValueDecomp (float *mat, int nRows, int nCols,
                          float *vals,
                          float *left,
                          float *right, // V^H
                          bool zeroSpace)
{
   char jobz = zeroSpace ? 'A'  // U= M x M, V=N x N
                         : 'S'; // U=M x s, V= s x N
   if (!left || !right)  {
      jobz = 'N';
      if (left || right) {
         // It is not possible to compute only left or only right vectors.
         std::cout << "Internal error in singular-value decomposition"
                   << std::endl;
         SX_EXIT;
      }
   }
   int minMN = nRows < nCols ? nRows : nCols;
   int ldvt = zeroSpace ? nCols : minMN;
#  if   defined (USE_VECLIB)
   SX_EXIT; // not yet implemented
#  elif defined (USE_ESSL)
   SX_EXIT; // not yet implemented
#  elif defined (USE_INTEL_MKL)
   // workspace query:
   int *iwork = new int[8 * minMN];
   int err = 0;
   int lwork = -1;
   float optwork;
   sgesdd(&jobz, &nRows, &nCols, mat, &nRows,
          vals, left, &nRows, right, &ldvt,
          &optwork, &lwork, iwork, &err );
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = int(lround(optwork));
   float* work = new float[lwork];
   // --- actual compute
   sgesdd(&jobz, &nRows, &nCols, mat, &nRows,
          vals, left, &nRows, right, &ldvt,
          work, &lwork, iwork, &err );
   delete[] work;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  elif defined (USE_ACML)
   SX_EXIT;
#  elif defined (USE_NETLIB)
   // workspace query:
   int *iwork = new int[8 * minMN];
   int err = 0;
   int lwork = -1;
   float optwork;
   err = LAPACKE_sgesdd_work(LAPACK_COL_MAJOR,
          jobz, nRows, nCols, mat, nRows,
          vals, left, nRows, right, ldvt,
          &optwork, lwork, iwork);
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = (int)optwork;
   float* work = new float[lwork];
   // --- actual compute
   err = LAPACKE_sgesdd_work(LAPACK_COL_MAJOR,
          jobz, nRows, nCols, mat, nRows,
          vals, left, nRows, right, ldvt,
          work, lwork, iwork);
   delete[] work;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  else
   integer *iwork = new integer[8 * minMN];
   integer lwork = -1;
   integer err = 0;
   // workspace query:
   float optwork;
   integer nr = nRows, nc = nCols, ldvt_ = ldvt;
   sgesdd_(&jobz, &nr, &nc, mat, &nr,
          vals, left, &nr, right, &ldvt_,
          &optwork, &lwork, iwork, &err );
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = (integer)optwork;
   float* work = new float[lwork];
   // --- actual compute
   sgesdd_(&jobz, &nr, &nc, mat, &nr,
          vals, left, &nr, right, &ldvt_,
          work, &lwork, iwork, &err );
   delete[] work;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  endif   
   delete[] iwork;
}

void singularValueDecomp (double *mat, int nRows, int nCols,
                          double *vals,
                          double *left,
                          double *right, // V^H
                          bool zeroSpace)
{
   char jobz = zeroSpace ? 'A'  // U= M x M, V=N x N
                         : 'S'; // U=M x s, V= s x N
   if (!left || !right)  {
      jobz = 'N';
      if (left || right) {
         // It is not possible to compute only left or only right vectors.
         std::cout << "Internal error in singular-value decomposition"
                   << std::endl;
         SX_EXIT;
      }
   }
   int minMN = nRows < nCols ? nRows : nCols;
   int ldvt = zeroSpace ? nCols : minMN;
#  if   defined (USE_VECLIB)
   SX_EXIT; // not yet implemented
#  elif defined (USE_ESSL)
   SX_EXIT; // not yet implemented
#  elif defined (USE_INTEL_MKL)
   // workspace query:
   int *iwork = new int[8 * minMN];
   int err = 0;
   int lwork = -1;
   double optwork;
   dgesdd(&jobz, &nRows, &nCols, mat, &nRows,
          vals, left, &nRows, right, &ldvt,
          &optwork, &lwork, iwork, &err );
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = int(lround(optwork));
   double* work = new double[lwork];
   // --- actual compute
   dgesdd(&jobz, &nRows, &nCols, mat, &nRows,
          vals, left, &nRows, right, &ldvt,
          work, &lwork, iwork, &err );
   delete[] work;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  elif defined (USE_ACML)
   SX_EXIT;
#  elif defined (USE_NETLIB)
   // workspace query:
   int *iwork = new int[8 * minMN];
   int err = 0;
   int lwork = -1;
   double optwork;
   err = LAPACKE_dgesdd_work(LAPACK_COL_MAJOR,
          jobz, nRows, nCols, mat, nRows,
          vals, left, nRows, right, ldvt,
          &optwork, lwork, iwork);
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = (int)optwork;
   double* work = new double[lwork];
   // --- actual compute
   err = LAPACKE_dgesdd_work(LAPACK_COL_MAJOR,
          jobz, nRows, nCols, mat, nRows,
          vals, left, nRows, right, ldvt,
          work, lwork, iwork);
   delete[] work;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  else
   integer *iwork = new integer[8 * minMN];
   integer lwork = -1;
   integer err = 0;
   // workspace query:
   double optwork;
   integer nr = nRows, nc = nCols, ldvt_ = ldvt;
   dgesdd_(&jobz, &nr, &nc, mat, &nr,
          vals, left, &nr, right, &ldvt_,
          &optwork, &lwork, iwork, &err );
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = (integer)optwork;
   double* work = new double[lwork];
   // --- actual compute
   dgesdd_(&jobz, &nr, &nc, mat, &nr,
          vals, left, &nr, right, &ldvt_,
          work, &lwork, iwork, &err );
   delete[] work;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  endif   
   delete[] iwork;
}
// --- SVDREAL

// --- The next 248 lines were generated from math/snippets/SxBlasLib.cpp snippet SVDCOMPLEX
void singularValueDecomp (SxComplex8 *mat, int nRows, int nCols,
                          float *vals,
                          SxComplex8 *left,
                          SxComplex8 *right, // V^H
                          bool zeroSpace)
{
   char jobz = zeroSpace ? 'A'  // U= M x M, V=N x N
                         : 'S'; // U=M x s, V= s x N
   if (!left || !right)  {
      jobz = 'N';
      if (left || right) {
         // It is not possible to compute only left or only right vectors.
         std::cout << "Internal error in singular-value decomposition"
                   << std::endl;
         SX_EXIT;
      }
   }
   int minMN = nRows < nCols ? nRows : nCols;
#  if   defined (USE_VECLIB)
   SX_EXIT; // not yet implemented
#  elif defined (USE_ESSL)
   SX_EXIT; // not yet implemented
#  elif defined (USE_INTEL_MKL)
   int ldvt = zeroSpace ? nCols : minMN;
   int lrwork = (jobz == 'N') ? (7 * minMN) : ( (5 * minMN + 7) * minMN);
   int err = 0;
   int lwork = -1;
   float *rwork = new float[lrwork];
   int *iwork = new int[8 * minMN];
   // workspace query:
   MKL_Complex8 optwork;
   cgesdd(&jobz, &nRows, &nCols, (MKL_Complex8*)mat, &nRows,
          vals, (MKL_Complex8*)left, &nRows, (MKL_Complex8*)right, &ldvt,
          &optwork, &lwork, rwork, iwork, &err );
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = int(lround(optwork.real));
   MKL_Complex8* work = new MKL_Complex8[lwork];
   // --- actual compute
   cgesdd(&jobz, &nRows, &nCols, (MKL_Complex8*)mat, &nRows,
          vals, (MKL_Complex8*)left, &nRows, (MKL_Complex8*)right, &ldvt,
          work, &lwork, rwork, iwork, &err );
   delete[] work;
   delete[] iwork;
   delete[] rwork;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  elif defined (USE_ACML)
   SX_EXIT;
#  elif defined (USE_NETLIB)
   int ldvt = zeroSpace ? nCols : minMN;
   int lrwork = (jobz == 'N') ? (7 * minMN) : ( (5 * minMN + 7) * minMN);
   int err = 0;
   int lwork = -1;
   float *rwork = new float[lrwork];
   int *iwork = new int[8 * minMN];
   // workspace query:
   SxComplex8 optwork;
   err = LAPACKE_cgesdd_work(LAPACK_COL_MAJOR,
          jobz, nRows, nCols, mat, nRows,
          vals, left, nRows, right, ldvt,
          &optwork, lwork, rwork, iwork);
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = int(lround(optwork.re));
   SxComplex8 *work = new SxComplex8[lwork];
   // --- actual compute
   err = LAPACKE_cgesdd_work(LAPACK_COL_MAJOR,
          jobz, nRows, nCols, mat, nRows,
          vals, left, nRows, right, ldvt,
          work, lwork, rwork, iwork);
   delete[] work;
   delete[] iwork;
   delete[] rwork;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  else
   integer ldvt = zeroSpace ? nCols : minMN;
   integer lrwork = (jobz == 'N') ? (7 * minMN) : ( (5 * minMN + 7) * minMN);
   integer err = 0;
   integer lwork = -1;
   float *rwork = new float[lrwork];
   integer *iwork = new integer[8 * minMN];
   // workspace query:
   integer nr = nRows, nc = nCols;
   complex optwork;
   SX_EXIT;
   /* missing prototype
   cgesdd_(&jobz, &nr, &nc, (complex*)mat, &nr,
          vals, (complex*)left, &nr, (complex*)right, &ldvt,
          &optwork, &lwork, rwork, iwork, &err );
          */
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = (integer)optwork.r;
   complex* work = new complex[lwork];
   // --- actual compute
   SX_EXIT;
   /* missing prototype
   cgesdd_(&jobz, &nr, &nc, (complex*)mat, &nr,
          vals, (complex*)left, &nr, (complex*)right, &ldvt,
          work, &lwork, rwork, iwork, &err );
          */
   delete[] work;
   delete[] iwork;
   delete[] rwork;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  endif   

}

void singularValueDecomp (SxComplex16 *mat, int nRows, int nCols,
                          double *vals,
                          SxComplex16 *left,
                          SxComplex16 *right, // V^H
                          bool zeroSpace)
{
   char jobz = zeroSpace ? 'A'  // U= M x M, V=N x N
                         : 'S'; // U=M x s, V= s x N
   if (!left || !right)  {
      jobz = 'N';
      if (left || right) {
         // It is not possible to compute only left or only right vectors.
         std::cout << "Internal error in singular-value decomposition"
                   << std::endl;
         SX_EXIT;
      }
   }
   int minMN = nRows < nCols ? nRows : nCols;
#  if   defined (USE_VECLIB)
   SX_EXIT; // not yet implemented
#  elif defined (USE_ESSL)
   SX_EXIT; // not yet implemented
#  elif defined (USE_INTEL_MKL)
   int ldvt = zeroSpace ? nCols : minMN;
   int lrwork = (jobz == 'N') ? (7 * minMN) : ( (5 * minMN + 7) * minMN);
   int err = 0;
   int lwork = -1;
   double *rwork = new double[lrwork];
   int *iwork = new int[8 * minMN];
   // workspace query:
   MKL_Complex16 optwork;
   zgesdd(&jobz, &nRows, &nCols, (MKL_Complex16*)mat, &nRows,
          vals, (MKL_Complex16*)left, &nRows, (MKL_Complex16*)right, &ldvt,
          &optwork, &lwork, rwork, iwork, &err );
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = int(lround(optwork.real));
   MKL_Complex16* work = new MKL_Complex16[lwork];
   // --- actual compute
   zgesdd(&jobz, &nRows, &nCols, (MKL_Complex16*)mat, &nRows,
          vals, (MKL_Complex16*)left, &nRows, (MKL_Complex16*)right, &ldvt,
          work, &lwork, rwork, iwork, &err );
   delete[] work;
   delete[] iwork;
   delete[] rwork;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  elif defined (USE_ACML)
   SX_EXIT;
#  elif defined (USE_NETLIB)
   int ldvt = zeroSpace ? nCols : minMN;
   int lrwork = (jobz == 'N') ? (7 * minMN) : ( (5 * minMN + 7) * minMN);
   int err = 0;
   int lwork = -1;
   double *rwork = new double[lrwork];
   int *iwork = new int[8 * minMN];
   // workspace query:
   SxComplex16 optwork;
   err = LAPACKE_zgesdd_work(LAPACK_COL_MAJOR,
          jobz, nRows, nCols, mat, nRows,
          vals, left, nRows, right, ldvt,
          &optwork, lwork, rwork, iwork);
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = int(lround(optwork.re));
   SxComplex16 *work = new SxComplex16[lwork];
   // --- actual compute
   err = LAPACKE_zgesdd_work(LAPACK_COL_MAJOR,
          jobz, nRows, nCols, mat, nRows,
          vals, left, nRows, right, ldvt,
          work, lwork, rwork, iwork);
   delete[] work;
   delete[] iwork;
   delete[] rwork;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  else
   integer ldvt = zeroSpace ? nCols : minMN;
   integer lrwork = (jobz == 'N') ? (7 * minMN) : ( (5 * minMN + 7) * minMN);
   integer err = 0;
   integer lwork = -1;
   double *rwork = new double[lrwork];
   integer *iwork = new integer[8 * minMN];
   // workspace query:
   integer nr = nRows, nc = nCols;
   doublecomplex optwork;
   SX_EXIT;
   /* missing prototype
   zgesdd_(&jobz, &nr, &nc, (doublecomplex*)mat, &nr,
          vals, (doublecomplex*)left, &nr, (doublecomplex*)right, &ldvt,
          &optwork, &lwork, rwork, iwork, &err );
          */
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
   lwork = (integer)optwork.r;
   doublecomplex* work = new doublecomplex[lwork];
   // --- actual compute
   SX_EXIT;
   /* missing prototype
   zgesdd_(&jobz, &nr, &nc, (doublecomplex*)mat, &nr,
          vals, (doublecomplex*)left, &nr, (doublecomplex*)right, &ldvt,
          work, &lwork, rwork, iwork, &err );
          */
   delete[] work;
   delete[] iwork;
   delete[] rwork;
   if (err)  {
      std::cout << "svd err = " << err << std::endl;
      SX_EXIT;
   }
#  endif   

}
// --- SVDCOMPLEX

//------------------------------------------------------------------------------
// Matrix inversion
//------------------------------------------------------------------------------

// --- The next 140 lines were generated from snippets/SxBlasLib.cpp snippet INVREAL
void matInverse (float *mat, int nRows, int nCols)
{
   // --- (1) factorization
#  if   defined (USE_VECLIB)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      sgetrf (&nRows, &nCols, (float *)mat, &r, pivots, &err);
#  elif defined (USE_ESSL)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      sgetrf (nRows, nCols, mat, r, pivots, err);
#  elif defined (USE_INTEL_MKL)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      sgetrf (&nRows, &nCols, mat, &r, pivots, &err);
#  elif defined (USE_ACML)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      sgetrf (nRows, nCols, mat, r, pivots, &err);
#  elif defined (USE_NETLIB)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      err = LAPACKE_sgetrf (LAPACK_COL_MAJOR, nRows, nCols, mat, r, pivots);
#  else
      integer r = (integer)nRows, c = (integer)nCols, err = 0;
      integer *pivots = new integer[r < c ? r : c];
      sgetrf_ (&r, &c, (real *)mat, &r, pivots, &err);
#  endif
   if (err)  {   // TODO: throw execption
      std::cout << "matrix inverse: Error in sgetrf: " << err <<std::endl;
      delete [] pivots;
      SX_EXIT;
   }

   // --- (2) construct inverse
#  if   defined (USE_VECLIB)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      float *work = new float [lWork];
      sgetri (&r, (float *)mat, &r, pivots, work, &lWork, &err);
#  elif defined (USE_ESSL)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      float *work = new float [lWork];
      sgetri (r, mat, r, pivots, work, lWork, err);
#  elif defined (USE_INTEL_MKL)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      float *work = new float [lWork];
      sgetri (&r, mat, &r, pivots, work, &lWork, &err);
#  elif defined (USE_ACML)
      float *work = NULL; // ACML doesn't use work
      sgetri (r, mat, r, pivots, &err);
#  elif defined (USE_NETLIB)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      float *work = new float [lWork];
      err = LAPACKE_sgetri_work (LAPACK_COL_MAJOR, r, mat, r, pivots,
                                 work, lWork);
#  else
      integer lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      real *work = new real [lWork];
      sgetri_ (&r, (real *)mat, &r, pivots, work, &lWork, &err);
#  endif

   if ( err )  {  // TODO: throw execption
      std::cout << "matrix inverse: Error in sgetri: " << err <<std::endl;
      delete [] work; delete [] pivots;
      SX_EXIT;
   }

   delete [] work; delete [] pivots;
}

void matInverse (double *mat, int nRows, int nCols)
{
   // --- (1) factorization
#  if   defined (USE_VECLIB)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      dgetrf (&nRows, &nCols, (double *)mat, &r, pivots, &err);
#  elif defined (USE_ESSL)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      dgetrf (nRows, nCols, mat, r, pivots, err);
#  elif defined (USE_INTEL_MKL)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      dgetrf (&nRows, &nCols, mat, &r, pivots, &err);
#  elif defined (USE_ACML)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      dgetrf (nRows, nCols, mat, r, pivots, &err);
#  elif defined (USE_NETLIB)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      err = LAPACKE_dgetrf (LAPACK_COL_MAJOR, nRows, nCols, mat, r, pivots);
#  else
      integer r = (integer)nRows, c = (integer)nCols, err = 0;
      integer *pivots = new integer[r < c ? r : c];
      dgetrf_ (&r, &c, (doublereal *)mat, &r, pivots, &err);
#  endif
   if (err)  {   // TODO: throw execption
      std::cout << "matrix inverse: Error in dgetrf: " << err <<std::endl;
      delete [] pivots;
      SX_EXIT;
   }

   // --- (2) construct inverse
#  if   defined (USE_VECLIB)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      double *work = new double [lWork];
      dgetri (&r, (double *)mat, &r, pivots, work, &lWork, &err);
#  elif defined (USE_ESSL)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      double *work = new double [lWork];
      dgetri (r, mat, r, pivots, work, lWork, err);
#  elif defined (USE_INTEL_MKL)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      double *work = new double [lWork];
      dgetri (&r, mat, &r, pivots, work, &lWork, &err);
#  elif defined (USE_ACML)
      double *work = NULL; // ACML doesn't use work
      dgetri (r, mat, r, pivots, &err);
#  elif defined (USE_NETLIB)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      double *work = new double [lWork];
      err = LAPACKE_dgetri_work (LAPACK_COL_MAJOR, r, mat, r, pivots,
                                 work, lWork);
#  else
      integer lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      doublereal *work = new doublereal [lWork];
      dgetri_ (&r, (doublereal *)mat, &r, pivots, work, &lWork, &err);
#  endif

   if ( err )  {  // TODO: throw execption
      std::cout << "matrix inverse: Error in dgetri: " << err <<std::endl;
      delete [] work; delete [] pivots;
      SX_EXIT;
   }

   delete [] work; delete [] pivots;
}
// --- INVREAL

// --- The next 226 lines were generated from snippets/SxBlasLib.cpp snippet INVCOMPLEX
void matInverse (SxComplex8 *mat, int nRows, int nCols)
{
   // --- (1) factorization
#  if   defined (USE_VECLIB)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      cgetrf (&r, &c, (complex8_t *)mat, &r, pivots, &err);
#  elif defined (USE_ESSL)
      int r, c, err=0;
      // map complex inversion to real inversion
      //             +----+----+     +-----+-----+
      // +----+      | Re |-Im |     | re1 | im1 |
      // |cmlx|  =>  +----+----+  =  +-----+-----+
      // +----+      | Im | Re |     | im2 | re2 |
      //             +----+----+     +-----+-----+
      SxComplex8 *matPtr = mat;
      float *realMat = new float [ (2*nRows)*(2*nCols) ];
      float *re1Ptr  =  realMat;
      float *im1Ptr  = &realMat[nCols];
      float *im2Ptr  = &realMat[2*nRows*nCols];
      float *re2Ptr  = &realMat[2*nRows*nCols + nCols];

      for (r=0; r<nRows; r++, re1Ptr+=nCols, re2Ptr+=nCols,
                              im1Ptr+=nCols, im2Ptr+=nCols)
      {
         for (c=0; c<nCols; c++, matPtr++)  {
            *re1Ptr++ = *re2Ptr++ = matPtr->re;
            *im1Ptr++ =            -matPtr->im;
            *im2Ptr++ =             matPtr->im;
         }
      }
      // --- compute inverse of real helper matrix
      matInverse (realMat, 2*nRows, 2*nCols);
      // construct complex result
      // +----+----+     +-----+-----+
      // | Re | *  |     | re1 |  *  |    +----+
      // +----+----+  =  +-----+-----+ => |cmlx|
      // | Im | *  |     | im2 |  *  |    +----+
      // +----+----+     +-----+-----+
      matPtr = mat;
      re1Ptr =  realMat;
      im2Ptr = &realMat[2*nRows*nCols];

      for (r=0; r<nRows; r++, re1Ptr+=nCols, im2Ptr+=nCols)  {
         for (c=0; c<nCols; c++, matPtr++)  {
            matPtr->re = *re1Ptr++;
            matPtr->im = *im2Ptr++;
         }
      }
      delete [] realMat;
#  elif defined (USE_INTEL_MKL)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      cgetrf (&r, &c, (MKL_Complex8 *)mat, &r, pivots, &err);
#  elif defined (USE_ACML)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      cgetrf (r, c, (complex *)mat, r, pivots, &err);
#  elif defined (USE_NETLIB)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      err = LAPACKE_cgetrf (LAPACK_COL_MAJOR,
               nRows, nCols, mat, r, pivots);
#  else
      integer r = (integer)nRows, c = (integer)nCols, err = 0;
      integer *pivots = new integer[r < c ? r : c];
      cgetrf_ (&r, &c, (complex *)mat, &r, pivots, &err);
#  endif

#  ifndef USE_ESSL
      if (err)  {   // TODO: throw execption
         std::cout << "matrix inverse: Error in cgetrf: "<<err<<std::endl;
         delete [] pivots;
         SX_EXIT;
      }
#  endif

   // --- (2) construct inverse
#  if   defined (USE_VECLIB)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      complex8_t *work = new complex8_t [lWork];
      cgetri (&r, (complex8_t *)mat, &r, pivots, work, &lWork, &err);
#  elif defined (USE_ESSL)
      // empty
#  elif defined (USE_INTEL_MKL)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      MKL_Complex8 *work = new MKL_Complex8 [lWork];
      cgetri (&r, (MKL_Complex8 *)mat, &r, pivots, work, &lWork, &err);
#  elif defined (USE_ACML)
      complex *work = NULL;  // ACML doesn't use work
      cgetri (r, (complex *)mat, r, pivots, &err);
#  elif defined (USE_NETLIB)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      SxComplex8 *work = new SxComplex8 [lWork];
      err = LAPACKE_cgetri_work (LAPACK_COL_MAJOR, r, mat, r, pivots,
                                 work, lWork);
#  else
      integer lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      complex *work = new complex [lWork];
      cgetri_ (&r, (complex *)mat, &r, pivots, work, &lWork, &err);
#  endif

#  ifndef USE_ESSL
      if ( err )  {  // TODO: throw execption
         std::cout << "matrix inverse: Error in cgetri: "<<err<<std::endl;
         delete [] work; delete [] pivots;
         SX_EXIT;
      }

      delete [] work; delete [] pivots;
#  endif
}

void matInverse (SxComplex16 *mat, int nRows, int nCols)
{
   // --- (1) factorization
#  if   defined (USE_VECLIB)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      zgetrf (&r, &c, (complex16_t *)mat, &r, pivots, &err);
#  elif defined (USE_ESSL)
      int r, c, err=0;
      // map complex inversion to real inversion
      //             +----+----+     +-----+-----+
      // +----+      | Re |-Im |     | re1 | im1 |
      // |cmlx|  =>  +----+----+  =  +-----+-----+
      // +----+      | Im | Re |     | im2 | re2 |
      //             +----+----+     +-----+-----+
      SxComplex16 *matPtr = mat;
      double *realMat = new double [ (2*nRows)*(2*nCols) ];
      double *re1Ptr  =  realMat;
      double *im1Ptr  = &realMat[nCols];
      double *im2Ptr  = &realMat[2*nRows*nCols];
      double *re2Ptr  = &realMat[2*nRows*nCols + nCols];

      for (r=0; r<nRows; r++, re1Ptr+=nCols, re2Ptr+=nCols,
                              im1Ptr+=nCols, im2Ptr+=nCols)
      {
         for (c=0; c<nCols; c++, matPtr++)  {
            *re1Ptr++ = *re2Ptr++ = matPtr->re;
            *im1Ptr++ =            -matPtr->im;
            *im2Ptr++ =             matPtr->im;
         }
      }
      // --- compute inverse of real helper matrix
      matInverse (realMat, 2*nRows, 2*nCols);
      // construct complex result
      // +----+----+     +-----+-----+
      // | Re | *  |     | re1 |  *  |    +----+
      // +----+----+  =  +-----+-----+ => |cmlx|
      // | Im | *  |     | im2 |  *  |    +----+
      // +----+----+     +-----+-----+
      matPtr = mat;
      re1Ptr =  realMat;
      im2Ptr = &realMat[2*nRows*nCols];

      for (r=0; r<nRows; r++, re1Ptr+=nCols, im2Ptr+=nCols)  {
         for (c=0; c<nCols; c++, matPtr++)  {
            matPtr->re = *re1Ptr++;
            matPtr->im = *im2Ptr++;
         }
      }
      delete [] realMat;
#  elif defined (USE_INTEL_MKL)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      zgetrf (&r, &c, (MKL_Complex16 *)mat, &r, pivots, &err);
#  elif defined (USE_ACML)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      zgetrf (r, c, (doublecomplex *)mat, r, pivots, &err);
#  elif defined (USE_NETLIB)
      int r = nRows, c = nCols, err = 0;
      int *pivots = new int [r < c ? r : c];
      err = LAPACKE_zgetrf (LAPACK_COL_MAJOR,
               nRows, nCols, mat, r, pivots);
#  else
      integer r = (integer)nRows, c = (integer)nCols, err = 0;
      integer *pivots = new integer[r < c ? r : c];
      zgetrf_ (&r, &c, (doublecomplex *)mat, &r, pivots, &err);
#  endif

#  ifndef USE_ESSL
      if (err)  {   // TODO: throw execption
         std::cout << "matrix inverse: Error in zgetrf: "<<err<<std::endl;
         delete [] pivots;
         SX_EXIT;
      }
#  endif

   // --- (2) construct inverse
#  if   defined (USE_VECLIB)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      complex16_t *work = new complex16_t [lWork];
      zgetri (&r, (complex16_t *)mat, &r, pivots, work, &lWork, &err);
#  elif defined (USE_ESSL)
      // empty
#  elif defined (USE_INTEL_MKL)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      MKL_Complex16 *work = new MKL_Complex16 [lWork];
      zgetri (&r, (MKL_Complex16 *)mat, &r, pivots, work, &lWork, &err);
#  elif defined (USE_ACML)
      doublecomplex *work = NULL;  // ACML doesn't use work
      zgetri (r, (doublecomplex *)mat, r, pivots, &err);
#  elif defined (USE_NETLIB)
      int lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      SxComplex16 *work = new SxComplex16 [lWork];
      err = LAPACKE_zgetri_work (LAPACK_COL_MAJOR, r, mat, r, pivots,
                                 work, lWork);
#  else
      integer lWork = r;   // for opt. performance = r * OPT_BLOCKSIZE (from ILAENV)
      doublecomplex *work = new doublecomplex [lWork];
      zgetri_ (&r, (doublecomplex *)mat, &r, pivots, work, &lWork, &err);
#  endif

#  ifndef USE_ESSL
      if ( err )  {  // TODO: throw execption
         std::cout << "matrix inverse: Error in zgetri: "<<err<<std::endl;
         delete [] work; delete [] pivots;
         SX_EXIT;
      }

      delete [] work; delete [] pivots;
#  endif
}
// --- INVCOMPLEX

void matInverseTri (float * /*mat*/, int /*nRows*/, enum UPLO /*uplo*/)
{
   SX_EXIT; // not yet implemented
}

void matInverseTri (double *mat, int nRows, enum UPLO uplo)
{
   char uploChar = (uplo == UpperRight ? 'U' : 'L');
   // --- (1) Factorization: A = U*D*U^t 
#  if   defined (USE_VECLIB)
      SX_EXIT; // not yet implemented
#  elif defined ( USE_ESSL)
      SX_EXIT; // not yet implemented
#  elif defined (USE_INTEL_MKL)
      int r = nRows, err = 0;
      int *pivots = new int[r];
      dsptrf (&uploChar, &r, mat, pivots, &err);
#  elif defined (USE_ACML)
      int r = nRows, err = 0;
      int *pivots = new int[r];
      dsptrf (uploChar, r, mat, pivots, &err);
#  elif defined (USE_NETLIB)
      int r = nRows, err = 0;
      int *pivots = new int[r];
      err = LAPACKE_dsptrf (LAPACK_COL_MAJOR, uploChar, r, mat, pivots);
#  else
      integer r = (integer)nRows, err = 0;
      integer *pivots = new integer[r];
      dsptrf_ (&uploChar, &r, (doublereal *)mat, pivots, &err);
#  endif    
      if ( err )  {  // TODO: throw execption
         std::cout << "packed matrix inverse: Error in DSPTRF: "
                   << err << std::endl;
         SX_EXIT;
      }

   // --- (2) Inverse of A
#  if   defined (USE_VECLIB)
       SX_EXIT; // not yet implemented
#  elif defined (USE_ESSL)
       SX_EXIT; // not yet implemented
#  elif defined (USE_INTEL_MKL)
      double *work = new double [r];
      dsptri (&uploChar, &r, mat, pivots, work, &err);
#  elif defined (USE_ACML)
      double *work = NULL; // ACML doesn't use work
      dsptri (uploChar, r, mat, pivots, &err);
#  elif defined (USE_NETLIB)
      double *work = new double [r];
      err = LAPACKE_dsptri_work (LAPACK_COL_MAJOR, uploChar, r, mat, pivots, work);
#  else
      doublereal *work = new doublereal [r];
      dsptri_ (&uploChar, &r, (doublereal *)mat, pivots, work, &err);
#  endif      

#  ifndef USE_ESSL      
      if ( err )  {  // TODO: throw execption
         std::cout << "packed matrix inverse: Error in DSPTRI: "
                   << err << std::endl;
         delete [] work; delete [] pivots;
         SX_EXIT;
      }

      delete [] work; delete [] pivots;
#  endif
}

void matInverseTri (SxComplex8 * /*mat*/, int /*nRows*/, enum UPLO /*uplo*/)
{
   SX_EXIT; // not yet implemented
}

void matInverseTri (SxComplex16 * /*mat*/, int /*nRows*/, enum UPLO /*uplo*/)
{
   SX_EXIT; // not yet implemented
}

//------------------------------------------------------------------------------
// Linear equation solver (least square based)
//------------------------------------------------------------------------------
// --- The next 318 lines were generated from math/snippets/SxBlasLib.cpp snippet SOLVELIN
void solveLinEq (float *mat, int nRows, int nCols, float *b, int bCols)
{
#  if   defined (USE_VECLIB)
      //TODO
      SX_EXIT;
#  elif defined (USE_ESSL)
      int iopt = 2; //compute singulare Values, V and U^TB
      int m = nRows;
      int n = nCols;
      int nrhs = bCols;
      int lda = m;
      int ldb = m;
      char transa = 'N';
      int naux = 0; //ESSL choose size of work array dynamically
      float *aux = NULL; //workarray is ignored
      int info;
      float *s = new float [n];
      float *x = new float [n*nrhs];
      float tau = 1e-10; // error tolerance for zero
      sgesvf(iopt,mat,lda,b,ldb,nrhs,s,m,n,aux,naux);
      sgesvs (mat,n,b,n,nrhs,s,x,n,m,n,tau);
      // results are now in x; copy to b
      for(int i = 0; i < n*nrhs; i++)   {
         b[i] = x[i];
      }

      delete [] s; delete [] x;

      //TODO TEST
      SX_EXIT;
#  elif defined (USE_INTEL_MKL)
      // for ilaenv
      MKL_INT iSpec = 9;
      char *name    = const_cast<char*>("SGELSD");
      char *opts    = const_cast<char*>(" ");
      MKL_INT n1    = 0;
      MKL_INT n2    = 0;
      MKL_INT n3    = 0;
      MKL_INT n4    = 0;
      // for sgelsd
      MKL_INT m          = nRows;
      MKL_INT n          = nCols;
      MKL_INT minmn      = n < m ? n : m;
      MKL_INT nrhs       = bCols;
      MKL_INT lda        = m;
      MKL_INT ldb        = m;
      float *s          = new float [minmn];
      float rcond       = -1.;
      MKL_INT rank;
      MKL_INT SMLSIZ     = ilaenv(&iSpec,name,opts,&n1,&n2,&n3,&n4);
      MKL_INT logVal     = MKL_INT(log( 1.0*n/(SMLSIZ+1))/log(2.0)) + 1;
      MKL_INT NLVL       = 0 < logVal ? logVal : 0;
      MKL_INT lwork      = -1;
      float autolwork;
      MKL_INT liwork     = 3 * minmn * NLVL + 11 * minmn;
      MKL_INT *iwork     = new MKL_INT [liwork];
      MKL_INT info;
      // determine optimal lwork
      sgelsd (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,&autolwork,&lwork,iwork,&info);
      lwork              = MKL_INT(autolwork+0.5);
      float *work       = new float [lwork];

      sgelsd (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,work,&lwork,iwork,&info);

      delete [] s; delete [] work; delete [] iwork;
#  elif defined (USE_ACML)
      int m        = nRows;
      int n        = nCols;
      int nrhs     = bCols;
      int lda      = m;
      int ldb      = m;
      int minmn    = n < m ? n : m;
      float *s    = new float [minmn];
      float rcond = -1.;
      int rank;
      int info;
      sgelsd (m,n,nrhs,mat,lda,b,ldb,s,rcond,&rank,&info);

      delete [] s;
      //TODO TEST
      SX_EXIT;
#   elif defined(USE_NETLIB)
      // for sgelsd
      int m          = nRows;
      int n          = nCols;
      int minmn      = n < m ? n : m;
      int nrhs       = bCols;
      int lda        = m;
      int ldb        = m;
      float *s    = new float [minmn];
      float rcond = -1.;
      int rank;
      int lwork = -1;
      float autolwork;
      int liwork = -1;
      int info;
      // determine optimal lwork
      info = LAPACKE_sgelsd_work (LAPACK_COL_MAJOR,
               m,n,nrhs,mat,lda,b,ldb,s,rcond,&rank,&autolwork,lwork,&liwork);
      int *iwork = new int [liwork];
      lwork = int(autolwork+0.5);
      float *work = new float [lwork];

      info = LAPACKE_sgelsd_work (LAPACK_COL_MAJOR,
            m,n,nrhs,mat,lda,b,ldb,s,rcond,&rank,work,lwork,iwork);

      delete [] s; delete [] work; delete [] iwork;
      if (info)  {
         std::cout << "solveLinEq: Error in sgelsd: " << info << std::endl;
         SX_EXIT;
      }
#  else
      // for ilaenv
      integer iSpec = 9;
      char *name    = const_cast<char*>("SGELSD");
      char *opts    = const_cast<char*>(" ");
      integer n1    = 0;
      integer n2    = 0;
      integer n3    = 0;
      integer n4    = 0;
      ftnlen lname  = 6;
      ftnlen lopts  = 1;
      // for sgelsd
      integer m          = nRows;
      integer n          = nCols;
      integer minmn      = n < m ? n : m;
      integer nrhs       = bCols;
      integer lda        = m;
      integer ldb        = m;
      float *s          = new float [minmn];
      float rcond       = -1.;
      integer rank;
#ifdef MACOSX
#  if ( DIST_VERSION_L >= 1070L )
      integer SMLSIZ     = ilaenv_(&iSpec,name,opts,&n1,&n2,&n3,&n4);
#  else
      integer SMLSIZ     = ilaenv_(&iSpec,name,opts,&n1,&n2,&n3,&n4,lname,lopts);
#  endif /* DIST_VERSION_L */
#else
      integer SMLSIZ     = ilaenv_(&iSpec,name,opts,&n1,&n2,&n3,&n4,lname,lopts);
#endif /* MACOSX */
      integer logVal     = integer(log( double(n)/double(SMLSIZ+1))/log(2.0)) + 1;
      integer NLVL       = 0 < logVal ? logVal : 0;
      integer lwork      = -1;
      float autolwork;
      integer liwork     = 3 * minmn * NLVL + 11 * minmn;
      integer *iwork     = new integer [liwork];
      integer info;
      // determine optimal lwork
      sgelsd_ (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,&autolwork,&lwork,iwork,&info);
      lwork              = integer(autolwork+0.5);
      float *work       = new float [lwork];

      sgelsd_ (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,work,&lwork,iwork,&info);

      delete [] s; delete [] work; delete [] iwork;
#  endif
}

void solveLinEq (double *mat, int nRows, int nCols, double *b, int bCols)
{
#  if   defined (USE_VECLIB)
      //TODO
      SX_EXIT;
#  elif defined (USE_ESSL)
      int iopt = 2; //compute singulare Values, V and U^TB
      int m = nRows;
      int n = nCols;
      int nrhs = bCols;
      int lda = m;
      int ldb = m;
      char transa = 'N';
      int naux = 0; //ESSL choose size of work array dynamically
      double *aux = NULL; //workarray is ignored
      int info;
      double *s = new double [n];
      double *x = new double [n*nrhs];
      double tau = 1e-10; // error tolerance for zero
      dgesvf(iopt,mat,lda,b,ldb,nrhs,s,m,n,aux,naux);
      dgesvs (mat,n,b,n,nrhs,s,x,n,m,n,tau);
      // results are now in x; copy to b
      for(int i = 0; i < n*nrhs; i++)   {
         b[i] = x[i];
      }

      delete [] s; delete [] x;

      //TODO TEST
      SX_EXIT;
#  elif defined (USE_INTEL_MKL)
      // for ilaenv
      MKL_INT iSpec = 9;
      char *name    = const_cast<char*>("DGELSD");
      char *opts    = const_cast<char*>(" ");
      MKL_INT n1    = 0;
      MKL_INT n2    = 0;
      MKL_INT n3    = 0;
      MKL_INT n4    = 0;
      // for dgelsd
      MKL_INT m          = nRows;
      MKL_INT n          = nCols;
      MKL_INT minmn      = n < m ? n : m;
      MKL_INT nrhs       = bCols;
      MKL_INT lda        = m;
      MKL_INT ldb        = m;
      double *s          = new double [minmn];
      double rcond       = -1.;
      MKL_INT rank;
      MKL_INT SMLSIZ     = ilaenv(&iSpec,name,opts,&n1,&n2,&n3,&n4);
      MKL_INT logVal     = MKL_INT(log( 1.0*n/(SMLSIZ+1))/log(2.0)) + 1;
      MKL_INT NLVL       = 0 < logVal ? logVal : 0;
      MKL_INT lwork      = -1;
      double autolwork;
      MKL_INT liwork     = 3 * minmn * NLVL + 11 * minmn;
      MKL_INT *iwork     = new MKL_INT [liwork];
      MKL_INT info;
      // determine optimal lwork
      dgelsd (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,&autolwork,&lwork,iwork,&info);
      lwork              = MKL_INT(autolwork+0.5);
      double *work       = new double [lwork];

      dgelsd (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,work,&lwork,iwork,&info);

      delete [] s; delete [] work; delete [] iwork;
#  elif defined (USE_ACML)
      int m        = nRows;
      int n        = nCols;
      int nrhs     = bCols;
      int lda      = m;
      int ldb      = m;
      int minmn    = n < m ? n : m;
      double *s    = new double [minmn];
      double rcond = -1.;
      int rank;
      int info;
      dgelsd (m,n,nrhs,mat,lda,b,ldb,s,rcond,&rank,&info);

      delete [] s;
      //TODO TEST
      SX_EXIT;
#   elif defined(USE_NETLIB)
      // for dgelsd
      int m          = nRows;
      int n          = nCols;
      int minmn      = n < m ? n : m;
      int nrhs       = bCols;
      int lda        = m;
      int ldb        = m;
      double *s    = new double [minmn];
      double rcond = -1.;
      int rank;
      int lwork = -1;
      double autolwork;
      int liwork = -1;
      int info;
      // determine optimal lwork
      info = LAPACKE_dgelsd_work (LAPACK_COL_MAJOR,
               m,n,nrhs,mat,lda,b,ldb,s,rcond,&rank,&autolwork,lwork,&liwork);
      int *iwork = new int [liwork];
      lwork = int(autolwork+0.5);
      double *work = new double [lwork];

      info = LAPACKE_dgelsd_work (LAPACK_COL_MAJOR,
            m,n,nrhs,mat,lda,b,ldb,s,rcond,&rank,work,lwork,iwork);

      delete [] s; delete [] work; delete [] iwork;
      if (info)  {
         std::cout << "solveLinEq: Error in dgelsd: " << info << std::endl;
         SX_EXIT;
      }
#  else
      // for ilaenv
      integer iSpec = 9;
      char *name    = const_cast<char*>("DGELSD");
      char *opts    = const_cast<char*>(" ");
      integer n1    = 0;
      integer n2    = 0;
      integer n3    = 0;
      integer n4    = 0;
      ftnlen lname  = 6;
      ftnlen lopts  = 1;
      // for dgelsd
      integer m          = nRows;
      integer n          = nCols;
      integer minmn      = n < m ? n : m;
      integer nrhs       = bCols;
      integer lda        = m;
      integer ldb        = m;
      double *s          = new double [minmn];
      double rcond       = -1.;
      integer rank;
#ifdef MACOSX
#  if ( DIST_VERSION_L >= 1070L )
      integer SMLSIZ     = ilaenv_(&iSpec,name,opts,&n1,&n2,&n3,&n4);
#  else
      integer SMLSIZ     = ilaenv_(&iSpec,name,opts,&n1,&n2,&n3,&n4,lname,lopts);
#  endif /* DIST_VERSION_L */
#else
      integer SMLSIZ     = ilaenv_(&iSpec,name,opts,&n1,&n2,&n3,&n4,lname,lopts);
#endif /* MACOSX */
      integer logVal     = integer(log( double(n)/double(SMLSIZ+1))/log(2.0)) + 1;
      integer NLVL       = 0 < logVal ? logVal : 0;
      integer lwork      = -1;
      double autolwork;
      integer liwork     = 3 * minmn * NLVL + 11 * minmn;
      integer *iwork     = new integer [liwork];
      integer info;
      // determine optimal lwork
      dgelsd_ (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,&autolwork,&lwork,iwork,&info);
      lwork              = integer(autolwork+0.5);
      double *work       = new double [lwork];

      dgelsd_ (&m,&n,&nrhs,mat,&lda,b,&ldb,s,&rcond,&rank,work,&lwork,iwork,&info);

      delete [] s; delete [] work; delete [] iwork;
#  endif
}
// --- SOLVELIN

void solveLinEq (SxComplex8 * /*mat*/, int /*nRows*/, int /*nCols*/,
                 SxComplex8 * /*b*/, int /*bCols*/)
{
#  if   defined (USE_VECLIB)
      //TODO
      SX_EXIT;
#  elif defined (USE_ESSL)
      //TODO 
      SX_EXIT;
#  elif defined (USE_INTEL_MKL)
      //TODO 
      SX_EXIT;
#  elif defined (USE_ACML)
      //TODO
      SX_EXIT;
#  endif      
}

void solveLinEq (SxComplex16 * /*mat*/, int /*nRows*/, int /*nCols*/,
                 SxComplex16 * /*b*/, int /*bCols*/)
{
#  if   defined (USE_VECLIB)
      //TODO
      SX_EXIT;
#  elif defined (USE_ESSL)
      //TODO 
      SX_EXIT;
#  elif defined (USE_INTEL_MKL)
      //TODO 
      SX_EXIT;
#  elif defined (USE_ACML)
      //TODO 
      SX_EXIT;
#  else 
      //TODO
      SX_EXIT;
#  endif      
}

//------------------------------------------------------------------------------
// Eigensolver
//------------------------------------------------------------------------------
// --- The next 370 lines were generated from snippets/SxBlasLib.cpp snippet EIGENREAL
int matEigensolver (SxComplex<float> *eigVals, float *eigVecs,
                    float *inMat, int n, EIGCMD cmd, int size)
{
#  if defined (USE_VECLIB) || defined (USE_INTEL_MKL) || defined (USE_NETLIB)
      int rank  = n, info = 0, ldVecLeft = 1;
      int lWork = !size ? 4*n : size;
      int workDim = lWork;
#  elif defined (USE_ACML)
      int rank  = n, info = 0, ldVecLeft = 1, lWork=size;
#  else
      integer rank  = (integer)n, info = 0, ldVecLeft = 1;
      integer lWork = !size ? 4*n : size;
      integer workDim = lWork;
#  endif
   if ( cmd == OptSize )  lWork = -1;
   char jobvl = 'N', jobvr = 'V';

#  if   defined (USE_VECLIB)
      float *work = new float[workDim];
      float *epsRe = new float[n], *epsIm = new float[n];
      float *eigVecLeft = new float[ldVecLeft];
      sgeev (&jobvl, &jobvr,
             &rank, (float *)inMat, &rank, epsRe, epsIm, eigVecLeft,
             &ldVecLeft, (float *)eigVecs, &rank, work, &lWork, &info, 0, 0);
#  elif defined (USE_ESSL)
      complex<float> *vecs = new complex<float>[rank*rank];
      int iOpt = 1;   // compute both vecs and vals
      if (cmd != OptSize)  {  // not supported by ESSL
         if (cmd == ValuesOnly)  iOpt = 0;
         sgeev (iOpt, inMat, rank,
                (complex<float> *)eigVals, vecs, rank,
                NULL, rank, NULL, 0);
         complex<float> *srcPtr = vecs;
         float *dstPtr = eigVecs;
         int i, len = n*n;
         for (i=0; i < len; i++, srcPtr++)  {
            *dstPtr++ = real(*srcPtr);
         }
         // --- normalize eigenvectors
         float c;
         for (i=0; i < n; i++)  {
            c = 1. / norm2 (&eigVecs[i*n], n);
            scale (&eigVecs[i*n], c, n);
         }
      }
      delete [] vecs;
      return 0;
#  elif defined (USE_INTEL_MKL)
      float *work = new float[workDim];
      float *epsRe = new float[n], *epsIm = new float[n];
      float *eigVecLeft = new float[ldVecLeft];
      sgeev (&jobvl, &jobvr,
             &rank, inMat, &rank, epsRe, epsIm, eigVecLeft,
             &ldVecLeft, eigVecs, &rank, work, &lWork, &info);
#  elif defined (USE_ACML)
      float *work = NULL; // ACML doesn't use work
      float *epsRe = new float [n], *epsIm = new float [n];
      float *eigVecLeft = new float [ldVecLeft];
      sgeev (jobvl, jobvr,
             rank, inMat, rank, epsRe, epsIm, eigVecLeft,
             ldVecLeft, eigVecs, rank, &info);
#  elif defined (USE_NETLIB)
      float *work = new float[workDim];
      float *epsRe = new float[n], *epsIm = new float[n];
      float *eigVecLeft = new float[ldVecLeft];
      info = LAPACKE_sgeev_work (LAPACK_COL_MAJOR,
             jobvl, jobvr,
             rank, inMat, rank, epsRe, epsIm, eigVecLeft,
             ldVecLeft, eigVecs, rank, work, lWork);
#  else
      real *work = new real [workDim];
      real *epsRe = new real [n], *epsIm = new real [n];
      real *eigVecLeft = new real [ldVecLeft];
      sgeev_ (&jobvl, &jobvr,
              &rank, (real *)inMat, &rank, epsRe, epsIm, eigVecLeft,
              &ldVecLeft, (real *)eigVecs, &rank, work, &lWork, &info);
#  endif

#  ifndef USE_ESSL
      float *ptr=(float *)eigVals;
      float *rePtr=(float *)epsRe, *imPtr=(float *)epsIm;
      for (int i=0; i < n; i++)  {
         *ptr++ = *rePtr++;  // eigVals[i].re = epsRe[i];
         *ptr++ = *imPtr++;  // eigVals[i].im = epsIm[i];
      }
      delete [] eigVecLeft; delete [] epsIm; delete [] epsRe; delete [] work;

      if ( info )  {  // TODO: throw exception
         std::cout << "matEigensolver: Error in sgeev: " << info << std::endl;
         SX_EXIT;
      }
      return (cmd == OptSize)  ?  (int)work[0] : 0;
#  endif

}

void matSymEigensolver (float *eigVals, float *eigVecs,
                        float *inMat, int n, EIGCMD cmd)
{
   char jobz = (cmd == ValuesOnly) ? 'N' : 'V';
   char range = 'A'; // all eigenvalues
   char uplo = 'U';  // use upper right half
   float vdummy = 0., abstol = 0.;
   // abstol = dlamch ('Safe minimum');
   // use automatic workspace LAPACKE interface where available
#  if defined (USE_NETLIB) || defined (USE_ESSL) || defined (USE_ACML) \
   || defined (USE_INTEL_MKL)
#     define SX_USE_LAPACKE
#  endif

   // --- workspace query
#  ifdef SX_USE_LAPACKE
      lapack_int idummy = 0;
      lapack_int nFound = 0;
      lapack_int info = 0;
# else
      int idummy = 0;
      int nFound = 0;
      int info = 0;
      int liwork = -1;
      int queryWork = -1;
      float lwork_ = -1.;
#     if defined (USE_VECLIB)
         SX_EXIT;
#     else
         ssyevr_ (&jobz, &range, &uplo,
             &n, (real *)inMat, &n,
             (real *)&vdummy, (real *)&vdummy,
             &idummy, &idummy, (real *)&abstol,
             &nFound, (real *)eigVals, (real *)eigVecs, &n,
             NULL, (real *)&lwork_, &queryWork, &liwork, &queryWork,
             &info);
#     endif
      if ( info )  {  // TODO: throw exception
         std::cout << "matSymEigensolver: Error in ssyevr: " << info << std::endl;
         SX_EXIT;
      }
      // --- allocate workspaces
      int *iwork = new int[liwork];
      int lwork = (int)lround(lwork_);
      float *work = new float[lwork];
#  endif

   // --- computational call
#  if   defined (USE_VECLIB)
      SX_EXIT;
#  elif defined (SX_USE_LAPACKE)
      // allocate iSuppZ workspace
      lapack_int *iSuppZ = new lapack_int[2 * n];
      info = LAPACKE_ssyevr (LAPACK_COL_MAJOR,
             jobz, range, uplo,
             n, inMat, n,
             vdummy, vdummy, idummy, idummy, abstol,
             &nFound, eigVals, eigVecs, n,
             iSuppZ);
#  else
      // allocate iSuppZ workspace
      int *iSuppZ = new int[2 * n];
      ssyevr_ (&jobz, &range, &uplo,
             &n, (real *)inMat, &n,
             (real *)&vdummy, (real *)&vdummy,
             &idummy, &idummy, (real *)&abstol,
             &nFound, (real *)eigVals, (real *)eigVecs, &n,
             iSuppZ, (real *)work, &lwork, iwork, &liwork,
             &info);
#  endif

   // --- delete workspaces
   delete [] iSuppZ;
#  ifndef SX_USE_LAPACKE
      delete [] work;
      delete [] iwork;
#  endif

   // --- check for success
   if ( info )  {  // TODO: throw exception
      std::cout << "matSymEigensolver: Error in ssyevr: "
                << info << std::endl;
      SX_EXIT;
   }
#  ifdef SX_USE_LAPACKE
#     undef SX_USE_LAPACKE
#  endif
}

int matEigensolver (SxComplex<double> *eigVals, double *eigVecs,
                    double *inMat, int n, EIGCMD cmd, int size)
{
#  if defined (USE_VECLIB) || defined (USE_INTEL_MKL) || defined (USE_NETLIB)
      int rank  = n, info = 0, ldVecLeft = 1;
      int lWork = !size ? 4*n : size;
      int workDim = lWork;
#  elif defined (USE_ACML)
      int rank  = n, info = 0, ldVecLeft = 1, lWork=size;
#  else
      integer rank  = (integer)n, info = 0, ldVecLeft = 1;
      integer lWork = !size ? 4*n : size;
      integer workDim = lWork;
#  endif
   if ( cmd == OptSize )  lWork = -1;
   char jobvl = 'N', jobvr = 'V';

#  if   defined (USE_VECLIB)
      double *work = new double[workDim];
      double *epsRe = new double[n], *epsIm = new double[n];
      double *eigVecLeft = new double[ldVecLeft];
      dgeev (&jobvl, &jobvr,
             &rank, (double *)inMat, &rank, epsRe, epsIm, eigVecLeft,
             &ldVecLeft, (double *)eigVecs, &rank, work, &lWork, &info, 0, 0);
#  elif defined (USE_ESSL)
      complex<double> *vecs = new complex<double>[rank*rank];
      int iOpt = 1;   // compute both vecs and vals
      if (cmd != OptSize)  {  // not supported by ESSL
         if (cmd == ValuesOnly)  iOpt = 0;
         dgeev (iOpt, inMat, rank,
                (complex<double> *)eigVals, vecs, rank,
                NULL, rank, NULL, 0);
         complex<double> *srcPtr = vecs;
         double *dstPtr = eigVecs;
         int i, len = n*n;
         for (i=0; i < len; i++, srcPtr++)  {
            *dstPtr++ = real(*srcPtr);
         }
         // --- normalize eigenvectors
         double c;
         for (i=0; i < n; i++)  {
            c = 1. / norm2 (&eigVecs[i*n], n);
            scale (&eigVecs[i*n], c, n);
         }
      }
      delete [] vecs;
      return 0;
#  elif defined (USE_INTEL_MKL)
      double *work = new double[workDim];
      double *epsRe = new double[n], *epsIm = new double[n];
      double *eigVecLeft = new double[ldVecLeft];
      dgeev (&jobvl, &jobvr,
             &rank, inMat, &rank, epsRe, epsIm, eigVecLeft,
             &ldVecLeft, eigVecs, &rank, work, &lWork, &info);
#  elif defined (USE_ACML)
      double *work = NULL; // ACML doesn't use work
      double *epsRe = new double [n], *epsIm = new double [n];
      double *eigVecLeft = new double [ldVecLeft];
      dgeev (jobvl, jobvr,
             rank, inMat, rank, epsRe, epsIm, eigVecLeft,
             ldVecLeft, eigVecs, rank, &info);
#  elif defined (USE_NETLIB)
      double *work = new double[workDim];
      double *epsRe = new double[n], *epsIm = new double[n];
      double *eigVecLeft = new double[ldVecLeft];
      info = LAPACKE_dgeev_work (LAPACK_COL_MAJOR,
             jobvl, jobvr,
             rank, inMat, rank, epsRe, epsIm, eigVecLeft,
             ldVecLeft, eigVecs, rank, work, lWork);
#  else
      doublereal *work = new doublereal [workDim];
      doublereal *epsRe = new doublereal [n], *epsIm = new doublereal [n];
      doublereal *eigVecLeft = new doublereal [ldVecLeft];
      dgeev_ (&jobvl, &jobvr,
              &rank, (doublereal *)inMat, &rank, epsRe, epsIm, eigVecLeft,
              &ldVecLeft, (doublereal *)eigVecs, &rank, work, &lWork, &info);
#  endif

#  ifndef USE_ESSL
      double *ptr=(double *)eigVals;
      double *rePtr=(double *)epsRe, *imPtr=(double *)epsIm;
      for (int i=0; i < n; i++)  {
         *ptr++ = *rePtr++;  // eigVals[i].re = epsRe[i];
         *ptr++ = *imPtr++;  // eigVals[i].im = epsIm[i];
      }
      delete [] eigVecLeft; delete [] epsIm; delete [] epsRe; delete [] work;

      if ( info )  {  // TODO: throw exception
         std::cout << "matEigensolver: Error in dgeev: " << info << std::endl;
         SX_EXIT;
      }
      return (cmd == OptSize)  ?  (int)work[0] : 0;
#  endif

}

void matSymEigensolver (double *eigVals, double *eigVecs,
                        double *inMat, int n, EIGCMD cmd)
{
   char jobz = (cmd == ValuesOnly) ? 'N' : 'V';
   char range = 'A'; // all eigenvalues
   char uplo = 'U';  // use upper right half
   double vdummy = 0., abstol = 0.;
   // abstol = dlamch ('Safe minimum');
   // use automatic workspace LAPACKE interface where available
#  if defined (USE_NETLIB) || defined (USE_ESSL) || defined (USE_ACML) \
   || defined (USE_INTEL_MKL)
#     define SX_USE_LAPACKE
#  endif

   // --- workspace query
#  ifdef SX_USE_LAPACKE
      lapack_int idummy = 0;
      lapack_int nFound = 0;
      lapack_int info = 0;
# else
      int idummy = 0;
      int nFound = 0;
      int info = 0;
      int liwork = -1;
      int queryWork = -1;
      double lwork_ = -1.;
#     if defined (USE_VECLIB)
         SX_EXIT;
#     else
         dsyevr_ (&jobz, &range, &uplo,
             &n, (doublereal *)inMat, &n,
             (doublereal *)&vdummy, (doublereal *)&vdummy,
             &idummy, &idummy, (doublereal *)&abstol,
             &nFound, (doublereal *)eigVals, (doublereal *)eigVecs, &n,
             NULL, (doublereal *)&lwork_, &queryWork, &liwork, &queryWork,
             &info);
#     endif
      if ( info )  {  // TODO: throw exception
         std::cout << "matSymEigensolver: Error in dsyevr: " << info << std::endl;
         SX_EXIT;
      }
      // --- allocate workspaces
      int *iwork = new int[liwork];
      int lwork = (int)lround(lwork_);
      double *work = new double[lwork];
#  endif

   // --- computational call
#  if   defined (USE_VECLIB)
      SX_EXIT;
#  elif defined (SX_USE_LAPACKE)
      // allocate iSuppZ workspace
      lapack_int *iSuppZ = new lapack_int[2 * n];
      info = LAPACKE_dsyevr (LAPACK_COL_MAJOR,
             jobz, range, uplo,
             n, inMat, n,
             vdummy, vdummy, idummy, idummy, abstol,
             &nFound, eigVals, eigVecs, n,
             iSuppZ);
#  else
      // allocate iSuppZ workspace
      int *iSuppZ = new int[2 * n];
      dsyevr_ (&jobz, &range, &uplo,
             &n, (doublereal *)inMat, &n,
             (doublereal *)&vdummy, (doublereal *)&vdummy,
             &idummy, &idummy, (doublereal *)&abstol,
             &nFound, (doublereal *)eigVals, (doublereal *)eigVecs, &n,
             iSuppZ, (doublereal *)work, &lwork, iwork, &liwork,
             &info);
#  endif

   // --- delete workspaces
   delete [] iSuppZ;
#  ifndef SX_USE_LAPACKE
      delete [] work;
      delete [] iwork;
#  endif

   // --- check for success
   if ( info )  {  // TODO: throw exception
      std::cout << "matSymEigensolver: Error in dsyevr: "
                << info << std::endl;
      SX_EXIT;
   }
#  ifdef SX_USE_LAPACKE
#     undef SX_USE_LAPACKE
#  endif
}
// --- EIGENREAL

// --- The next 408 lines were generated from snippets/SxBlasLib.cpp snippet EIGENCOMPLEX
int matEigensolver (SxComplex8 *eigVals, SxComplex8 *eigVecs,
                    SxComplex8 *inMat, int n, EIGCMD cmd, int size)
{
#  if defined (USE_VECLIB) || defined (USE_INTEL_MKL) || defined (USE_NETLIB)
      int rank  = n, info = 0, ldVecLeft = 1;
      int lWork = !size ? 2*n : size;
      int workDim = lWork;
#  elif defined (USE_ESSL)
      int rank  = n, ldVecLeft = 1;
      int lWork = !size ? 2*n : size;
      int workDim = lWork;
#  elif defined (USE_ACML)
      int rank  = n, info = 0, ldVecLeft = 1, lWork=size;
#  else
      integer rank  = (integer)n, info = 0, ldVecLeft = 1;
      integer lWork = !size ? 2*n : size;
      integer workDim = lWork;
#  endif
   if ( cmd == OptSize )  lWork = -1;
   char jobvl = 'N', jobvr = 'V';

#  if   defined (USE_VECLIB)
      complex8_t *eigVecLeft = new complex8_t [ldVecLeft];
      complex8_t *work       = new complex8_t [workDim];
      float      *rWork      = new float      [workDim];
      cgeev (&jobvl, &jobvr,
             &rank, (complex8_t *)inMat, &rank,
             (complex8_t *)eigVals, eigVecLeft,
             &ldVecLeft, (complex8_t *)eigVecs, &rank, work, &lWork,
             rWork, &info, 0, 0);
      if (cmd == OptSize) lWork = (int)work[0].re;
      delete [] rWork;  delete [] work; delete [] eigVecLeft;
#  elif defined (USE_ESSL)
      int iOpt = 1;  // compute both vecs and vals
      if (cmd != OptSize)  {  // not supported by ESSL
         if (cmd == ValuesOnly)  iOpt = 0;
         cgeev (iOpt, (complex<float> *)inMat, rank,
                (complex<float> *)eigVals, (complex<float> *)eigVecs, rank,
                NULL, rank, NULL, 0);
         // --- normalize eigenvectors
         double c;
         for (int i=0; i < n; i++)  {
            c = 1. / norm2 (&eigVecs[i*n], n);
            scale (&eigVecs[i*n], c, n);
         }
      }
#  elif defined (USE_INTEL_MKL)
      MKL_Complex8 *eigVecLeft = new MKL_Complex8 [ldVecLeft];
      MKL_Complex8 *work       = new MKL_Complex8 [workDim];
      float        *rWork      = new float        [workDim];
      cgeev (&jobvl, &jobvr,
             &rank, (MKL_Complex8 *)inMat, &rank,
             (MKL_Complex8 *)eigVals, eigVecLeft,
             &ldVecLeft, (MKL_Complex8 *)eigVecs, &rank, work, &lWork,
             rWork, &info);
      if (cmd == OptSize) lWork = (int)lround(work[0].real);
      delete [] rWork;
      delete [] work;
      delete [] eigVecLeft;
#  elif defined (USE_ACML)
      complex *eigVecLeft = new complex [ldVecLeft];
      cgeev (jobvl, jobvr,
             rank, (complex *)inMat, rank,
             (complex *)eigVals, eigVecLeft,
             ldVecLeft, (complex *)eigVecs, rank, &info);
      delete [] eigVecLeft;
#  elif defined (USE_NETLIB)
      SxComplex8 *eigVecLeft = new SxComplex8[ldVecLeft];
      SxComplex8 *work       = new SxComplex8[workDim];
      float *rWork = new float[workDim];
      info = LAPACKE_cgeev_work (LAPACK_COL_MAJOR,
             jobvl, jobvr,
             rank, inMat, rank, eigVals, eigVecLeft,
             ldVecLeft, eigVecs, rank, work, lWork, rWork);
      if (cmd == OptSize) lWork = (int)lround(work[0].re);
      delete [] rWork;
      delete [] work;
      delete [] eigVecLeft;
#  else
      complex *eigVecLeft = new complex[ldVecLeft];
      complex *work       = new complex[workDim];
      float *rWork = new float[workDim];
      cgeev_ (&jobvl, &jobvr,
              &rank, (complex *)inMat, &rank,
              (complex *)eigVals, eigVecLeft,
              &ldVecLeft, (complex *)eigVecs, &rank, work, &lWork,
              rWork, &info);
      if (cmd == OptSize) lWork = (int)work[0].r;
      delete [] rWork;  delete [] work; delete [] eigVecLeft;
#  endif

#  ifndef USE_ESSL
      if ( info )  {  // TODO: throw exception
         std::cout << "matEigensolver: Error in cgeev: " << info << std::endl;
         SX_EXIT;
      }
#  endif

#  if   defined (USE_VECLIB)
      return (cmd == OptSize)  ?  lWork : 0;
#  elif defined (USE_ESSL)
      return 0;
#  elif defined (USE_INTEL_MKL)
      return (cmd == OptSize)  ?  lWork : 0;
#  elif defined (USE_ACML)
      return 0;
#  else
      return (cmd == OptSize)  ?  (int)lWork : 0;
#  endif
}

void matSymEigensolver (float *eigVals, SxComplex8 *eigVecs,
                        SxComplex8 *inMat, int n, EIGCMD cmd)
{
   char jobz = (cmd == ValuesOnly) ? 'N' : 'V';
   char range = 'A'; // all eigenvalues
   char uplo = 'U';  // use upper right half
   float vdummy = 0., abstol = 0.;
   // abstol = dlamch ('Safe minimum');
   // use automatic workspace LAPACKE interface where available
#  if defined (USE_NETLIB) || defined (USE_ESSL) || defined (USE_ACML) \
   || defined (USE_INTEL_MKL)
#     define SX_USE_LAPACKE
#  endif

   // --- workspace query
#  ifdef SX_USE_LAPACKE
      lapack_int idummy = 0;
      lapack_int nFound = 0;
      lapack_int info = 0;
# else
      int idummy = 0;
      int nFound = 0;
      int info = 0;
      int liwork = -1;
      int queryWork = -1;
      SxComplex8 lwork_ = -1.;
      float lrwork_ = -1.;
#     if defined (USE_VECLIB)
         SX_EXIT;
#     else
         cheevr_ (&jobz, &range, &uplo,
             &n, (complex *)inMat, &n,
             (real *)&vdummy, (real *)&vdummy,
             &idummy, &idummy, (real *)&abstol,
             &nFound, (real *)eigVals, (complex *)eigVecs, &n,
             NULL, (complex *)&lwork_, &queryWork,
             (real*)&lrwork_, &queryWork, &liwork, &queryWork, &info);
#     endif
      if ( info )  {  // TODO: throw exception
         std::cout << "matSymEigensolver: Error in csyevr: " << info << std::endl;
         SX_EXIT;
      }
      // --- allocate workspaces
      int *iwork = new int[liwork];
      int lwork = (int)lround(lwork_.re);
      SxComplex8 *work = new SxComplex8[lwork];
      int lrwork = (int)lround(lrwork_);
      float *rwork = new float[lwork];
#  endif

   // --- computational call
#  if   defined (USE_VECLIB)
      SX_EXIT;
#  elif defined (SX_USE_LAPACKE)
      // allocate iSuppZ workspace
      lapack_int *iSuppZ = new lapack_int[2 * n];
      info = LAPACKE_cheevr (LAPACK_COL_MAJOR,
             jobz, range, uplo,
             n, (lapack_complex_float *)inMat, n,
             vdummy, vdummy, idummy, idummy, abstol,
             &nFound, eigVals, (lapack_complex_float *) eigVecs, n,
             iSuppZ);
#  else
      // allocate iSuppZ workspace
      int *iSuppZ = new int[2 * n];
      cheevr_ (&jobz, &range, &uplo,
             &n, (complex *)inMat, &n,
             (real *)&vdummy, (real *)&vdummy,
             &idummy, &idummy, (real *)&abstol,
             &nFound, (real *)eigVals, (complex *)eigVecs, &n,
             iSuppZ, (complex *)work, &lwork, (real *)rwork, &lrwork,
             iwork, &liwork, &info);
#  endif

   // --- delete workspaces
   delete [] iSuppZ;
#  ifndef SX_USE_LAPACKE
      delete [] rwork;
      delete [] work;
      delete [] iwork;
#  endif

   // --- check for success
   if ( info )  {  // TODO: throw exception
      std::cout << "matSymEigensolver: Error in csyevr: "
                << info << std::endl;
      SX_EXIT;
   }
#  ifdef SX_USE_LAPACKE
#     undef SX_USE_LAPACKE
#  endif
}

int matEigensolver (SxComplex16 *eigVals, SxComplex16 *eigVecs,
                    SxComplex16 *inMat, int n, EIGCMD cmd, int size)
{
#  if defined (USE_VECLIB) || defined (USE_INTEL_MKL) || defined (USE_NETLIB)
      int rank  = n, info = 0, ldVecLeft = 1;
      int lWork = !size ? 2*n : size;
      int workDim = lWork;
#  elif defined (USE_ESSL)
      int rank  = n, ldVecLeft = 1;
      int lWork = !size ? 2*n : size;
      int workDim = lWork;
#  elif defined (USE_ACML)
      int rank  = n, info = 0, ldVecLeft = 1, lWork=size;
#  else
      integer rank  = (integer)n, info = 0, ldVecLeft = 1;
      integer lWork = !size ? 2*n : size;
      integer workDim = lWork;
#  endif
   if ( cmd == OptSize )  lWork = -1;
   char jobvl = 'N', jobvr = 'V';

#  if   defined (USE_VECLIB)
      complex16_t *eigVecLeft = new complex16_t [ldVecLeft];
      complex16_t *work       = new complex16_t [workDim];
      double      *rWork      = new double      [workDim];
      zgeev (&jobvl, &jobvr,
             &rank, (complex16_t *)inMat, &rank,
             (complex16_t *)eigVals, eigVecLeft,
             &ldVecLeft, (complex16_t *)eigVecs, &rank, work, &lWork,
             rWork, &info, 0, 0);
      if (cmd == OptSize) lWork = (int)work[0].re;
      delete [] rWork;  delete [] work; delete [] eigVecLeft;
#  elif defined (USE_ESSL)
      int iOpt = 1;  // compute both vecs and vals
      if (cmd != OptSize)  {  // not supported by ESSL
         if (cmd == ValuesOnly)  iOpt = 0;
         zgeev (iOpt, (complex<double> *)inMat, rank,
                (complex<double> *)eigVals, (complex<double> *)eigVecs, rank,
                NULL, rank, NULL, 0);
         // --- normalize eigenvectors
         double c;
         for (int i=0; i < n; i++)  {
            c = 1. / norm2 (&eigVecs[i*n], n);
            scale (&eigVecs[i*n], c, n);
         }
      }
#  elif defined (USE_INTEL_MKL)
      MKL_Complex16 *eigVecLeft = new MKL_Complex16 [ldVecLeft];
      MKL_Complex16 *work       = new MKL_Complex16 [workDim];
      double        *rWork      = new double        [workDim];
      zgeev (&jobvl, &jobvr,
             &rank, (MKL_Complex16 *)inMat, &rank,
             (MKL_Complex16 *)eigVals, eigVecLeft,
             &ldVecLeft, (MKL_Complex16 *)eigVecs, &rank, work, &lWork,
             rWork, &info);
      if (cmd == OptSize) lWork = (int)lround(work[0].real);
      delete [] rWork;
      delete [] work;
      delete [] eigVecLeft;
#  elif defined (USE_ACML)
      doublecomplex *eigVecLeft = new doublecomplex [ldVecLeft];
      zgeev (jobvl, jobvr,
             rank, (doublecomplex *)inMat, rank,
             (doublecomplex *)eigVals, eigVecLeft,
             ldVecLeft, (doublecomplex *)eigVecs, rank, &info);
      delete [] eigVecLeft;
#  elif defined (USE_NETLIB)
      SxComplex16 *eigVecLeft = new SxComplex16[ldVecLeft];
      SxComplex16 *work       = new SxComplex16[workDim];
      double *rWork = new double[workDim];
      info = LAPACKE_zgeev_work (LAPACK_COL_MAJOR,
             jobvl, jobvr,
             rank, inMat, rank, eigVals, eigVecLeft,
             ldVecLeft, eigVecs, rank, work, lWork, rWork);
      if (cmd == OptSize) lWork = (int)lround(work[0].re);
      delete [] rWork;
      delete [] work;
      delete [] eigVecLeft;
#  else
      doublecomplex *eigVecLeft = new doublecomplex[ldVecLeft];
      doublecomplex *work       = new doublecomplex[workDim];
      double *rWork = new double[workDim];
      zgeev_ (&jobvl, &jobvr,
              &rank, (doublecomplex *)inMat, &rank,
              (doublecomplex *)eigVals, eigVecLeft,
              &ldVecLeft, (doublecomplex *)eigVecs, &rank, work, &lWork,
              rWork, &info);
      if (cmd == OptSize) lWork = (int)work[0].r;
      delete [] rWork;  delete [] work; delete [] eigVecLeft;
#  endif

#  ifndef USE_ESSL
      if ( info )  {  // TODO: throw exception
         std::cout << "matEigensolver: Error in zgeev: " << info << std::endl;
         SX_EXIT;
      }
#  endif

#  if   defined (USE_VECLIB)
      return (cmd == OptSize)  ?  lWork : 0;
#  elif defined (USE_ESSL)
      return 0;
#  elif defined (USE_INTEL_MKL)
      return (cmd == OptSize)  ?  lWork : 0;
#  elif defined (USE_ACML)
      return 0;
#  else
      return (cmd == OptSize)  ?  (int)lWork : 0;
#  endif
}

void matSymEigensolver (double *eigVals, SxComplex16 *eigVecs,
                        SxComplex16 *inMat, int n, EIGCMD cmd)
{
   char jobz = (cmd == ValuesOnly) ? 'N' : 'V';
   char range = 'A'; // all eigenvalues
   char uplo = 'U';  // use upper right half
   double vdummy = 0., abstol = 0.;
   // abstol = dlamch ('Safe minimum');
   // use automatic workspace LAPACKE interface where available
#  if defined (USE_NETLIB) || defined (USE_ESSL) || defined (USE_ACML) \
   || defined (USE_INTEL_MKL)
#     define SX_USE_LAPACKE
#  endif

   // --- workspace query
#  ifdef SX_USE_LAPACKE
      lapack_int idummy = 0;
      lapack_int nFound = 0;
      lapack_int info = 0;
# else
      int idummy = 0;
      int nFound = 0;
      int info = 0;
      int liwork = -1;
      int queryWork = -1;
      SxComplex16 lwork_ = -1.;
      double lrwork_ = -1.;
#     if defined (USE_VECLIB)
         SX_EXIT;
#     else
         zheevr_ (&jobz, &range, &uplo,
             &n, (doublecomplex *)inMat, &n,
             (doublereal *)&vdummy, (doublereal *)&vdummy,
             &idummy, &idummy, (doublereal *)&abstol,
             &nFound, (doublereal *)eigVals, (doublecomplex *)eigVecs, &n,
             NULL, (doublecomplex *)&lwork_, &queryWork,
             (doublereal*)&lrwork_, &queryWork, &liwork, &queryWork, &info);
#     endif
      if ( info )  {  // TODO: throw exception
         std::cout << "matSymEigensolver: Error in zsyevr: " << info << std::endl;
         SX_EXIT;
      }
      // --- allocate workspaces
      int *iwork = new int[liwork];
      int lwork = (int)lround(lwork_.re);
      SxComplex16 *work = new SxComplex16[lwork];
      int lrwork = (int)lround(lrwork_);
      double *rwork = new double[lwork];
#  endif

   // --- computational call
#  if   defined (USE_VECLIB)
      SX_EXIT;
#  elif defined (SX_USE_LAPACKE)
      // allocate iSuppZ workspace
      lapack_int *iSuppZ = new lapack_int[2 * n];
      info = LAPACKE_zheevr (LAPACK_COL_MAJOR,
             jobz, range, uplo,
             n, (lapack_complex_double *)inMat, n,
             vdummy, vdummy, idummy, idummy, abstol,
             &nFound, eigVals, (lapack_complex_double *) eigVecs, n,
             iSuppZ);
#  else
      // allocate iSuppZ workspace
      int *iSuppZ = new int[2 * n];
      zheevr_ (&jobz, &range, &uplo,
             &n, (doublecomplex *)inMat, &n,
             (doublereal *)&vdummy, (doublereal *)&vdummy,
             &idummy, &idummy, (doublereal *)&abstol,
             &nFound, (doublereal *)eigVals, (doublecomplex *)eigVecs, &n,
             iSuppZ, (doublecomplex *)work, &lwork, (doublereal *)rwork, &lrwork,
             iwork, &liwork, &info);
#  endif

   // --- delete workspaces
   delete [] iSuppZ;
#  ifndef SX_USE_LAPACKE
      delete [] rwork;
      delete [] work;
      delete [] iwork;
#  endif

   // --- check for success
   if ( info )  {  // TODO: throw exception
      std::cout << "matSymEigensolver: Error in zsyevr: "
                << info << std::endl;
      SX_EXIT;
   }
#  ifdef SX_USE_LAPACKE
#     undef SX_USE_LAPACKE
#  endif
}
// --- EIGENCOMPLEX

//------------------------------------------------------------------------------
// Eigensolver - trigonal matrices
//------------------------------------------------------------------------------
// --- The next 150 lines were generated from math/snippets/SxBlasLib.cpp snippet EIGEN3_REAL
void matEigensolverTri (float *eigVals, float *eigVecs,
                        float *inMat, int n, enum UPLO uplo,
                        EIGCMD
#                       ifdef USE_ESSL
                           cmd
#                       endif
                        )
{
   char jobz = 'V';
   char uploChar = (uplo == UpperRight ? 'U' : 'L');
#  if   defined (USE_VECLIB)
      int rank = n;
      int ldEigVecs = rank;
      int info = 0;
      int workDim  = 3*n;
      float *work  = new float [workDim];
      sspev (&jobz, &uploChar,
             &rank, (float *)inMat, (float *)eigVals,
             (float *)eigVecs, &ldEigVecs,
              work, &info, 0, 0);
      delete [] work;
#  elif defined (USE_ESSL)
      int info = 0, iOpt = 1;
      if (cmd  == VectorsOnly)  iOpt = 0;
      if (uplo == UpperRight)   iOpt += 20;
      int rank = n;
      sspev (iOpt, inMat, eigVals,
             eigVecs, rank, n, NULL, 0);
#  elif defined (USE_INTEL_MKL)
      int rank = n;
      int ldEigVecs = rank;
      int info = 0;
      int workDim  = 3*n;
      float *work  = new float [workDim];
      sspev (&jobz, &uploChar,
             &rank, (float *)inMat, (float *)eigVals,
             (float *)eigVecs, &ldEigVecs,
              work, &info);
      delete [] work;
#  elif defined (USE_ACML)
      int rank = n;
      int ldEigVecs = rank;
      int info = 0;
      sspev (jobz, uploChar,
             rank, inMat, eigVals,
             eigVecs, ldEigVecs, &info);
#  elif defined (USE_NETLIB)
      int rank = n;
      int ldEigVecs = rank;
      int info = 0;
      int workDim  = 3*n;
      float *work  = new float [workDim];
      info = LAPACKE_sspev_work (LAPACK_COL_MAJOR,
               jobz, uploChar, rank, inMat,
               eigVals, eigVecs, ldEigVecs, work);
      delete [] work;
#  else
      integer rank = (integer)n;
      integer ldEigVecs = rank;
      integer info = 0;
      integer workDim  = 3*n;
      real *work  = new real [workDim];
      sspev_ (&jobz, &uploChar,
              &rank, (real *)inMat, (real *)eigVals,
              (real *)eigVecs, &ldEigVecs,
               work, &info);
      delete [] work;
#  endif

   if ( info )  {  // TODO: throw exception
      std::cout << "matEigensolverHerm: Error in sspev: " << info << std::endl;
      SX_EXIT;
   }
}

void matEigensolverTri (double *eigVals, double *eigVecs,
                        double *inMat, int n, enum UPLO uplo,
                        EIGCMD
#                       ifdef USE_ESSL
                           cmd
#                       endif
                        )
{
   char jobz = 'V';
   char uploChar = (uplo == UpperRight ? 'U' : 'L');
#  if   defined (USE_VECLIB)
      int rank = n;
      int ldEigVecs = rank;
      int info = 0;
      int workDim  = 3*n;
      double *work  = new double [workDim];
      dspev (&jobz, &uploChar,
             &rank, (double *)inMat, (double *)eigVals,
             (double *)eigVecs, &ldEigVecs,
              work, &info, 0, 0);
      delete [] work;
#  elif defined (USE_ESSL)
      int info = 0, iOpt = 1;
      if (cmd  == VectorsOnly)  iOpt = 0;
      if (uplo == UpperRight)   iOpt += 20;
      int rank = n;
      dspev (iOpt, inMat, eigVals,
             eigVecs, rank, n, NULL, 0);
#  elif defined (USE_INTEL_MKL)
      int rank = n;
      int ldEigVecs = rank;
      int info = 0;
      int workDim  = 3*n;
      double *work  = new double [workDim];
      dspev (&jobz, &uploChar,
             &rank, (double *)inMat, (double *)eigVals,
             (double *)eigVecs, &ldEigVecs,
              work, &info);
      delete [] work;
#  elif defined (USE_ACML)
      int rank = n;
      int ldEigVecs = rank;
      int info = 0;
      dspev (jobz, uploChar,
             rank, inMat, eigVals,
             eigVecs, ldEigVecs, &info);
#  elif defined (USE_NETLIB)
      int rank = n;
      int ldEigVecs = rank;
      int info = 0;
      int workDim  = 3*n;
      double *work  = new double [workDim];
      info = LAPACKE_dspev_work (LAPACK_COL_MAJOR,
               jobz, uploChar, rank, inMat,
               eigVals, eigVecs, ldEigVecs, work);
      delete [] work;
#  else
      integer rank = (integer)n;
      integer ldEigVecs = rank;
      integer info = 0;
      integer workDim  = 3*n;
      doublereal *work  = new doublereal [workDim];
      dspev_ (&jobz, &uploChar,
              &rank, (doublereal *)inMat, (doublereal *)eigVals,
              (doublereal *)eigVecs, &ldEigVecs,
               work, &info);
      delete [] work;
#  endif

   if ( info )  {  // TODO: throw exception
      std::cout << "matEigensolverHerm: Error in dspev: " << info << std::endl;
      SX_EXIT;
   }
}
// --- EIGEN3_REAL

// --- The next 166 lines were generated from snippets/SxBlasLib.cpp snippet EIGEN3_COMPLEX
void matEigensolverTri (float *eigVals, SxComplex8 *eigVecs,
                        SxComplex8 *inMat, int n, enum UPLO uplo,
                        EIGCMD
#                       ifdef USE_ESSL
                           cmd
#                       endif
                        )
{
   char jobz = 'V';
   char uploChar = (uplo == UpperRight ? 'U' : 'L');
#  if   defined (USE_VECLIB)
      int rank = n;
      int ldEigVecs = rank;
      int info = 0;
      int workDim  = 2*n-1;
      int rWorkDim = 3*n-2;
      complex8_t *work = new complex8_t [workDim];
      float    *rWork  = new float[rWorkDim];
      chpev (&jobz, &uploChar,
             &rank, (complex8_t *)inMat, (float *)eigVals,
             (complex8_t *)eigVecs, &ldEigVecs,
              work, rWork, &info, 0, 0);
      delete [] rWork; delete [] work;
#  elif defined (USE_ESSL)
      int info = 0, iOpt = 1;
      if (cmd  == VectorsOnly)  iOpt = 0;
      if (uplo == UpperRight)   iOpt += 20;
      int rank = n;
      chpev (iOpt, (complex<float> *)inMat, eigVals,
             (complex<float> *)eigVecs, rank, n, NULL, 0);
#  elif defined (USE_INTEL_MKL)
      int rank = n;
      int ldEigVecs = rank;
      int info = 0;
      int workDim  = 2*n-1;
      int rWorkDim = 3*n-2;
      MKL_Complex8 *work = new MKL_Complex8 [workDim];
      float       *rWork = new float[rWorkDim];
      chpev (&jobz, &uploChar,
             &rank, (MKL_Complex8 *)inMat, (float *)eigVals,
             (MKL_Complex8 *)eigVecs, &ldEigVecs,
              work, rWork, &info);
      delete [] rWork; delete [] work;
#  elif defined (USE_ACML)
      int rank = n;
      int ldEigVecs = rank;
      int info = 0;
      chpev (jobz, uploChar,
             rank, (complex *)inMat, (float *)eigVals,
             (complex *)eigVecs, ldEigVecs, &info);
#  elif defined (USE_NETLIB)
      int rank = n;
      int ldEigVecs = rank;
      int info = 0;
      int workDim  = 2*n-1;
      int rWorkDim = 3*n-2;
      SxComplex8 *work = new SxComplex8[workDim];
      float *rWork = new float[rWorkDim];
      info = LAPACKE_chpev_work (LAPACK_COL_MAJOR,
               jobz, uploChar, rank, inMat,
               eigVals, eigVecs, ldEigVecs, work, rWork);
      delete [] rWork; delete [] work;
#  else
      integer rank = (integer)n;
      integer ldEigVecs = rank;
      integer info = 0;
      integer workDim  = 2*n-1;
      integer rWorkDim = 3*n-2;
      complex *work       = new complex [workDim];
      real *rWork      = new real [rWorkDim];
      chpev_ (&jobz, &uploChar,
              &rank, (complex *)inMat, (real *)eigVals,
              (complex *)eigVecs, &ldEigVecs,
               work, rWork, &info);
      delete [] rWork; delete [] work;
#  endif

   if ( info )  {  // TODO: throw exception
      std::cout << "matEigensolverHerm: Error in chpev: " << info << std::endl;
      SX_EXIT;
   }
}

void matEigensolverTri (double *eigVals, SxComplex16 *eigVecs,
                        SxComplex16 *inMat, int n, enum UPLO uplo,
                        EIGCMD
#                       ifdef USE_ESSL
                           cmd
#                       endif
                        )
{
   char jobz = 'V';
   char uploChar = (uplo == UpperRight ? 'U' : 'L');
#  if   defined (USE_VECLIB)
      int rank = n;
      int ldEigVecs = rank;
      int info = 0;
      int workDim  = 2*n-1;
      int rWorkDim = 3*n-2;
      complex16_t *work = new complex16_t [workDim];
      double    *rWork  = new double[rWorkDim];
      zhpev (&jobz, &uploChar,
             &rank, (complex16_t *)inMat, (double *)eigVals,
             (complex16_t *)eigVecs, &ldEigVecs,
              work, rWork, &info, 0, 0);
      delete [] rWork; delete [] work;
#  elif defined (USE_ESSL)
      int info = 0, iOpt = 1;
      if (cmd  == VectorsOnly)  iOpt = 0;
      if (uplo == UpperRight)   iOpt += 20;
      int rank = n;
      zhpev (iOpt, (complex<double> *)inMat, eigVals,
             (complex<double> *)eigVecs, rank, n, NULL, 0);
#  elif defined (USE_INTEL_MKL)
      int rank = n;
      int ldEigVecs = rank;
      int info = 0;
      int workDim  = 2*n-1;
      int rWorkDim = 3*n-2;
      MKL_Complex16 *work = new MKL_Complex16 [workDim];
      double       *rWork = new double[rWorkDim];
      zhpev (&jobz, &uploChar,
             &rank, (MKL_Complex16 *)inMat, (double *)eigVals,
             (MKL_Complex16 *)eigVecs, &ldEigVecs,
              work, rWork, &info);
      delete [] rWork; delete [] work;
#  elif defined (USE_ACML)
      int rank = n;
      int ldEigVecs = rank;
      int info = 0;
      zhpev (jobz, uploChar,
             rank, (doublecomplex *)inMat, (double *)eigVals,
             (doublecomplex *)eigVecs, ldEigVecs, &info);
#  elif defined (USE_NETLIB)
      int rank = n;
      int ldEigVecs = rank;
      int info = 0;
      int workDim  = 2*n-1;
      int rWorkDim = 3*n-2;
      SxComplex16 *work = new SxComplex16[workDim];
      double *rWork = new double[rWorkDim];
      info = LAPACKE_zhpev_work (LAPACK_COL_MAJOR,
               jobz, uploChar, rank, inMat,
               eigVals, eigVecs, ldEigVecs, work, rWork);
      delete [] rWork; delete [] work;
#  else
      integer rank = (integer)n;
      integer ldEigVecs = rank;
      integer info = 0;
      integer workDim  = 2*n-1;
      integer rWorkDim = 3*n-2;
      doublecomplex *work       = new doublecomplex [workDim];
      doublereal *rWork      = new doublereal [rWorkDim];
      zhpev_ (&jobz, &uploChar,
              &rank, (doublecomplex *)inMat, (doublereal *)eigVals,
              (doublecomplex *)eigVecs, &ldEigVecs,
               work, rWork, &info);
      delete [] rWork; delete [] work;
#  endif

   if ( info )  {  // TODO: throw exception
      std::cout << "matEigensolverHerm: Error in zhpev: " << info << std::endl;
      SX_EXIT;
   }
}
// --- EIGEN3_COMPLEX

// --- multiply packed symmetric matrix with a vector
//     target vector is always stride 1
// --- The next 144 lines were generated from snippets/SxBlasLib.cpp snippet symmatmult
void symmatmult (int n, float *resVec, const float *aPacked, enum UPLO uplo,
                 const float *xVec, int xStride)
{
   float alpha = 1.0;
   float beta = 0.0;
#  if defined USE_VECLIB || defined USE_ESSL || defined USE_INTEL_MKL || defined USE_ACML
   char uploChar = (uplo == UpperRight ? 'U' : 'L');
#  endif
   int resStride = 1;
#  if   defined (USE_VECLIB)
      sspmv (&uploChar, &n, &alpha,
             const_cast<float*>(aPacked),
             const_cast<float*>(xVec), &xStride,
             &beta, resVec, &resStride);
#  elif defined (USE_ESSL)
      sspmv (&uploChar, &n, &alpha,
             const_cast<float*>(aPacked),
             const_cast<float*>(xVec), &xStride,
             &beta, resVec, &resStride);
#  elif defined (USE_INTEL_MKL)
      sspmv (&uploChar, &n, &alpha,
             const_cast<float*>(aPacked),
             const_cast<float*>(xVec), &xStride,
             &beta, resVec, &resStride);
#  elif defined (USE_ACML)
      sspmv (uploChar, n, alpha,
             const_cast<float*>(aPacked),
             const_cast<float*>(xVec), xStride,
             beta, resVec, resStride);
#  else /* ATLAS / GOTO / NETLIB */
      cblas_sspmv (CblasColMajor, (uplo == UpperRight ? CblasUpper : CblasLower),
                   n, alpha, aPacked, xVec, xStride,
                   beta, resVec, resStride);
#  endif
}

void symmatmult (int n, double *resVec, const double *aPacked, enum UPLO uplo,
                 const double *xVec, int xStride)
{
   double alpha = 1.0;
   double beta = 0.0;
#  if defined USE_VECLIB || defined USE_ESSL || defined USE_INTEL_MKL || defined USE_ACML
   char uploChar = (uplo == UpperRight ? 'U' : 'L');
#  endif
   int resStride = 1;
#  if   defined (USE_VECLIB)
      dspmv (&uploChar, &n, &alpha,
             const_cast<double*>(aPacked),
             const_cast<double*>(xVec), &xStride,
             &beta, resVec, &resStride);
#  elif defined (USE_ESSL)
      dspmv (&uploChar, &n, &alpha,
             const_cast<double*>(aPacked),
             const_cast<double*>(xVec), &xStride,
             &beta, resVec, &resStride);
#  elif defined (USE_INTEL_MKL)
      dspmv (&uploChar, &n, &alpha,
             const_cast<double*>(aPacked),
             const_cast<double*>(xVec), &xStride,
             &beta, resVec, &resStride);
#  elif defined (USE_ACML)
      dspmv (uploChar, n, alpha,
             const_cast<double*>(aPacked),
             const_cast<double*>(xVec), xStride,
             beta, resVec, resStride);
#  else /* ATLAS / GOTO / NETLIB */
      cblas_dspmv (CblasColMajor, (uplo == UpperRight ? CblasUpper : CblasLower),
                   n, alpha, aPacked, xVec, xStride,
                   beta, resVec, resStride);
#  endif
}

void symmatmult (int n, SxComplex8 *resVec, const SxComplex8 *aPacked, enum UPLO uplo,
                 const SxComplex8 *xVec, int xStride)
{
   SxComplex8 alpha = 1.0;
   SxComplex8 beta = 0.0;
#  if defined USE_VECLIB || defined USE_ESSL || defined USE_INTEL_MKL || defined USE_ACML
   char uploChar = (uplo == UpperRight ? 'U' : 'L');
#  endif
   int resStride = 1;
#  if   defined (USE_VECLIB)
      chpmv (&uploChar, &n, &alpha,
             (complex8_t*)const_cast<SxComplex8*>(aPacked),
             (complex8_t*)const_cast<SxComplex8*>(xVec), &xStride,
             &beta, (complex8_t*)resVec, &resStride);
#  elif defined (USE_ESSL)
      chpmv (&uploChar, &n, &alpha,
             (complex<float>*)const_cast<SxComplex8*>(aPacked),
             (complex<float>*)const_cast<SxComplex8*>(xVec), &xStride,
             &beta, (complex<float>*)resVec, &resStride);
#  elif defined (USE_INTEL_MKL)
      chpmv (&uploChar, &n, (MKL_Complex8*)&alpha,
             (MKL_Complex8*)const_cast<SxComplex8*>(aPacked),
             (MKL_Complex8*)const_cast<SxComplex8*>(xVec), &xStride,
             (MKL_Complex8*)&beta, (MKL_Complex8*)resVec, &resStride);
#  elif defined (USE_ACML)
      chpmv (uploChar, n, (complex*)&alpha,
             (complex*)const_cast<SxComplex8*>(aPacked),
             (complex*)const_cast<SxComplex8*>(xVec), xStride,
             (complex*)&beta, (complex*)resVec, resStride);
#  else /* ATLAS / GOTO / NETLIB */
      cblas_chpmv (CblasColMajor, (uplo == UpperRight ? CblasUpper : CblasLower),
                   n, &alpha, aPacked, xVec, xStride,
                   &beta, resVec, resStride);
#  endif
}

void symmatmult (int n, SxComplex16 *resVec, const SxComplex16 *aPacked, enum UPLO uplo,
                 const SxComplex16 *xVec, int xStride)
{
   SxComplex16 alpha = 1.0;
   SxComplex16 beta = 0.0;
#  if defined USE_VECLIB || defined USE_ESSL || defined USE_INTEL_MKL || defined USE_ACML
   char uploChar = (uplo == UpperRight ? 'U' : 'L');
#  endif
   int resStride = 1;
#  if   defined (USE_VECLIB)
      zhpmv (&uploChar, &n, &alpha,
             (complex16_t*)const_cast<SxComplex16*>(aPacked),
             (complex16_t*)const_cast<SxComplex16*>(xVec), &xStride,
             &beta, (complex16_t*)resVec, &resStride);
#  elif defined (USE_ESSL)
      zhpmv (&uploChar, &n, &alpha,
             (complex<double>*)const_cast<SxComplex16*>(aPacked),
             (complex<double>*)const_cast<SxComplex16*>(xVec), &xStride,
             &beta, (complex<double>*)resVec, &resStride);
#  elif defined (USE_INTEL_MKL)
      zhpmv (&uploChar, &n, (MKL_Complex16*)&alpha,
             (MKL_Complex16*)const_cast<SxComplex16*>(aPacked),
             (MKL_Complex16*)const_cast<SxComplex16*>(xVec), &xStride,
             (MKL_Complex16*)&beta, (MKL_Complex16*)resVec, &resStride);
#  elif defined (USE_ACML)
      zhpmv (uploChar, n, (doublecomplex*)&alpha,
             (doublecomplex*)const_cast<SxComplex16*>(aPacked),
             (doublecomplex*)const_cast<SxComplex16*>(xVec), xStride,
             (doublecomplex*)&beta, (doublecomplex*)resVec, resStride);
#  else /* ATLAS / GOTO / NETLIB */
      cblas_zhpmv (CblasColMajor, (uplo == UpperRight ? CblasUpper : CblasLower),
                   n, &alpha, aPacked, xVec, xStride,
                   &beta, resVec, resStride);
#  endif
}
// --- symmatmult

// -------------------------- SX_IGNORE_THE_REST_OF_THE_FILE
# endif
// --------------------------
