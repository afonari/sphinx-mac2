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

#define SXVEC_INSTANCE template
#include <SxVector.h>

// --- The next 300 lines were generated from snippets/SxVecCompute.cpp snippet Simd_abs_real
// Absolute value
SxVector<int> SxVecCompute<Simd>::abs (const SxVecRef<int,Compact> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(int) * n);

   int* resPtr = (int*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         int x = vec.elements[i];
         resPtr[i] = ::abs(x) ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         int x = vec.elements[i];
         resPtr[i] = ::abs(x) ;
      }
   return SxVector<int>(mem, vec, vec.auxData, vec.auxData);
}

// Absolute value
SxVector<int>
SxVecCompute<Simd>::absInPlace (SxVecRef<int,Compact> &&vec)
{
   SX_CHECK (vec.allocMem);
   if (vec.allocMem->refCounter > 1) return abs(vec);
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         int x = vec.elements[i];
         vec.elements[i] = ::abs(x) ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         int x = vec.elements[i];
         vec.elements[i] = ::abs(x) ;
      }
   return std::move(vec);
}

// Absolute value squared
SxVector<int> SxVecCompute<Simd>::absSqr (const SxVecRef<int,Compact> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(int) * n);

   int* resPtr = (int*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         int x = vec.elements[i];
         resPtr[i] = x*x ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         int x = vec.elements[i];
         resPtr[i] = x*x ;
      }
   return SxVector<int>(mem, vec, vec.auxData, vec.auxData);
}

// Absolute value squared
SxVector<int>
SxVecCompute<Simd>::absSqrInPlace (SxVecRef<int,Compact> &&vec)
{
   SX_CHECK (vec.allocMem);
   if (vec.allocMem->refCounter > 1) return absSqr(vec);
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         int x = vec.elements[i];
         vec.elements[i] = x*x ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         int x = vec.elements[i];
         vec.elements[i] = x*x ;
      }
   return std::move(vec);
}

// Absolute value
SxVector<float> SxVecCompute<Simd>::abs (const SxVecRef<float,Compact> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(float) * n);

   float* resPtr = (float*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         float x = vec.elements[i];
         resPtr[i] = ::fabsf(x) ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         float x = vec.elements[i];
         resPtr[i] = ::fabsf(x) ;
      }
   return SxVector<float>(mem, vec, vec.auxData, vec.auxData);
}

// Absolute value
SxVector<float>
SxVecCompute<Simd>::absInPlace (SxVecRef<float,Compact> &&vec)
{
   SX_CHECK (vec.allocMem);
   if (vec.allocMem->refCounter > 1) return abs(vec);
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         float x = vec.elements[i];
         vec.elements[i] = ::fabsf(x) ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         float x = vec.elements[i];
         vec.elements[i] = ::fabsf(x) ;
      }
   return std::move(vec);
}

// Absolute value squared
SxVector<float> SxVecCompute<Simd>::absSqr (const SxVecRef<float,Compact> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(float) * n);

   float* resPtr = (float*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         float x = vec.elements[i];
         resPtr[i] = x*x ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         float x = vec.elements[i];
         resPtr[i] = x*x ;
      }
   return SxVector<float>(mem, vec, vec.auxData, vec.auxData);
}

// Absolute value squared
SxVector<float>
SxVecCompute<Simd>::absSqrInPlace (SxVecRef<float,Compact> &&vec)
{
   SX_CHECK (vec.allocMem);
   if (vec.allocMem->refCounter > 1) return absSqr(vec);
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         float x = vec.elements[i];
         vec.elements[i] = x*x ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         float x = vec.elements[i];
         vec.elements[i] = x*x ;
      }
   return std::move(vec);
}

// Absolute value
SxVector<double> SxVecCompute<Simd>::abs (const SxVecRef<double,Compact> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(double) * n);

   double* resPtr = (double*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         double x = vec.elements[i];
         resPtr[i] = ::fabs(x) ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         double x = vec.elements[i];
         resPtr[i] = ::fabs(x) ;
      }
   return SxVector<double>(mem, vec, vec.auxData, vec.auxData);
}

// Absolute value
SxVector<double>
SxVecCompute<Simd>::absInPlace (SxVecRef<double,Compact> &&vec)
{
   SX_CHECK (vec.allocMem);
   if (vec.allocMem->refCounter > 1) return abs(vec);
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         double x = vec.elements[i];
         vec.elements[i] = ::fabs(x) ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         double x = vec.elements[i];
         vec.elements[i] = ::fabs(x) ;
      }
   return std::move(vec);
}

// Absolute value squared
SxVector<double> SxVecCompute<Simd>::absSqr (const SxVecRef<double,Compact> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(double) * n);

   double* resPtr = (double*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         double x = vec.elements[i];
         resPtr[i] = x*x ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         double x = vec.elements[i];
         resPtr[i] = x*x ;
      }
   return SxVector<double>(mem, vec, vec.auxData, vec.auxData);
}

// Absolute value squared
SxVector<double>
SxVecCompute<Simd>::absSqrInPlace (SxVecRef<double,Compact> &&vec)
{
   SX_CHECK (vec.allocMem);
   if (vec.allocMem->refCounter > 1) return absSqr(vec);
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         double x = vec.elements[i];
         vec.elements[i] = x*x ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         double x = vec.elements[i];
         vec.elements[i] = x*x ;
      }
   return std::move(vec);
}
// --- Simd_abs_real

#ifdef __GNUG__
#pragma GCC optimize "3"
#include <immintrin.h>
#endif

#if defined __GNUG__ && defined __SSE2__
static inline __m128d sum_simd (const double *elements, ssize_t n)
{
#ifdef __AVX__
   SX_CHECK (n % 16 == 0, n);
   __m256d r1 = {0., 0., 0., 0.},
           r2 = {0., 0., 0., 0.},
           r3 = {0., 0., 0., 0.},
           r4 = {0., 0., 0., 0.};
   for (ssize_t i = 0; i < n; i+= 16)  {
      r1 = _mm256_add_pd(r1, _mm256_loadu_pd(elements + i));
      r2 = _mm256_add_pd(r2, _mm256_loadu_pd(elements + i + 4));
      r3 = _mm256_add_pd(r3, _mm256_loadu_pd(elements + i + 8));
      r4 = _mm256_add_pd(r4, _mm256_loadu_pd(elements + i + 12));
   }
   r1 = _mm256_add_pd(r1,r2);
   r3 = _mm256_add_pd(r3,r4);
   r1 = _mm256_add_pd(r1,r3);
   return _mm_add_pd(_mm256_castpd256_pd128(r1),_mm256_extractf128_pd(r1,1));
#else // only SSE2
   SX_CHECK (n % 8 == 0, n);
   __m128d r1 = {0., 0.}, r2 = {0., 0.},
           r3 = {0., 0.}, r4 = {0., 0.};
   if (((size_t)elements & 15) == 0)  {
      /// SSE2 aligned
      const __m128d *ptr = (const __m128d*)elements;
      n = n/2;
      for (ssize_t i = 0; i < n; i+= 4)  {
         r1 = _mm_add_pd(r1, ptr[i]);
         r2 = _mm_add_pd(r2, ptr[i + 1]);
         r3 = _mm_add_pd(r3, ptr[i + 2]);
         r4 = _mm_add_pd(r4, ptr[i + 3]);
      }
   } else {
      /// SSE2 unaligned (could be slower)
      for (ssize_t i = 0; i < n; i+= 8)  {
         r1 = _mm_add_pd(r1, _mm_loadu_pd(elements + i));
         r2 = _mm_add_pd(r2, _mm_loadu_pd(elements + i + 2));
         r3 = _mm_add_pd(r3, _mm_loadu_pd(elements + i + 4));
         r4 = _mm_add_pd(r4, _mm_loadu_pd(elements + i + 6));
      }
   }
   r1 = _mm_add_pd(r1,r2);
   r3 = _mm_add_pd(r3,r4);
   r1 = _mm_add_pd(r1,r3);
   return r1;
#endif
}
#   ifdef USE_OPENMP
#      pragma omp declare reduction \
         (+ : __m128d : omp_out = _mm_add_pd(omp_in, omp_out)) \
         initializer(omp_priv = {0., 0. })
#   endif
#endif

double SxVecCompute<Simd>::sum (const double* elements, ssize_t n)
{
#if defined __GNUG__ && defined __AVX__
   double res = 0.;
   // --- clean-up: make n compatible with loop stride (16)
   while (n & 15)  {
      res += *elements++;
      n--;
   }
#     ifdef USE_OPENMP
   if (n > sxChunkSize * 16)  {
      double res_ = 0.;
#     pragma omp parallel reduction(+:res_)
      {
         ssize_t n16 = n / 16;
         ssize_t nChunk = sx_omp_chunk (n16);
         ssize_t iStart = sx_omp_start (nChunk, n16);
         __m128d rx = sum_simd(elements + iStart * 16, nChunk * 16);
         res_ = _mm_cvtsd_f64(_mm_hadd_pd(rx, rx));
      }
      return res + res_;
   }
#     endif
   // no openmp parallelization
   __m128d rx = sum_simd(elements, n);
   return res + _mm_cvtsd_f64(_mm_hadd_pd(rx, rx));

#elif defined __GNUG__ && defined __SSE2__
   if (n<=0) return 0.;
   double res = 0.;
   // SSE2 alignment
   if (((size_t)elements & 15) == 8)  {
      res = *elements++;
      n--;
   }
   // now odd number of elements? => add last one later
   bool addLast = false;
   if (n & 1)  {
      addLast = true;
      n--;
   }
   // => n should be even and elements aligned to 16 bytes
   // --- clean-up: make n compatible with loop stride (8)
   while (n & 7)  {
      res += *elements++;
      n--;
   }
#     ifdef USE_OPENMP
   if (n > sxChunkSize * 8)  {
      double res_ = 0.;
#     pragma omp parallel reduction(+:res_)
      {
         ssize_t n8 = n / 8;
         ssize_t nChunk = sx_omp_chunk (n8);
         ssize_t iStart = sx_omp_start (nChunk, n8);
         __m128d rx = sum_simd(elements + iStart * 8, nChunk * 8);
         res_ = _mm_cvtsd_f64(_mm_add_pd(rx, _mm_unpackhi_pd(rx,rx)));
      }
      res += res_;
   } else
#     endif
   {
      __m128d rx = sum_simd(elements, n);
      res += _mm_cvtsd_f64(_mm_add_pd(rx, _mm_unpackhi_pd(rx,rx)));
   }
   if (addLast) res += elements[n];
   return res;
#else
   double res = 0.;
#ifdef USE_OPENMP
# pragma omp parallel for reduction(+:res) if (n > sxChunkSize)
#endif
   for (ssize_t i = 0; i < n; i++)
      res += elements[i];
   return res;
#endif
}

SxComplex16 SxVecCompute<Simd>::sum (const SxComplex16* elements, ssize_t n)
{
#if defined __GNUG__ && defined __AVX__
   __m128d res = { 0., 0. };
   // --- clean-up: make n compatible with loop stride (8)
   while (n & 7)  {
      res = _mm_add_pd(res, _mm_loadu_pd((const double*)(elements++)));
      n--;
   }
   __m128d res_ = { 0., 0. };
#     ifdef USE_OPENMP
   if (n > sxChunkSize * 8)  {
#     pragma omp parallel reduction(+:res_)
      {
         ssize_t n8 = n / 8;
         ssize_t nChunk = sx_omp_chunk (n8);
         ssize_t iStart = sx_omp_start (nChunk, n8);
         res_ = sum_simd((const double*)(elements + iStart * 8), nChunk * 16);
      }
   } else
#     endif
   {
      res_ = sum_simd((const double*)elements, n * 2);
   }
   SxComplex16 result;
   _mm_storeu_pd((double*)&result, _mm_add_pd(res,res_));
   return result;
#elif defined __GNUG__ && defined __SSE2__
   __m128d res = {0., 0.};
   // --- clean-up: make n compatible with loop stride (4)
   while (n & 3)  {
      res = _mm_add_pd(res, _mm_loadu_pd((const double*)(elements++)));
      n--;
   }
   __m128d res_ = {0., 0.};
#     ifdef USE_OPENMP
   if (n > sxChunkSize * 4)  {
#     pragma omp parallel reduction(+:res_)
      {
         ssize_t n4 = n / 4;
         ssize_t nChunk = sx_omp_chunk (n4);
         ssize_t iStart = sx_omp_start (nChunk, n4);
         res_ = sum_simd((const double*)(elements + iStart * 4), nChunk * 8);
      }
   } else
#     endif
   {
      res_ = sum_simd((const double*)elements, n * 2);
   }
   SxComplex16 result;
   _mm_storeu_pd((double*)(&result), _mm_add_pd(res,res_));
   return result;
#else
   SxComplex16 res = 0.;
#ifdef USE_OPENMP
# pragma omp parallel for reduction(+:res) if (n > sxChunkSize)
#endif
   for (ssize_t i = 0; i < n; i++)
      res += elements[i];
   return res;
#endif
}

#if defined __GNUG__ && defined __SSE2__
/// Sum squares of elements
static inline __m128d normSqr_simd (const double *elements, ssize_t n)
{
#ifdef __AVX__
   SX_CHECK (n % 16 == 0, n);
   __m256d r1 = {0., 0., 0., 0.},
           r2 = {0., 0., 0., 0.},
           r3 = {0., 0., 0., 0.},
           r4 = {0., 0., 0., 0.};
   for (ssize_t i = 0; i < n; i+= 16)  {
      __m256d x1 = _mm256_loadu_pd(elements + i);
      __m256d x2 = _mm256_loadu_pd(elements + i + 4);
#     ifdef __AVX2__
      r1 = _mm256_fmadd_pd(x1, x1, r1);
      r2 = _mm256_fmadd_pd(x2, x2, r2);
#     else
      r1 = _mm256_add_pd(r1, _mm256_mul_pd(x1, x1));
      r2 = _mm256_add_pd(r2, _mm256_mul_pd(x2, x2));
#     endif
      __m256d x3 = _mm256_loadu_pd(elements + i + 8);
      __m256d x4 = _mm256_loadu_pd(elements + i + 12);
#     ifdef __AVX2__
      r3 = _mm256_fmadd_pd(x3, x3, r3);
      r4 = _mm256_fmadd_pd(x4, x4, r4);
#     else
      r3 = _mm256_add_pd(r3, _mm256_mul_pd(x3, x3));
      r4 = _mm256_add_pd(r4, _mm256_mul_pd(x4, x4));
#     endif
   }
   r1 = _mm256_add_pd(r1,r2);
   r3 = _mm256_add_pd(r3,r4);
   r1 = _mm256_add_pd(r1,r3);
   return _mm_add_pd(_mm256_castpd256_pd128(r1),_mm256_extractf128_pd(r1,1));
#else // only SSE2
   SX_CHECK (n % 8 == 0, n);
   __m128d r1 = {0., 0.}, r2 = {0., 0.},
           r3 = {0., 0.}, r4 = {0., 0.};
   if (((size_t)elements & 15) == 0)  {
      /// SSE2 aligned
      const __m128d *ptr = (const __m128d*)elements;
      n = n/2;
      for (ssize_t i = 0; i < n; i+= 4)  {
         __m128d x1 = ptr[i];
         __m128d x2 = ptr[i+1];
         r1 = _mm_add_pd(r1, _mm_mul_pd(x1, x1));
         r2 = _mm_add_pd(r2, _mm_mul_pd(x2, x2));
         __m128d x3 = ptr[i + 2];
         __m128d x4 = ptr[i + 3];
         r3 = _mm_add_pd(r3, _mm_mul_pd(x3, x3));
         r4 = _mm_add_pd(r4, _mm_mul_pd(x4, x4));
      }
   } else {
      /// SSE2 unaligned (could be slower)
      for (ssize_t i = 0; i < n; i+= 4)  {
         __m128d x1 = _mm_loadu_pd(elements + i);
         __m128d x2 = _mm_loadu_pd(elements + i + 2);
         r1 = _mm_add_pd(r1, _mm_mul_pd(x1, x1));
         r2 = _mm_add_pd(r2, _mm_mul_pd(x2, x2));
         __m128d x3 = _mm_loadu_pd(elements + i + 4);
         __m128d x4 = _mm_loadu_pd(elements + i + 6);
         r3 = _mm_add_pd(r3, _mm_mul_pd(x3, x3));
         r4 = _mm_add_pd(r4, _mm_mul_pd(x4, x4));
      }
   }
   r1 = _mm_add_pd(r1,r2);
   r3 = _mm_add_pd(r3,r4);
   r1 = _mm_add_pd(r1,r3);
   return r1;
#endif
}
#endif

/// Sum squares of elements
double SxVecCompute<Simd>::normSqr (const double *elements, ssize_t n)
{
#if defined __GNUG__ && defined __AVX__
   double res = 0.;
   // --- clean-up: make n compatible with loop stride (16)
   while (n & 15)  {
      double x = *elements++;
      res += x*x;
      n--;
   }
#     ifdef USE_OPENMP
   if (n > sxChunkSize * 16)  {
      double res_ = 0.;
#     pragma omp parallel reduction(+:res_)
      {
         ssize_t n16 = n / 16;
         ssize_t nChunk = sx_omp_chunk (n16);
         ssize_t iStart = sx_omp_start (nChunk, n16);
         __m128d rx = normSqr_simd(elements + iStart * 16, nChunk * 16);
         res_ = _mm_cvtsd_f64(_mm_hadd_pd(rx, rx));
      }
      return res + res_;
   }
#     endif
   // no openmp parallelization
   __m128d rx = normSqr_simd(elements, n);
   return res + _mm_cvtsd_f64(_mm_hadd_pd(rx, rx));

#elif defined __GNUG__ && defined __SSE2__
   if (n<=0) return 0.;
   double res = 0.;
   // SSE2 alignment
   if (((size_t)elements & 15) == 8)  {
      res = (*elements) * (*elements);
      elements++;
      n--;
   }
   // now odd number of elements? => add last one later
   bool addLast = false;
   if (n & 1)  {
      addLast = true;
      n--;
   }
   // => n should be even and elements aligned to 16 bytes
   // --- clean-up: make n compatible with loop stride (8)
   while (n & 7)  {
      double x = *elements++;
      res += x*x;
      n--;
   }
#     ifdef USE_OPENMP
   if (n > sxChunkSize * 8)  {
      double res_ = 0.;
#     pragma omp parallel reduction(+:res_)
      {
         ssize_t n8 = n / 8;
         ssize_t nChunk = sx_omp_chunk (n8);
         ssize_t iStart = sx_omp_start (nChunk, n8);
         __m128d rx = normSqr_simd(elements + iStart * 8, nChunk * 8);
         res_ = _mm_cvtsd_f64(_mm_add_pd(rx, _mm_unpackhi_pd(rx,rx)));
      }
      res += res_;
   } else
#     endif
   {
      __m128d rx = normSqr_simd(elements, n);
      res += _mm_cvtsd_f64(_mm_add_pd(rx, _mm_unpackhi_pd(rx,rx)));
   }
   if (addLast) res += elements[1] * elements[1];
   return res;
#else
   double res = 0.;
#ifdef USE_OPENMP
# pragma omp parallel for reduction(+:res) if (n > sxChunkSize)
#endif
   for (ssize_t i = 0; i < n; i++)
      res += elements[i] * elements[i];
   return res;
#endif
}

// --- The next 132 lines were generated from snippets/SxVecCompute.cpp snippet isHermitian
bool SxVecCompute<Simd>::isHermitian (const SxVector<float> &H)
{
   if (H.getNRows () != H.getNCols ()) return false;

   float colLimit = 0.;
   for (ssize_t i=0; i < H.getNRows (); i++)  {
      bool haveNorm = false;
      // off-diag elements:  M_ij = M_ij
      for (ssize_t j=0; j < i; j++)  {
         float Hij = H(i,j);
         float Hji = H(j,i);
         float diff = fabsf(Hij - Hji);
         if (diff > 1e-6f)  {
            // check relative deviation of less than 1e-6f
            if (diff > 1e-6f * fabsf(Hij)) {
               // --- check relative deviation compared to column norm
               if (!haveNorm)  {
                  float colNorm = norm2(H.elements + i * H.getColStride (), H.getNRows ());
                  colLimit = 1e-6f * colNorm;
                  haveNorm = true;
               }
               if (diff > colLimit)
                  return false;
            }
         }
      }
   }
   return true;
}

bool SxVecCompute<Simd>::isHermitian (const SxVector<SxComplex8> &H)
{
   if (H.getNRows () != H.getNCols ()) return false;

   float colLimit = 0.;
   for (ssize_t i=0; i < H.getNRows (); i++)  {
      // diag must be real
      if (fabsf(H(i,i).im) > 1e-6f)  {
         if (fabsf(H(i,i).im > 1e-6f * fabsf(H(i,i).re)))
            return false;
      }
      bool haveNorm = false;
      // off-diag elements:  M_ij = M~_ij
      for (ssize_t j=0; j < i; j++)  {
         const SxComplex8 &Hij = H(i,j);
         const SxComplex8 &Hji = H(j,i);
         if (   fabsf(Hij.re - Hji.re) > 1e-6f
             || fabsf(Hij.im + Hji.im) > 1e-6f)  {
            // check relative deviation of less than 1e-12
            float diff2 = (Hij - Hji.conj ()).absSqr ();
            if (diff2 > 1e-12f * Hij.absSqr ())  {
               // --- check relative deviation compared to column norm
               if (!haveNorm)  {
                  float colNorm = norm2(H.elements + i * H.getColStride (), H.getNRows ());
                  colLimit = 1e-12f * colNorm * colNorm;
                  haveNorm = true;
               }
               if (diff2 > colLimit)
                  return false;
            }
         }
      }
   }
   return true;
}

bool SxVecCompute<Simd>::isHermitian (const SxVector<double> &H)
{
   if (H.getNRows () != H.getNCols ()) return false;

   double colLimit = 0.;
   for (ssize_t i=0; i < H.getNRows (); i++)  {
      bool haveNorm = false;
      // off-diag elements:  M_ij = M_ij
      for (ssize_t j=0; j < i; j++)  {
         double Hij = H(i,j);
         double Hji = H(j,i);
         double diff = fabs(Hij - Hji);
         if (diff > 1e-12)  {
            // check relative deviation of less than 1e-12
            if (diff > 1e-12 * fabs(Hij)) {
               // --- check relative deviation compared to column norm
               if (!haveNorm)  {
                  double colNorm = norm2(H.elements + i * H.getColStride (), H.getNRows ());
                  colLimit = 1e-12 * colNorm;
                  haveNorm = true;
               }
               if (diff > colLimit)
                  return false;
            }
         }
      }
   }
   return true;
}

bool SxVecCompute<Simd>::isHermitian (const SxVector<SxComplex16> &H)
{
   if (H.getNRows () != H.getNCols ()) return false;

   double colLimit = 0.;
   for (ssize_t i=0; i < H.getNRows (); i++)  {
      // diag must be real
      if (fabs(H(i,i).im) > 1e-12)  {
         if (fabs(H(i,i).im > 1e-12 * fabs(H(i,i).re)))
            return false;
      }
      bool haveNorm = false;
      // off-diag elements:  M_ij = M~_ij
      for (ssize_t j=0; j < i; j++)  {
         const SxComplex16 &Hij = H(i,j);
         const SxComplex16 &Hji = H(j,i);
         if (   fabs(Hij.re - Hji.re) > 1e-12
             || fabs(Hij.im + Hji.im) > 1e-12)  {
            // check relative deviation of less than 1e-12
            double diff2 = (Hij - Hji.conj ()).absSqr ();
            if (diff2 > 1e-24 * Hij.absSqr ())  {
               // --- check relative deviation compared to column norm
               if (!haveNorm)  {
                  double colNorm = norm2(H.elements + i * H.getColStride (), H.getNRows ());
                  colLimit = 1e-24 * colNorm * colNorm;
                  haveNorm = true;
               }
               if (diff2 > colLimit)
                  return false;
            }
         }
      }
   }
   return true;
}
// --- isHermitian
