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

#ifndef _SX_VEC_COMPUTE_HPP
#define _SX_VEC_COMPUTE_HPP
#include <SxCacheOblivious.h>

#ifdef USE_OPENMP
inline ssize_t sx_omp_chunk(ssize_t n)
{
   int nThread = omp_get_num_threads ();
   ssize_t nChunk = n / nThread;
   if (nChunk * nThread < n) nChunk++;
   return nChunk;
}

inline ssize_t sx_omp_start(ssize_t &nChunk, ssize_t n)
{
   int iThread = omp_get_thread_num ();
   ssize_t iStart = iThread * nChunk;
   if (nChunk + iStart > n) nChunk = n - iStart;
   return iStart;
}
#endif

// --- absSqr implementations for real and complex numbers
inline float  sxAbsSqr(const SxComplex<float> &x) { return x.absSqr (); }
inline double sxAbsSqr(const SxComplex<double> &x) { return x.absSqr (); }
template <class T>
inline T sxAbsSqr(const T x) { return x*x; }

// --- The next 232 lines were generated from snippets/SxVecCompute.hpp snippet Simd_OP
// vector x + y (x, y both contiguous, same type)
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::add (const SxVecRef<T, Layout> &x, const SxVecRef<T, Layout> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

   add ((T*)mem->mem, n, x.elements, y.elements);
   return typename SxVecResult<T,Layout>::VecType(mem, x, x.auxData, y.auxData);
}

template<class T, class T2>
void SxVecCompute<Simd>::add (decltype(T(0) + T2(0)) *resPtr, ssize_t n,
                               const T* x, const T2* y)
{
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x[i] + y[i];
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x[i] + y[i];
}

template<class T, StrideType Layout>
inline typename SxVecResult<T,Layout>::VecType
SxVecCompute<Simd>::add (const SxVecRef<T,Layout> &x, SxVecRef<T,Layout> &&y)
{
   // add in place of memory of vector y
   SX_CHECK (y.isOnlyRef ());
   SX_CHECK (x.getSize () == y.getSize (), x.getSize (), y.getSize ());
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxVecCompute<Simd>::add (y.elements, x.getSize (), x.elements, y.elements);
   SxAllocMem* alloc = y.allocMem;
   alloc->ref ();

   // --- make sure that y cannot access the data any longer
   y.unref_ ();
   y.allocMem = NULL;
   y.elements = NULL;
   // .. but keep y's layout (in case x==y) and auxData (still needed)

   return typename SxVecResult<T,Layout>::VecType(alloc, x, x.auxData, y.auxData);
}

// vector x - y (x, y both contiguous, same type)
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::subtract (const SxVecRef<T, Layout> &x, const SxVecRef<T, Layout> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

   subtract ((T*)mem->mem, n, x.elements, y.elements);
   return typename SxVecResult<T,Layout>::VecType(mem, x, x.auxData, y.auxData);
}

template<class T, class T2>
void SxVecCompute<Simd>::subtract (decltype(T(0) - T2(0)) *resPtr, ssize_t n,
                               const T* x, const T2* y)
{
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x[i] - y[i];
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x[i] - y[i];
}

template<class T, StrideType Layout>
inline typename SxVecResult<T,Layout>::VecType
SxVecCompute<Simd>::subtract (const SxVecRef<T,Layout> &x, SxVecRef<T,Layout> &&y)
{
   // subtract in place of memory of vector y
   SX_CHECK (y.isOnlyRef ());
   SX_CHECK (x.getSize () == y.getSize (), x.getSize (), y.getSize ());
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxVecCompute<Simd>::subtract (y.elements, x.getSize (), x.elements, y.elements);
   SxAllocMem* alloc = y.allocMem;
   alloc->ref ();

   // --- make sure that y cannot access the data any longer
   y.unref_ ();
   y.allocMem = NULL;
   y.elements = NULL;
   // .. but keep y's layout (in case x==y) and auxData (still needed)

   return typename SxVecResult<T,Layout>::VecType(alloc, x, x.auxData, y.auxData);
}

// vector x * y (x, y both contiguous, same type)
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::multiply (const SxVecRef<T, Layout> &x, const SxVecRef<T, Layout> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

   multiply ((T*)mem->mem, n, x.elements, y.elements);
   return typename SxVecResult<T,Layout>::VecType(mem, x, x.auxData, y.auxData);
}

template<class T, class T2>
void SxVecCompute<Simd>::multiply (decltype(T(0) * T2(0)) *resPtr, ssize_t n,
                               const T* x, const T2* y)
{
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x[i] * y[i];
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x[i] * y[i];
}

template<class T, StrideType Layout>
inline typename SxVecResult<T,Layout>::VecType
SxVecCompute<Simd>::multiply (const SxVecRef<T,Layout> &x, SxVecRef<T,Layout> &&y)
{
   // multiply in place of memory of vector y
   SX_CHECK (y.isOnlyRef ());
   SX_CHECK (x.getSize () == y.getSize (), x.getSize (), y.getSize ());
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxVecCompute<Simd>::multiply (y.elements, x.getSize (), x.elements, y.elements);
   SxAllocMem* alloc = y.allocMem;
   alloc->ref ();

   // --- make sure that y cannot access the data any longer
   y.unref_ ();
   y.allocMem = NULL;
   y.elements = NULL;
   // .. but keep y's layout (in case x==y) and auxData (still needed)

   return typename SxVecResult<T,Layout>::VecType(alloc, x, x.auxData, y.auxData);
}

// vector x / y (x, y both contiguous, same type)
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::divide (const SxVecRef<T, Layout> &x, const SxVecRef<T, Layout> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

   divide ((T*)mem->mem, n, x.elements, y.elements);
   return typename SxVecResult<T,Layout>::VecType(mem, x, x.auxData, y.auxData);
}

template<class T, class T2>
void SxVecCompute<Simd>::divide (decltype(T(0) / T2(0)) *resPtr, ssize_t n,
                               const T* x, const T2* y)
{
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x[i] / y[i];
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x[i] / y[i];
}

template<class T, StrideType Layout>
inline typename SxVecResult<T,Layout>::VecType
SxVecCompute<Simd>::divide (const SxVecRef<T,Layout> &x, SxVecRef<T,Layout> &&y)
{
   // divide in place of memory of vector y
   SX_CHECK (y.isOnlyRef ());
   SX_CHECK (x.getSize () == y.getSize (), x.getSize (), y.getSize ());
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxVecCompute<Simd>::divide (y.elements, x.getSize (), x.elements, y.elements);
   SxAllocMem* alloc = y.allocMem;
   alloc->ref ();

   // --- make sure that y cannot access the data any longer
   y.unref_ ();
   y.allocMem = NULL;
   y.elements = NULL;
   // .. but keep y's layout (in case x==y) and auxData (still needed)

   return typename SxVecResult<T,Layout>::VecType(alloc, x, x.auxData, y.auxData);
}
// --- Simd_OP

template<class T1, class T2>
inline
void SxVecCompute<Simd>::multiplyMx1 (decltype(T1(0) * T2(0)) *resPtr,
                                      /* ssize_t colStrideRes, */
                                      ssize_t n, ssize_t m,
                                      const T1* x,
                                      ssize_t colStrideX,
                                      const T2* y)
{
   typedef decltype(T1(0) * T2(0)) TRes;
   ssize_t colStrideRes = n;
#ifdef USE_OPENMP
#  pragma omp parallel for if (n > sxChunkSize)
#endif
   // chunking the rows (to keep y chunk in cache)
   for (ssize_t i0 = 0; i0 < n; i0 += 128)  {
      ssize_t n1 = n - i0;
      if (n1 > 128) n1 = 128;
      const T2 *y1 = y + i0;
      // chunking the columns (to economize y loads)
      for (ssize_t j0 = 0; j0 < m; j0 += 4)  {
         // hardcoded loop implementations for j1 = 1, 2, 3, 4
         switch ( (m - j0) > 4 ? 4 : (m -j0) )  {
            case 1:
               {
                  TRes *dest1 = resPtr + i0 + j0 * colStrideRes;
                  const T1 *x1 = x + i0 + j0 * colStrideX;
#ifdef USE_OPENMP
#                 pragma omp simd
#endif
                  for (ssize_t i1 = 0; i1 < n1; ++i1)
                     dest1[i1] = x1[i1] * y1[i1];
                  break;
               }
            case 2:
               {
                  TRes *dest1 = resPtr + i0 + j0       * colStrideRes,
                       *dest2 = resPtr + i0 + (j0 + 1) * colStrideRes;
                  const T1 *x1 = x + i0 + j0       * colStrideX,
                           *x2 = x + i0 + (j0 + 1) * colStrideX;
#ifdef USE_OPENMP
#                 pragma omp simd
#endif
                  for (ssize_t i1 = 0; i1 < n1; ++i1)  {
                     dest1[i1] = x1[i1] * y1[i1];
                     dest2[i1] = x2[i1] * y1[i1];
                  }
                  break;
               }
            case 3:
               {
                  TRes *dest1 = resPtr + i0 + j0       * colStrideRes,
                       *dest2 = resPtr + i0 + (j0 + 1) * colStrideRes,
                       *dest3 = resPtr + i0 + (j0 + 2) * colStrideRes;
                  const T1 *x1 = x + i0 + j0       * colStrideX,
                           *x2 = x + i0 + (j0 + 1) * colStrideX,
                           *x3 = x + i0 + (j0 + 2) * colStrideX;
#ifdef USE_OPENMP
#                 pragma omp simd
#endif
                  for (ssize_t i1 = 0; i1 < n1; ++i1)  {
                     dest1[i1] = x1[i1] * y1[i1];
                     dest2[i1] = x2[i1] * y1[i1];
                     dest3[i1] = x3[i1] * y1[i1];
                  }
                  break;
               }
            case 4:
               {
                  TRes *dest1 = resPtr + i0 + j0       * colStrideRes,
                       *dest2 = resPtr + i0 + (j0 + 1) * colStrideRes,
                       *dest3 = resPtr + i0 + (j0 + 2) * colStrideRes,
                       *dest4 = resPtr + i0 + (j0 + 3) * colStrideRes;
                  const T1 *x1 = x + i0 + j0       * colStrideX,
                           *x2 = x + i0 + (j0 + 1) * colStrideX,
                           *x3 = x + i0 + (j0 + 2) * colStrideX,
                           *x4 = x + i0 + (j0 + 3) * colStrideX;
#ifdef USE_OPENMP
#                 pragma omp simd
#endif
                  for (ssize_t i1 = 0; i1 < n1; ++i1)  {
                     dest1[i1] = x1[i1] * y1[i1];
                     dest2[i1] = x2[i1] * y1[i1];
                     dest3[i1] = x3[i1] * y1[i1];
                     dest4[i1] = x4[i1] * y1[i1];
                  }
                  break;
               }
         }
      }
   }
}

// --- The next 100 lines were generated from snippets/SxVecCompute.hpp snippet Simd_OPEQ
// vector x+= y (x, y both contiguous, same type)
template<class T, StrideType Layout>
void
SxVecCompute<Simd>::addInPlace (SxVecRef<T, Layout> &x, const SxVecRef<T, Layout> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         x.elements[i] += y.elements[i];
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         x.elements[i] += y.elements[i];
   SX_VALIDATE_VECTOR(x);
}

// vector x-= y (x, y both contiguous, same type)
template<class T, StrideType Layout>
void
SxVecCompute<Simd>::subtractInPlace (SxVecRef<T, Layout> &x, const SxVecRef<T, Layout> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         x.elements[i] -= y.elements[i];
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         x.elements[i] -= y.elements[i];
   SX_VALIDATE_VECTOR(x);
}

// vector x*= y (x, y both contiguous, same type)
template<class T, StrideType Layout>
void
SxVecCompute<Simd>::multiplyInPlace (SxVecRef<T, Layout> &x, const SxVecRef<T, Layout> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         x.elements[i] *= y.elements[i];
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         x.elements[i] *= y.elements[i];
   SX_VALIDATE_VECTOR(x);
}

// vector x/= y (x, y both contiguous, same type)
template<class T, StrideType Layout>
void
SxVecCompute<Simd>::divideInPlace (SxVecRef<T, Layout> &x, const SxVecRef<T, Layout> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         x.elements[i] /= y.elements[i];
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         x.elements[i] /= y.elements[i];
   SX_VALIDATE_VECTOR(x);
}
// --- Simd_OPEQ

// --- The next 256 lines were generated from snippets/SxVecCompute.hpp snippet SimdCol_OP
template<class T, StrideType Layout, StrideType Layout2>
SxVector<T>
SxVecCompute<SimdCol>::add (const SxVecRef<T, Layout> &x, const SxVecRef<T, Layout2> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SX_CHECK (x.getRowStride () == 1);
   SX_CHECK (y.getRowStride () == 1);
   ssize_t n = x.getSize ();
   ssize_t nRows = x.getNRows ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t lenChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (lenChunk, n);
      ssize_t iEnd = iStart + lenChunk - 1;
      ssize_t iCol = iStart / nRows;
      ssize_t iRow = iStart - iCol * nRows;
      ssize_t iCEnd = iEnd / nRows;
      T* xPtr = x.elements + iCol * x.getColStride ();
      T* yPtr = y.elements + iCol * y.getColStride ();
      T* resPtr = (T*)mem->mem + iStart - iRow;
      for ( ; iCol < iCEnd ; iCol++)  {
         #pragma omp simd
         for (ssize_t i = iRow ; i < nRows; i++)
            resPtr[i] = xPtr[i] + yPtr[i];
         iRow = 0;
         xPtr += x.getColStride ();
         yPtr += y.getColStride ();
         resPtr += nRows;
      }
      // last column
      ssize_t rLast = iEnd - iCEnd * nRows;
      #pragma omp simd
      for (ssize_t i = 0 ; i <= rLast; i++)
         resPtr[i] = xPtr[i] + yPtr[i];
   } else
#endif
   {
      ssize_t nCols = x.getNCols ();
      T* xPtr = x.elements;
      T* yPtr = y.elements;
      T* resPtr = (T*)mem->mem;
      for (ssize_t c = 0 ; c < nCols; c++)  {
#ifdef USE_OPENMP
         #pragma omp simd
#endif
         for (ssize_t r = 0; r < nRows; r++)
            resPtr[r] = xPtr[r] + yPtr[r];
         xPtr += x.getColStride ();
         yPtr += y.getColStride ();
         resPtr += nRows;
      }
   }
   return SxVector<T>(mem, x, x.auxData, y.auxData);
}

template<class T, StrideType Layout, StrideType Layout2>
SxVector<T>
SxVecCompute<SimdCol>::subtract (const SxVecRef<T, Layout> &x, const SxVecRef<T, Layout2> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SX_CHECK (x.getRowStride () == 1);
   SX_CHECK (y.getRowStride () == 1);
   ssize_t n = x.getSize ();
   ssize_t nRows = x.getNRows ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t lenChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (lenChunk, n);
      ssize_t iEnd = iStart + lenChunk - 1;
      ssize_t iCol = iStart / nRows;
      ssize_t iRow = iStart - iCol * nRows;
      ssize_t iCEnd = iEnd / nRows;
      T* xPtr = x.elements + iCol * x.getColStride ();
      T* yPtr = y.elements + iCol * y.getColStride ();
      T* resPtr = (T*)mem->mem + iStart - iRow;
      for ( ; iCol < iCEnd ; iCol++)  {
         #pragma omp simd
         for (ssize_t i = iRow ; i < nRows; i++)
            resPtr[i] = xPtr[i] - yPtr[i];
         iRow = 0;
         xPtr += x.getColStride ();
         yPtr += y.getColStride ();
         resPtr += nRows;
      }
      // last column
      ssize_t rLast = iEnd - iCEnd * nRows;
      #pragma omp simd
      for (ssize_t i = 0 ; i <= rLast; i++)
         resPtr[i] = xPtr[i] - yPtr[i];
   } else
#endif
   {
      ssize_t nCols = x.getNCols ();
      T* xPtr = x.elements;
      T* yPtr = y.elements;
      T* resPtr = (T*)mem->mem;
      for (ssize_t c = 0 ; c < nCols; c++)  {
#ifdef USE_OPENMP
         #pragma omp simd
#endif
         for (ssize_t r = 0; r < nRows; r++)
            resPtr[r] = xPtr[r] - yPtr[r];
         xPtr += x.getColStride ();
         yPtr += y.getColStride ();
         resPtr += nRows;
      }
   }
   return SxVector<T>(mem, x, x.auxData, y.auxData);
}

template<class T, StrideType Layout, StrideType Layout2>
SxVector<T>
SxVecCompute<SimdCol>::multiply (const SxVecRef<T, Layout> &x, const SxVecRef<T, Layout2> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SX_CHECK (x.getRowStride () == 1);
   SX_CHECK (y.getRowStride () == 1);
   ssize_t n = x.getSize ();
   ssize_t nRows = x.getNRows ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t lenChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (lenChunk, n);
      ssize_t iEnd = iStart + lenChunk - 1;
      ssize_t iCol = iStart / nRows;
      ssize_t iRow = iStart - iCol * nRows;
      ssize_t iCEnd = iEnd / nRows;
      T* xPtr = x.elements + iCol * x.getColStride ();
      T* yPtr = y.elements + iCol * y.getColStride ();
      T* resPtr = (T*)mem->mem + iStart - iRow;
      for ( ; iCol < iCEnd ; iCol++)  {
         #pragma omp simd
         for (ssize_t i = iRow ; i < nRows; i++)
            resPtr[i] = xPtr[i] * yPtr[i];
         iRow = 0;
         xPtr += x.getColStride ();
         yPtr += y.getColStride ();
         resPtr += nRows;
      }
      // last column
      ssize_t rLast = iEnd - iCEnd * nRows;
      #pragma omp simd
      for (ssize_t i = 0 ; i <= rLast; i++)
         resPtr[i] = xPtr[i] * yPtr[i];
   } else
#endif
   {
      ssize_t nCols = x.getNCols ();
      T* xPtr = x.elements;
      T* yPtr = y.elements;
      T* resPtr = (T*)mem->mem;
      for (ssize_t c = 0 ; c < nCols; c++)  {
#ifdef USE_OPENMP
         #pragma omp simd
#endif
         for (ssize_t r = 0; r < nRows; r++)
            resPtr[r] = xPtr[r] * yPtr[r];
         xPtr += x.getColStride ();
         yPtr += y.getColStride ();
         resPtr += nRows;
      }
   }
   return SxVector<T>(mem, x, x.auxData, y.auxData);
}

template<class T, StrideType Layout, StrideType Layout2>
SxVector<T>
SxVecCompute<SimdCol>::divide (const SxVecRef<T, Layout> &x, const SxVecRef<T, Layout2> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SX_CHECK (x.getRowStride () == 1);
   SX_CHECK (y.getRowStride () == 1);
   ssize_t n = x.getSize ();
   ssize_t nRows = x.getNRows ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t lenChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (lenChunk, n);
      ssize_t iEnd = iStart + lenChunk - 1;
      ssize_t iCol = iStart / nRows;
      ssize_t iRow = iStart - iCol * nRows;
      ssize_t iCEnd = iEnd / nRows;
      T* xPtr = x.elements + iCol * x.getColStride ();
      T* yPtr = y.elements + iCol * y.getColStride ();
      T* resPtr = (T*)mem->mem + iStart - iRow;
      for ( ; iCol < iCEnd ; iCol++)  {
         #pragma omp simd
         for (ssize_t i = iRow ; i < nRows; i++)
            resPtr[i] = xPtr[i] / yPtr[i];
         iRow = 0;
         xPtr += x.getColStride ();
         yPtr += y.getColStride ();
         resPtr += nRows;
      }
      // last column
      ssize_t rLast = iEnd - iCEnd * nRows;
      #pragma omp simd
      for (ssize_t i = 0 ; i <= rLast; i++)
         resPtr[i] = xPtr[i] / yPtr[i];
   } else
#endif
   {
      ssize_t nCols = x.getNCols ();
      T* xPtr = x.elements;
      T* yPtr = y.elements;
      T* resPtr = (T*)mem->mem;
      for (ssize_t c = 0 ; c < nCols; c++)  {
#ifdef USE_OPENMP
         #pragma omp simd
#endif
         for (ssize_t r = 0; r < nRows; r++)
            resPtr[r] = xPtr[r] / yPtr[r];
         xPtr += x.getColStride ();
         yPtr += y.getColStride ();
         resPtr += nRows;
      }
   }
   return SxVector<T>(mem, x, x.auxData, y.auxData);
}
// --- SimdCol_OP

// --- The next 464 lines were generated from snippets/SxVecCompute.hpp snippet SimdIt_OP
template<class T, class T2, StrideType Layout>
inline void
SxVecCompute<SimdIt>::add (const T* xPtr, SxVecConstIt<T2, Layout> it,
                            ssize_t n, decltype(T() + T2()) *resPtr)
{
   if (n & 1)  {
      *resPtr++ = *xPtr++ + *it;
      ++it;
      n--;
   }
#  ifdef USE_OPENMP
#  pragma omp simd
#  endif
   for (ssize_t i = 0; i < n; i+=2) {
      T2 y1 = *it;
      ++it;
      T2 y2 = *it;
      ++it;
      resPtr[0] = xPtr[0] + y1;
      resPtr[1] = xPtr[1] + y2;
      resPtr += 2;
      xPtr += 2;
   }
}

// vector x + y (x contiguous, y is strided, same type)
template<class T, StrideType Layout>
SxVector<T>
SxVecCompute<SimdIt>::add (const SxVecRef<T, Compact> &x,
                            const SxVecRef<T, Layout>  &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      if (nChunk & 1) nChunk++;
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> it = y.begin ();
      it += iStart;
      add (x.elements + iStart, it, nChunk, (T*)mem->mem + iStart);
   } else
#endif
   {
      add (x.elements, y.begin (), n, (T*)mem->mem);
   }
   return SxVector<T>(mem, x, x.auxData, y.auxData);
}

template<class T, class T2, StrideType Layout>
inline void
SxVecCompute<SimdIt>::add (SxVecConstIt<T, Layout> it, const T2* xPtr,
                            ssize_t n, decltype(T() + T2()) *resPtr)
{
   if (n & 1)  {
      *resPtr++ =  *it + *xPtr++;
      ++it;
      n--;
   }
#  ifdef USE_OPENMP
#  pragma omp simd
#  endif
   for (ssize_t i = 0; i < n; i+=2) {
      T y1 = *it;
      ++it;
      T y2 = *it;
      ++it;
      resPtr[0] = y1 + xPtr[0];
      resPtr[1] = y2 + xPtr[1];
      resPtr += 2;
      xPtr += 2;
   }
}

// vector x + y (x strided, y is contiguous, same type)
template<class T, StrideType Layout>
SxVector<T>
SxVecCompute<SimdIt>::add (const SxVecRef<T, Layout> &x,
                            const SxVecRef<T, Compact> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      if (nChunk & 1) nChunk++;
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> it = x.begin ();
      it += iStart;
      add (it, x.elements + iStart, nChunk, (T*)mem->mem + iStart);
   } else
#endif
   {
      add (x.begin (), y.elements, n, (T*)mem->mem);
   }
   return SxVector<T>(mem, x, x.auxData, y.auxData);
}

template<class T, class T2, StrideType Layout>
inline void
SxVecCompute<SimdIt>::subtract (const T* xPtr, SxVecConstIt<T2, Layout> it,
                            ssize_t n, decltype(T() - T2()) *resPtr)
{
   if (n & 1)  {
      *resPtr++ = *xPtr++ - *it;
      ++it;
      n--;
   }
#  ifdef USE_OPENMP
#  pragma omp simd
#  endif
   for (ssize_t i = 0; i < n; i+=2) {
      T2 y1 = *it;
      ++it;
      T2 y2 = *it;
      ++it;
      resPtr[0] = xPtr[0] - y1;
      resPtr[1] = xPtr[1] - y2;
      resPtr += 2;
      xPtr += 2;
   }
}

// vector x + y (x contiguous, y is strided, same type)
template<class T, StrideType Layout>
SxVector<T>
SxVecCompute<SimdIt>::subtract (const SxVecRef<T, Compact> &x,
                            const SxVecRef<T, Layout>  &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      if (nChunk & 1) nChunk++;
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> it = y.begin ();
      it += iStart;
      subtract (x.elements + iStart, it, nChunk, (T*)mem->mem + iStart);
   } else
#endif
   {
      subtract (x.elements, y.begin (), n, (T*)mem->mem);
   }
   return SxVector<T>(mem, x, x.auxData, y.auxData);
}

template<class T, class T2, StrideType Layout>
inline void
SxVecCompute<SimdIt>::subtract (SxVecConstIt<T, Layout> it, const T2* xPtr,
                            ssize_t n, decltype(T() - T2()) *resPtr)
{
   if (n & 1)  {
      *resPtr++ =  *it - *xPtr++;
      ++it;
      n--;
   }
#  ifdef USE_OPENMP
#  pragma omp simd
#  endif
   for (ssize_t i = 0; i < n; i+=2) {
      T y1 = *it;
      ++it;
      T y2 = *it;
      ++it;
      resPtr[0] = y1 - xPtr[0];
      resPtr[1] = y2 - xPtr[1];
      resPtr += 2;
      xPtr += 2;
   }
}

// vector x + y (x strided, y is contiguous, same type)
template<class T, StrideType Layout>
SxVector<T>
SxVecCompute<SimdIt>::subtract (const SxVecRef<T, Layout> &x,
                            const SxVecRef<T, Compact> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      if (nChunk & 1) nChunk++;
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> it = x.begin ();
      it += iStart;
      subtract (it, x.elements + iStart, nChunk, (T*)mem->mem + iStart);
   } else
#endif
   {
      subtract (x.begin (), y.elements, n, (T*)mem->mem);
   }
   return SxVector<T>(mem, x, x.auxData, y.auxData);
}

template<class T, class T2, StrideType Layout>
inline void
SxVecCompute<SimdIt>::multiply (const T* xPtr, SxVecConstIt<T2, Layout> it,
                            ssize_t n, decltype(T() * T2()) *resPtr)
{
   if (n & 1)  {
      *resPtr++ = *xPtr++ * *it;
      ++it;
      n--;
   }
#  ifdef USE_OPENMP
#  pragma omp simd
#  endif
   for (ssize_t i = 0; i < n; i+=2) {
      T2 y1 = *it;
      ++it;
      T2 y2 = *it;
      ++it;
      resPtr[0] = xPtr[0] * y1;
      resPtr[1] = xPtr[1] * y2;
      resPtr += 2;
      xPtr += 2;
   }
}

// vector x + y (x contiguous, y is strided, same type)
template<class T, StrideType Layout>
SxVector<T>
SxVecCompute<SimdIt>::multiply (const SxVecRef<T, Compact> &x,
                            const SxVecRef<T, Layout>  &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      if (nChunk & 1) nChunk++;
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> it = y.begin ();
      it += iStart;
      multiply (x.elements + iStart, it, nChunk, (T*)mem->mem + iStart);
   } else
#endif
   {
      multiply (x.elements, y.begin (), n, (T*)mem->mem);
   }
   return SxVector<T>(mem, x, x.auxData, y.auxData);
}

template<class T, class T2, StrideType Layout>
inline void
SxVecCompute<SimdIt>::multiply (SxVecConstIt<T, Layout> it, const T2* xPtr,
                            ssize_t n, decltype(T() * T2()) *resPtr)
{
   if (n & 1)  {
      *resPtr++ =  *it * *xPtr++;
      ++it;
      n--;
   }
#  ifdef USE_OPENMP
#  pragma omp simd
#  endif
   for (ssize_t i = 0; i < n; i+=2) {
      T y1 = *it;
      ++it;
      T y2 = *it;
      ++it;
      resPtr[0] = y1 * xPtr[0];
      resPtr[1] = y2 * xPtr[1];
      resPtr += 2;
      xPtr += 2;
   }
}

// vector x + y (x strided, y is contiguous, same type)
template<class T, StrideType Layout>
SxVector<T>
SxVecCompute<SimdIt>::multiply (const SxVecRef<T, Layout> &x,
                            const SxVecRef<T, Compact> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      if (nChunk & 1) nChunk++;
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> it = x.begin ();
      it += iStart;
      multiply (it, x.elements + iStart, nChunk, (T*)mem->mem + iStart);
   } else
#endif
   {
      multiply (x.begin (), y.elements, n, (T*)mem->mem);
   }
   return SxVector<T>(mem, x, x.auxData, y.auxData);
}

template<class T, class T2, StrideType Layout>
inline void
SxVecCompute<SimdIt>::divide (const T* xPtr, SxVecConstIt<T2, Layout> it,
                            ssize_t n, decltype(T() / T2()) *resPtr)
{
   if (n & 1)  {
      *resPtr++ = *xPtr++ / *it;
      ++it;
      n--;
   }
#  ifdef USE_OPENMP
#  pragma omp simd
#  endif
   for (ssize_t i = 0; i < n; i+=2) {
      T2 y1 = *it;
      ++it;
      T2 y2 = *it;
      ++it;
      resPtr[0] = xPtr[0] / y1;
      resPtr[1] = xPtr[1] / y2;
      resPtr += 2;
      xPtr += 2;
   }
}

// vector x + y (x contiguous, y is strided, same type)
template<class T, StrideType Layout>
SxVector<T>
SxVecCompute<SimdIt>::divide (const SxVecRef<T, Compact> &x,
                            const SxVecRef<T, Layout>  &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      if (nChunk & 1) nChunk++;
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> it = y.begin ();
      it += iStart;
      divide (x.elements + iStart, it, nChunk, (T*)mem->mem + iStart);
   } else
#endif
   {
      divide (x.elements, y.begin (), n, (T*)mem->mem);
   }
   return SxVector<T>(mem, x, x.auxData, y.auxData);
}

template<class T, class T2, StrideType Layout>
inline void
SxVecCompute<SimdIt>::divide (SxVecConstIt<T, Layout> it, const T2* xPtr,
                            ssize_t n, decltype(T() / T2()) *resPtr)
{
   if (n & 1)  {
      *resPtr++ =  *it / *xPtr++;
      ++it;
      n--;
   }
#  ifdef USE_OPENMP
#  pragma omp simd
#  endif
   for (ssize_t i = 0; i < n; i+=2) {
      T y1 = *it;
      ++it;
      T y2 = *it;
      ++it;
      resPtr[0] = y1 / xPtr[0];
      resPtr[1] = y2 / xPtr[1];
      resPtr += 2;
      xPtr += 2;
   }
}

// vector x + y (x strided, y is contiguous, same type)
template<class T, StrideType Layout>
SxVector<T>
SxVecCompute<SimdIt>::divide (const SxVecRef<T, Layout> &x,
                            const SxVecRef<T, Compact> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      if (nChunk & 1) nChunk++;
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> it = x.begin ();
      it += iStart;
      divide (it, x.elements + iStart, nChunk, (T*)mem->mem + iStart);
   } else
#endif
   {
      divide (x.begin (), y.elements, n, (T*)mem->mem);
   }
   return SxVector<T>(mem, x, x.auxData, y.auxData);
}
// --- SimdIt_OP

// --- The next 164 lines were generated from snippets/SxVecCompute.hpp snippet UseIterator_OP
// vector x + y (x, y can be any compatible type)
template<class T, StrideType Layout, class T2, StrideType Layout2>
SxVector<decltype(T() + T2())>
SxVecCompute<UseIterator>::add (const SxVecRef<T, Layout>  &x,
                                const SxVecRef<T2, Layout2> &y)
{
   typedef decltype(T() + T2()) TRes;
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(TRes) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout>   xIt = x.begin ();
      SxVecConstIt<T2,Layout2> yIt = y.begin ();
      xIt += iStart;
      yIt += iStart;
      TRes* resPtr = (TRes*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt, ++yIt)
         resPtr[i] = *xIt + *yIt;
   } else
#endif
   {
      SxVecConstIt<T,Layout>   xIt = x.begin ();
      SxVecConstIt<T2,Layout2> yIt = y.begin ();
      TRes* resPtr = (TRes*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++xIt, ++yIt)
         resPtr[i] = *xIt + *yIt;
   }
   return SxVector<TRes>(mem, x, x.auxData, y.auxData);
}

// vector x - y (x, y can be any compatible type)
template<class T, StrideType Layout, class T2, StrideType Layout2>
SxVector<decltype(T() + T2())>
SxVecCompute<UseIterator>::subtract (const SxVecRef<T, Layout>  &x,
                                const SxVecRef<T2, Layout2> &y)
{
   typedef decltype(T() + T2()) TRes;
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(TRes) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout>   xIt = x.begin ();
      SxVecConstIt<T2,Layout2> yIt = y.begin ();
      xIt += iStart;
      yIt += iStart;
      TRes* resPtr = (TRes*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt, ++yIt)
         resPtr[i] = *xIt - *yIt;
   } else
#endif
   {
      SxVecConstIt<T,Layout>   xIt = x.begin ();
      SxVecConstIt<T2,Layout2> yIt = y.begin ();
      TRes* resPtr = (TRes*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++xIt, ++yIt)
         resPtr[i] = *xIt - *yIt;
   }
   return SxVector<TRes>(mem, x, x.auxData, y.auxData);
}

// vector x * y (x, y can be any compatible type)
template<class T, StrideType Layout, class T2, StrideType Layout2>
SxVector<decltype(T() + T2())>
SxVecCompute<UseIterator>::multiply (const SxVecRef<T, Layout>  &x,
                                const SxVecRef<T2, Layout2> &y)
{
   typedef decltype(T() + T2()) TRes;
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(TRes) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout>   xIt = x.begin ();
      SxVecConstIt<T2,Layout2> yIt = y.begin ();
      xIt += iStart;
      yIt += iStart;
      TRes* resPtr = (TRes*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt, ++yIt)
         resPtr[i] = *xIt * *yIt;
   } else
#endif
   {
      SxVecConstIt<T,Layout>   xIt = x.begin ();
      SxVecConstIt<T2,Layout2> yIt = y.begin ();
      TRes* resPtr = (TRes*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++xIt, ++yIt)
         resPtr[i] = *xIt * *yIt;
   }
   return SxVector<TRes>(mem, x, x.auxData, y.auxData);
}

// vector x / y (x, y can be any compatible type)
template<class T, StrideType Layout, class T2, StrideType Layout2>
SxVector<decltype(T() + T2())>
SxVecCompute<UseIterator>::divide (const SxVecRef<T, Layout>  &x,
                                const SxVecRef<T2, Layout2> &y)
{
   typedef decltype(T() + T2()) TRes;
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   ssize_t n = x.getSize ();
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(TRes) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout>   xIt = x.begin ();
      SxVecConstIt<T2,Layout2> yIt = y.begin ();
      xIt += iStart;
      yIt += iStart;
      TRes* resPtr = (TRes*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt, ++yIt)
         resPtr[i] = *xIt / *yIt;
   } else
#endif
   {
      SxVecConstIt<T,Layout>   xIt = x.begin ();
      SxVecConstIt<T2,Layout2> yIt = y.begin ();
      TRes* resPtr = (TRes*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++xIt, ++yIt)
         resPtr[i] = *xIt / *yIt;
   }
   return SxVector<TRes>(mem, x, x.auxData, y.auxData);
}
// --- UseIterator_OP

// --- The next 156 lines were generated from snippets/SxVecCompute.hpp snippet UseIterator_OPEQ
template<class T, StrideType Layout, class T2, StrideType Layout2>
void
SxVecCompute<UseIterator>::addInPlace (      SxVecRef<T, Layout> &x,
                                       const SxVecRef<T2, Layout2> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());

#ifdef USE_OPENMP
   ssize_t n = x.getSize ();
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecIt<T,Layout>       xIt = x.begin ();
      SxVecConstIt<T2,Layout2> yIt = y.begin ();
      xIt += iStart;
      yIt += iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt, ++yIt)  {
         *xIt += *yIt;
         SX_CHECK_NUM (*xIt);
      }
   } else
#endif
   {
      SxVecIt<T,Layout>       xIt = x.begin ();
      SxVecConstIt<T2,Layout2> yIt = y.begin ();
      for ( ; xIt != x.end (); ++xIt, ++yIt) {
         *xIt += *yIt;
         SX_CHECK_NUM (*xIt);
      }
   }
}

template<class T, StrideType Layout, class T2, StrideType Layout2>
void
SxVecCompute<UseIterator>::subtractInPlace (      SxVecRef<T, Layout> &x,
                                       const SxVecRef<T2, Layout2> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());

#ifdef USE_OPENMP
   ssize_t n = x.getSize ();
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecIt<T,Layout>       xIt = x.begin ();
      SxVecConstIt<T2,Layout2> yIt = y.begin ();
      xIt += iStart;
      yIt += iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt, ++yIt)  {
         *xIt -= *yIt;
         SX_CHECK_NUM (*xIt);
      }
   } else
#endif
   {
      SxVecIt<T,Layout>       xIt = x.begin ();
      SxVecConstIt<T2,Layout2> yIt = y.begin ();
      for ( ; xIt != x.end (); ++xIt, ++yIt) {
         *xIt -= *yIt;
         SX_CHECK_NUM (*xIt);
      }
   }
}

template<class T, StrideType Layout, class T2, StrideType Layout2>
void
SxVecCompute<UseIterator>::multiplyInPlace (      SxVecRef<T, Layout> &x,
                                       const SxVecRef<T2, Layout2> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());

#ifdef USE_OPENMP
   ssize_t n = x.getSize ();
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecIt<T,Layout>       xIt = x.begin ();
      SxVecConstIt<T2,Layout2> yIt = y.begin ();
      xIt += iStart;
      yIt += iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt, ++yIt)  {
         *xIt *= *yIt;
         SX_CHECK_NUM (*xIt);
      }
   } else
#endif
   {
      SxVecIt<T,Layout>       xIt = x.begin ();
      SxVecConstIt<T2,Layout2> yIt = y.begin ();
      for ( ; xIt != x.end (); ++xIt, ++yIt) {
         *xIt *= *yIt;
         SX_CHECK_NUM (*xIt);
      }
   }
}

template<class T, StrideType Layout, class T2, StrideType Layout2>
void
SxVecCompute<UseIterator>::divideInPlace (      SxVecRef<T, Layout> &x,
                                       const SxVecRef<T2, Layout2> &y)
{
   SX_CHECK(x.getSize () == y.getSize (),
            x.getSize (), y.getSize ());
   SX_CHECK (x.getNRows () == y.getNRows (),
             x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (),
             x.getNCols (), y.getNCols ());

#ifdef USE_OPENMP
   ssize_t n = x.getSize ();
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecIt<T,Layout>       xIt = x.begin ();
      SxVecConstIt<T2,Layout2> yIt = y.begin ();
      xIt += iStart;
      yIt += iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt, ++yIt)  {
         *xIt /= *yIt;
         SX_CHECK_NUM (*xIt);
      }
   } else
#endif
   {
      SxVecIt<T,Layout>       xIt = x.begin ();
      SxVecConstIt<T2,Layout2> yIt = y.begin ();
      for ( ; xIt != x.end (); ++xIt, ++yIt) {
         *xIt /= *yIt;
         SX_CHECK_NUM (*xIt);
      }
   }
}
// --- UseIterator_OPEQ

// y += a * x
template<class T, StrideType Layout, StrideType Layout2>
void SxVecCompute<UseIterator>::axpy(SxVecRef<T,Layout> &y, const T &a,
                                     const SxVecRef<T, Layout2> &x)
{
   SX_CHECK (y.getNCols () == x.getNCols (), y.getNCols (), x.getNCols ());
   // no support for strides beyond "int" range
   SX_CHECK (int(y.getRowStride ()) == y.getRowStride (), y.getRowStride ());
   SX_CHECK (int(x.getRowStride ()) == x.getRowStride (), x.getRowStride ());
   ssize_t nRow = x.getNRows ();
   if (   nRow * y.getRowStride () == y.getColStride ()
       && nRow * x.getRowStride () == x.getColStride ())
   {
      // arbitrary row stride
      ::axpy (y.elements, int(y.getRowStride ()), a,
              x.elements, int(x.getRowStride ()), x.getSize ());
   } else {
      // columnwise with arbitrary row/col strides
      ssize_t nCol = x.getNCols ();
#     ifdef USE_OPENMP
#     pragma omp parallel for \
         if (x.getSize () > sxChunkSize && nCol > omp_get_max_threads ())
#     endif
      for (ssize_t ic = 0; ic < nCol; ++ic)  {
         ::axpy (y.elements + ic * y.getColStride (), int(y.getRowStride ()), a,
                 x.elements + ic * x.getColStride (), int(x.getRowStride ()), nRow);
      }
   }
}

// y += a * x
template<class T, StrideType Layout, class Ta, class T2, StrideType Layout2>
void SxVecCompute<UseIterator>::axpy(SxVecRef<T,Layout> &y, const Ta &a,
                                     const SxVecRef<T2, Layout2> &x)
{
   // no implicit down conversion!
   SX_CHECK((std::is_same<T, decltype(T(0) + Ta(0) * T2(0))>::value));
   SX_CHECK(x.getSize () == y.getSize (), x.getSize (), y.getSize ());
   SX_CHECK (x.getNRows () == y.getNRows (), x.getNRows (), y.getNRows ());
   SX_CHECK (x.getNCols () == y.getNCols (), x.getNCols (), y.getNCols ());

#ifdef USE_OPENMP
   ssize_t n = x.getSize ();
   if (n > sxChunkSize)
#     pragma omp parallel if (n > sxChunkSize)
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecIt<T,Layout>        yIt = y.begin ();
      SxVecConstIt<T2,Layout2> xIt = x.begin ();
      xIt += iStart;
      yIt += iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt, ++yIt)  {
         *yIt += a * *xIt;
         SX_CHECK_NUM (*yIt);
      }
   } else
#endif
   {
      SxVecIt<T,Layout>        yIt = y.begin ();
      SxVecConstIt<T2,Layout2> xIt = x.begin ();
      for ( ; xIt != x.end (); ++xIt, ++yIt) {
         *yIt += a * *xIt;
         SX_CHECK_NUM (*yIt);
      }
   }
}

// --- The next 180 lines were generated from snippets/SxVecCompute.hpp snippet Simd_OP_scalar
// vector + scalar (vector contiguous, same type)
template<class T, StrideType Layout, class T2>
typename SxVecResult<decltype(T(0) + T2(0)), Layout>::VecType
SxVecCompute<Simd>::add (const SxVecRef<T, Layout> &x, const T2 &y)
{
   typedef decltype(T() + T2()) TRes;
   SX_CHECK (x.getSize () > 0);
   ssize_t n = x.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(TRes) * n);

   TRes* resPtr = (TRes*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x.elements[i] + y;
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x.elements[i] + y;
   return typename SxVecResult<TRes,Layout>::VecType(mem, x, x.auxData, x.auxData);
}

// vector + scalar in place (vector contiguous, same type)
template<class T, StrideType Layout, class T2>
void
SxVecCompute<Simd>::addInPlace (SxVecRef<T, Layout> &x, const T2 &y)
{
   SX_CHECK (x.getSize () > 0);
   ssize_t n = x.getSize ();

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         x.elements[i] +=  y;
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         x.elements[i] +=  y;
   SX_VALIDATE_VECTOR(x);
}

// vector - scalar (vector contiguous, same type)
template<class T, StrideType Layout, class T2>
typename SxVecResult<decltype(T(0) + T2(0)), Layout>::VecType
SxVecCompute<Simd>::subtract (const SxVecRef<T, Layout> &x, const T2 &y)
{
   typedef decltype(T() + T2()) TRes;
   SX_CHECK (x.getSize () > 0);
   ssize_t n = x.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(TRes) * n);

   TRes* resPtr = (TRes*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x.elements[i] - y;
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x.elements[i] - y;
   return typename SxVecResult<TRes,Layout>::VecType(mem, x, x.auxData, x.auxData);
}

// vector - scalar in place (vector contiguous, same type)
template<class T, StrideType Layout, class T2>
void
SxVecCompute<Simd>::subtractInPlace (SxVecRef<T, Layout> &x, const T2 &y)
{
   SX_CHECK (x.getSize () > 0);
   ssize_t n = x.getSize ();

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         x.elements[i] -=  y;
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         x.elements[i] -=  y;
   SX_VALIDATE_VECTOR(x);
}

// vector * scalar (vector contiguous, same type)
template<class T, StrideType Layout, class T2>
typename SxVecResult<decltype(T(0) + T2(0)), Layout>::VecType
SxVecCompute<Simd>::multiply (const SxVecRef<T, Layout> &x, const T2 &y)
{
   typedef decltype(T() + T2()) TRes;
   SX_CHECK (x.getSize () > 0);
   ssize_t n = x.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(TRes) * n);

   TRes* resPtr = (TRes*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x.elements[i] * y;
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x.elements[i] * y;
   return typename SxVecResult<TRes,Layout>::VecType(mem, x, x.auxData, x.auxData);
}

// vector * scalar in place (vector contiguous, same type)
template<class T, StrideType Layout, class T2>
void
SxVecCompute<Simd>::multiplyInPlace (SxVecRef<T, Layout> &x, const T2 &y)
{
   SX_CHECK (x.getSize () > 0);
   ssize_t n = x.getSize ();

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         x.elements[i] *=  y;
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         x.elements[i] *=  y;
   SX_VALIDATE_VECTOR(x);
}

// vector / scalar (vector contiguous, same type)
template<class T, StrideType Layout, class T2>
typename SxVecResult<decltype(T(0) + T2(0)), Layout>::VecType
SxVecCompute<Simd>::Divide<true>::divide (const SxVecRef<T, Layout> &x, const T2 &y)
{
   typedef decltype(T() + T2()) TRes;
   SX_CHECK (x.getSize () > 0);
   ssize_t n = x.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(TRes) * n);

   TRes* resPtr = (TRes*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x.elements[i] / y;
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x.elements[i] / y;
   return typename SxVecResult<TRes,Layout>::VecType(mem, x, x.auxData, x.auxData);
}

// vector / scalar in place (vector contiguous, same type)
template<class T, StrideType Layout, class T2>
void
SxVecCompute<Simd>::Divide<true>::divideInPlace (SxVecRef<T, Layout> &x, const T2 &y)
{
   SX_CHECK (x.getSize () > 0);
   ssize_t n = x.getSize ();

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         x.elements[i] /=  y;
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         x.elements[i] /=  y;
   SX_VALIDATE_VECTOR(x);
}
// --- Simd_OP_scalar

// vector * scalar (vector contiguous)
template<class T, class T2>
void SxVecCompute<Simd>::multiply (decltype(T(0) + T2(0)) *resPtr, ssize_t n,
                                   const T *x, const T2 &y)
{
   SX_CHECK (n > 0, n);
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x[i] * y;
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = x[i] * y;
}


// --- The next 48 lines were generated from snippets/SxVecCompute.hpp snippet Simd_scalar_vec
// scalar - vector (vector contiguous)
template<class T, StrideType Layout, class T2>
typename SxVecResult<decltype(T(0) + T2(0)), Layout>::VecType
SxVecCompute<Simd>::subtract (const T2 &y, const SxVecRef<T, Layout> &x)
{
   typedef decltype(T() + T2()) TRes;
   SX_CHECK (x.getSize () > 0);
   ssize_t n = x.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(TRes) * n);

   TRes* resPtr = (TRes*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = y - x.elements[i];
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = y - x.elements[i];
   return typename SxVecResult<TRes,Layout>::VecType(mem, x, x.auxData, x.auxData);
}

// scalar / vector (vector contiguous)
template<class T, StrideType Layout, class T2>
typename SxVecResult<decltype(T(0) + T2(0)), Layout>::VecType
SxVecCompute<Simd>::divide (const T2 &y, const SxVecRef<T, Layout> &x)
{
   typedef decltype(T() + T2()) TRes;
   SX_CHECK (x.getSize () > 0);
   ssize_t n = x.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(TRes) * n);

   TRes* resPtr = (TRes*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = y / x.elements[i];
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i)
         resPtr[i] = y / x.elements[i];
   return typename SxVecResult<TRes,Layout>::VecType(mem, x, x.auxData, x.auxData);
}
// --- Simd_scalar_vec

// --- The next 432 lines were generated from snippets/SxVecCompute.hpp snippet Simd_vecfunc
// Unary minus
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::minus (const SxVecRef<T, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

   T* resPtr = (T*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         resPtr[i] = -x ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         resPtr[i] = -x ;
      }
   return typename SxVecResult<T,Layout>::VecType(mem, vec, vec.auxData, vec.auxData);
}

// Unary minus
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::minusInPlace (SxVecRef<T, Layout> &&vec)
{
   SX_CHECK (vec.allocMem);
   if (!vec.isOnlyRef ()) return minus(vec);
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         vec.elements[i] = -x ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         vec.elements[i] = -x ;
      }
   SX_VALIDATE_VECTOR(vec);
   return std::move(vec);
}

// Square
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::sqr (const SxVecRef<T, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

   T* resPtr = (T*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         resPtr[i] = x * x ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         resPtr[i] = x * x ;
      }
   return typename SxVecResult<T,Layout>::VecType(mem, vec, vec.auxData, vec.auxData);
}

// Square
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::sqrInPlace (SxVecRef<T, Layout> &&vec)
{
   SX_CHECK (vec.allocMem);
   if (!vec.isOnlyRef ()) return sqr(vec);
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         vec.elements[i] = x * x ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         vec.elements[i] = x * x ;
      }
   SX_VALIDATE_VECTOR(vec);
   return std::move(vec);
}

// Cube (x*x*x)
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::cub (const SxVecRef<T, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

   T* resPtr = (T*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         resPtr[i] = x*x*x ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         resPtr[i] = x*x*x ;
      }
   return typename SxVecResult<T,Layout>::VecType(mem, vec, vec.auxData, vec.auxData);
}

// Cube (x*x*x)
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::cubInPlace (SxVecRef<T, Layout> &&vec)
{
   SX_CHECK (vec.allocMem);
   if (!vec.isOnlyRef ()) return cub(vec);
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         vec.elements[i] = x*x*x ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         vec.elements[i] = x*x*x ;
      }
   SX_VALIDATE_VECTOR(vec);
   return std::move(vec);
}

// Exponential
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::exp (const SxVecRef<T, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

   T* resPtr = (T*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         resPtr[i] = ::exp(x) ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         resPtr[i] = ::exp(x) ;
      }
   return typename SxVecResult<T,Layout>::VecType(mem, vec, vec.auxData, vec.auxData);
}

// Exponential
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::expInPlace (SxVecRef<T, Layout> &&vec)
{
   SX_CHECK (vec.allocMem);
   if (!vec.isOnlyRef ()) return exp(vec);
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         vec.elements[i] = ::exp(x) ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         vec.elements[i] = ::exp(x) ;
      }
   SX_VALIDATE_VECTOR(vec);
   return std::move(vec);
}

// Square root
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::sqrt (const SxVecRef<T, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

   T* resPtr = (T*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         resPtr[i] = ::sqrt(x) ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         resPtr[i] = ::sqrt(x) ;
      }
   return typename SxVecResult<T,Layout>::VecType(mem, vec, vec.auxData, vec.auxData);
}

// Square root
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::sqrtInPlace (SxVecRef<T, Layout> &&vec)
{
   SX_CHECK (vec.allocMem);
   if (!vec.isOnlyRef ()) return sqrt(vec);
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         vec.elements[i] = ::sqrt(x) ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         vec.elements[i] = ::sqrt(x) ;
      }
   SX_VALIDATE_VECTOR(vec);
   return std::move(vec);
}

// Logarithm
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::log (const SxVecRef<T, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

   T* resPtr = (T*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         resPtr[i] = ::log(x) ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         resPtr[i] = ::log(x) ;
      }
   return typename SxVecResult<T,Layout>::VecType(mem, vec, vec.auxData, vec.auxData);
}

// Logarithm
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::logInPlace (SxVecRef<T, Layout> &&vec)
{
   SX_CHECK (vec.allocMem);
   if (!vec.isOnlyRef ()) return log(vec);
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         vec.elements[i] = ::log(x) ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         vec.elements[i] = ::log(x) ;
      }
   SX_VALIDATE_VECTOR(vec);
   return std::move(vec);
}

// Error function
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::erf (const SxVecRef<T, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

   T* resPtr = (T*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         resPtr[i] = ::derf(x) ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         resPtr[i] = ::derf(x) ;
      }
   return typename SxVecResult<T,Layout>::VecType(mem, vec, vec.auxData, vec.auxData);
}

// Error function
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::erfInPlace (SxVecRef<T, Layout> &&vec)
{
   SX_CHECK (vec.allocMem);
   if (!vec.isOnlyRef ()) return erf(vec);
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         vec.elements[i] = ::derf(x) ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         vec.elements[i] = ::derf(x) ;
      }
   SX_VALIDATE_VECTOR(vec);
   return std::move(vec);
}

// Complementary error function
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::erfc (const SxVecRef<T, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

   T* resPtr = (T*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         resPtr[i] = ::derfc(x) ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         resPtr[i] = ::derfc(x) ;
      }
   return typename SxVecResult<T,Layout>::VecType(mem, vec, vec.auxData, vec.auxData);
}

// Complementary error function
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::erfcInPlace (SxVecRef<T, Layout> &&vec)
{
   SX_CHECK (vec.allocMem);
   if (!vec.isOnlyRef ()) return erfc(vec);
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         vec.elements[i] = ::derfc(x) ;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         T x = vec.elements[i];
         vec.elements[i] = ::derfc(x) ;
      }
   SX_VALIDATE_VECTOR(vec);
   return std::move(vec);
}
// --- Simd_vecfunc
// --- The next 52 lines were generated from snippets/SxVecCompute.hpp snippet Simd_abs_complex
template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::abs (const SxVecRef<SxComplex<T>, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

   T* resPtr = (T*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         const SxComplex<T> &x = vec.elements[i];
         resPtr[i] = sqrt(x.absSqr ());
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         const SxComplex<T> &x = vec.elements[i];
         resPtr[i] = sqrt(x.absSqr ());
      }
   return typename SxVecResult<T,Layout>::VecType(mem, vec, vec.auxData, vec.auxData);
}

template<class T, StrideType Layout>
typename SxVecResult<T, Layout>::VecType
SxVecCompute<Simd>::absSqr (const SxVecRef<SxComplex<T>, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

   T* resPtr = (T*)mem->mem;
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         const SxComplex<T> &x = vec.elements[i];
         resPtr[i] = x.absSqr ();
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         const SxComplex<T> &x = vec.elements[i];
         resPtr[i] = x.absSqr ();
      }
   return typename SxVecResult<T,Layout>::VecType(mem, vec, vec.auxData, vec.auxData);
}
// --- Simd_abs_complex

// --- The next 252 lines were generated from snippets/SxVecCompute.hpp snippet UseIterator_OP_scalar
// vector x + scalar y
template<class T, StrideType Layout, class T2>
SxVector<decltype(T(0) + T2(0))>
SxVecCompute<UseIterator>::add (const SxVecRef<T, Layout>  &x,
                                const T2 &y)
{
   typedef decltype(T() + T2()) TRes;
   ssize_t n = x.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(TRes) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> xIt = x.begin ();
      xIt += iStart;
      TRes* resPtr = (TRes*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt)
         resPtr[i] = *xIt + y;
   } else
#endif
   {
      SxVecConstIt<T,Layout>   xIt = x.begin ();
      TRes* resPtr = (TRes*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++xIt)
         resPtr[i] = *xIt + y;
   }
   return SxVector<TRes>(mem, x, x.auxData, x.auxData);
}

// vector x + scalar y (in place)
template<class T, StrideType Layout, class T2>
void
SxVecCompute<UseIterator>::addInPlace (SxVecRef<T, Layout>  &x,
                                           const T2 &y)
{
   ssize_t n = x.getSize ();

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecIt<T,Layout> xIt = x.begin ();
      xIt += iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt) {
         *xIt += y;
         SX_CHECK_NUM (*xIt);
      }
   } else
#endif
   {
      SxVecIt<T,Layout>   xIt = x.begin ();
      for (ssize_t i = 0 ; i < n; ++i, ++xIt) {
         *xIt += y;
         SX_CHECK_NUM (*xIt);
      }
   }
}

// vector x - scalar y
template<class T, StrideType Layout, class T2>
SxVector<decltype(T(0) + T2(0))>
SxVecCompute<UseIterator>::subtract (const SxVecRef<T, Layout>  &x,
                                const T2 &y)
{
   typedef decltype(T() + T2()) TRes;
   ssize_t n = x.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(TRes) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> xIt = x.begin ();
      xIt += iStart;
      TRes* resPtr = (TRes*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt)
         resPtr[i] = *xIt - y;
   } else
#endif
   {
      SxVecConstIt<T,Layout>   xIt = x.begin ();
      TRes* resPtr = (TRes*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++xIt)
         resPtr[i] = *xIt - y;
   }
   return SxVector<TRes>(mem, x, x.auxData, x.auxData);
}

// vector x - scalar y (in place)
template<class T, StrideType Layout, class T2>
void
SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<T, Layout>  &x,
                                           const T2 &y)
{
   ssize_t n = x.getSize ();

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecIt<T,Layout> xIt = x.begin ();
      xIt += iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt) {
         *xIt -= y;
         SX_CHECK_NUM (*xIt);
      }
   } else
#endif
   {
      SxVecIt<T,Layout>   xIt = x.begin ();
      for (ssize_t i = 0 ; i < n; ++i, ++xIt) {
         *xIt -= y;
         SX_CHECK_NUM (*xIt);
      }
   }
}

// vector x * scalar y
template<class T, StrideType Layout, class T2>
SxVector<decltype(T(0) + T2(0))>
SxVecCompute<UseIterator>::multiply (const SxVecRef<T, Layout>  &x,
                                const T2 &y)
{
   typedef decltype(T() + T2()) TRes;
   ssize_t n = x.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(TRes) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> xIt = x.begin ();
      xIt += iStart;
      TRes* resPtr = (TRes*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt)
         resPtr[i] = *xIt * y;
   } else
#endif
   {
      SxVecConstIt<T,Layout>   xIt = x.begin ();
      TRes* resPtr = (TRes*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++xIt)
         resPtr[i] = *xIt * y;
   }
   return SxVector<TRes>(mem, x, x.auxData, x.auxData);
}

// vector x * scalar y (in place)
template<class T, StrideType Layout, class T2>
void
SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<T, Layout>  &x,
                                           const T2 &y)
{
   ssize_t n = x.getSize ();

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecIt<T,Layout> xIt = x.begin ();
      xIt += iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt) {
         *xIt *= y;
         SX_CHECK_NUM (*xIt);
      }
   } else
#endif
   {
      SxVecIt<T,Layout>   xIt = x.begin ();
      for (ssize_t i = 0 ; i < n; ++i, ++xIt) {
         *xIt *= y;
         SX_CHECK_NUM (*xIt);
      }
   }
}

// vector x / scalar y
template<class T, StrideType Layout, class T2>
SxVector<decltype(T(0) + T2(0))>
SxVecCompute<UseIterator>::Divide<true>::divide (const SxVecRef<T, Layout>  &x,
                                const T2 &y)
{
   typedef decltype(T() + T2()) TRes;
   ssize_t n = x.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(TRes) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> xIt = x.begin ();
      xIt += iStart;
      TRes* resPtr = (TRes*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt)
         resPtr[i] = *xIt / y;
   } else
#endif
   {
      SxVecConstIt<T,Layout>   xIt = x.begin ();
      TRes* resPtr = (TRes*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++xIt)
         resPtr[i] = *xIt / y;
   }
   return SxVector<TRes>(mem, x, x.auxData, x.auxData);
}

// vector x / scalar y (in place)
template<class T, StrideType Layout, class T2>
void
SxVecCompute<UseIterator>::Divide<true>::divideInPlace (SxVecRef<T, Layout>  &x,
                                           const T2 &y)
{
   ssize_t n = x.getSize ();

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecIt<T,Layout> xIt = x.begin ();
      xIt += iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt) {
         *xIt /= y;
         SX_CHECK_NUM (*xIt);
      }
   } else
#endif
   {
      SxVecIt<T,Layout>   xIt = x.begin ();
      for (ssize_t i = 0 ; i < n; ++i, ++xIt) {
         *xIt /= y;
         SX_CHECK_NUM (*xIt);
      }
   }
}
// --- UseIterator_OP_scalar

// --- The next 64 lines were generated from snippets/SxVecCompute.hpp snippet UseIterator_scalar_vec
// scalar y - vector x
template<class T, StrideType Layout, class T2>
SxVector<decltype(T(0) + T2(0))>
SxVecCompute<UseIterator>::subtract (const T2 &y,
                                const SxVecRef<T, Layout>  &x)
{
   typedef decltype(T() + T2()) TRes;
   ssize_t n = x.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(TRes) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout>   xIt = x.begin ();
      xIt += iStart;
      TRes* resPtr = (TRes*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt)
         resPtr[i] = y - *xIt;
   } else
#endif
   {
      SxVecConstIt<T,Layout>   xIt = x.begin ();
      TRes* resPtr = (TRes*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++xIt)
         resPtr[i] = y - *xIt;
   }
   return SxVector<TRes>(mem, x, x.auxData, x.auxData);
}

// scalar y / vector x
template<class T, StrideType Layout, class T2>
SxVector<decltype(T(0) + T2(0))>
SxVecCompute<UseIterator>::divide (const T2 &y,
                                const SxVecRef<T, Layout>  &x)
{
   typedef decltype(T() + T2()) TRes;
   ssize_t n = x.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(TRes) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout>   xIt = x.begin ();
      xIt += iStart;
      TRes* resPtr = (TRes*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++xIt)
         resPtr[i] = y / *xIt;
   } else
#endif
   {
      SxVecConstIt<T,Layout>   xIt = x.begin ();
      TRes* resPtr = (TRes*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++xIt)
         resPtr[i] = y / *xIt;
   }
   return SxVector<TRes>(mem, x, x.auxData, x.auxData);
}
// --- UseIterator_scalar_vec

// --- The next 272 lines were generated from snippets/SxVecCompute.hpp snippet UseIterator_vecfunc
// Unary minus
template<class T, StrideType Layout>
SxVector<T> SxVecCompute<UseIterator>::minus (const SxVecRef<T, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> it = vec.begin ();
      it += iStart;
      T* resPtr = (T*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         T x = *it;
         resPtr[i] = -x ;
      }
   } else
#endif
   {
      SxVecConstIt<T,Layout> it = vec.begin ();
      T* resPtr = (T*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         T x = *it;
         resPtr[i] = -x ;
      }
   }
   return SxVector<T>(mem, vec, vec.auxData, vec.auxData);
}

// Square
template<class T, StrideType Layout>
SxVector<T> SxVecCompute<UseIterator>::sqr (const SxVecRef<T, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> it = vec.begin ();
      it += iStart;
      T* resPtr = (T*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         T x = *it;
         resPtr[i] = x * x ;
      }
   } else
#endif
   {
      SxVecConstIt<T,Layout> it = vec.begin ();
      T* resPtr = (T*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         T x = *it;
         resPtr[i] = x * x ;
      }
   }
   return SxVector<T>(mem, vec, vec.auxData, vec.auxData);
}

// Cube (x*x*x)
template<class T, StrideType Layout>
SxVector<T> SxVecCompute<UseIterator>::cub (const SxVecRef<T, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> it = vec.begin ();
      it += iStart;
      T* resPtr = (T*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         T x = *it;
         resPtr[i] = x*x*x ;
      }
   } else
#endif
   {
      SxVecConstIt<T,Layout> it = vec.begin ();
      T* resPtr = (T*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         T x = *it;
         resPtr[i] = x*x*x ;
      }
   }
   return SxVector<T>(mem, vec, vec.auxData, vec.auxData);
}

// Exponential
template<class T, StrideType Layout>
SxVector<T> SxVecCompute<UseIterator>::exp (const SxVecRef<T, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> it = vec.begin ();
      it += iStart;
      T* resPtr = (T*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         T x = *it;
         resPtr[i] = ::exp(x) ;
      }
   } else
#endif
   {
      SxVecConstIt<T,Layout> it = vec.begin ();
      T* resPtr = (T*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         T x = *it;
         resPtr[i] = ::exp(x) ;
      }
   }
   return SxVector<T>(mem, vec, vec.auxData, vec.auxData);
}

// Square root
template<class T, StrideType Layout>
SxVector<T> SxVecCompute<UseIterator>::sqrt (const SxVecRef<T, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> it = vec.begin ();
      it += iStart;
      T* resPtr = (T*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         T x = *it;
         resPtr[i] = ::sqrt(x) ;
      }
   } else
#endif
   {
      SxVecConstIt<T,Layout> it = vec.begin ();
      T* resPtr = (T*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         T x = *it;
         resPtr[i] = ::sqrt(x) ;
      }
   }
   return SxVector<T>(mem, vec, vec.auxData, vec.auxData);
}

// Logarithm
template<class T, StrideType Layout>
SxVector<T> SxVecCompute<UseIterator>::log (const SxVecRef<T, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> it = vec.begin ();
      it += iStart;
      T* resPtr = (T*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         T x = *it;
         resPtr[i] = ::log(x) ;
      }
   } else
#endif
   {
      SxVecConstIt<T,Layout> it = vec.begin ();
      T* resPtr = (T*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         T x = *it;
         resPtr[i] = ::log(x) ;
      }
   }
   return SxVector<T>(mem, vec, vec.auxData, vec.auxData);
}

// Error function
template<class T, StrideType Layout>
SxVector<T> SxVecCompute<UseIterator>::erf (const SxVecRef<T, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> it = vec.begin ();
      it += iStart;
      T* resPtr = (T*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         T x = *it;
         resPtr[i] = ::derf(x) ;
      }
   } else
#endif
   {
      SxVecConstIt<T,Layout> it = vec.begin ();
      T* resPtr = (T*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         T x = *it;
         resPtr[i] = ::derf(x) ;
      }
   }
   return SxVector<T>(mem, vec, vec.auxData, vec.auxData);
}

// Complementary error function
template<class T, StrideType Layout>
SxVector<T> SxVecCompute<UseIterator>::erfc (const SxVecRef<T, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<T,Layout> it = vec.begin ();
      it += iStart;
      T* resPtr = (T*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         T x = *it;
         resPtr[i] = ::derfc(x) ;
      }
   } else
#endif
   {
      SxVecConstIt<T,Layout> it = vec.begin ();
      T* resPtr = (T*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         T x = *it;
         resPtr[i] = ::derfc(x) ;
      }
   }
   return SxVector<T>(mem, vec, vec.auxData, vec.auxData);
}
// --- UseIterator_vecfunc

template<StrideType Layout>
SxVector<double> SxVecCompute<UseIterator>::pow (const SxVecRef<double, Layout> &vec, double a)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(double) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<double,Layout> it = vec.begin ();
      it += iStart;
      double* resPtr = (double*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         resPtr[i] = ::pow(*it, a);
      }
   } else
#endif
   {
      SxVecConstIt<double,Layout> it = vec.begin ();
      double* resPtr = (double*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         resPtr[i] = ::pow(*it, a);
      }
   }
   return SxVector<double>(mem, vec, vec.auxData, vec.auxData);
}


// --- The next 246 lines were generated from snippets/SxVecCompute.hpp snippet UseIterator_abs_real
// Absolute value
template<StrideType Layout>
SxVector<int> SxVecCompute<UseIterator>::abs (const SxVecRef<int,Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(int) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<int,Layout> it = vec.begin ();
      it += iStart;
      int* resPtr = (int*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         int x = *it;
         resPtr[i] = ::abs(x) ;
      }
   } else
#endif
   {
      SxVecConstIt<int,Layout> it = vec.begin ();
      int* resPtr = (int*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         int x = *it;
         resPtr[i] = ::abs(x) ;
      }
   }
   return SxVector<int>(mem, vec, vec.auxData, vec.auxData);
}

template<StrideType Layout>
inline SxVector<int>
SxVecCompute<UseIterator>::absInPlace (SxVecRef<int,Layout> &&x)
{
   return abs(x);
}

// Absolute value squared
template<StrideType Layout>
SxVector<int> SxVecCompute<UseIterator>::absSqr (const SxVecRef<int,Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(int) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<int,Layout> it = vec.begin ();
      it += iStart;
      int* resPtr = (int*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         int x = *it;
         resPtr[i] = x*x ;
      }
   } else
#endif
   {
      SxVecConstIt<int,Layout> it = vec.begin ();
      int* resPtr = (int*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         int x = *it;
         resPtr[i] = x*x ;
      }
   }
   return SxVector<int>(mem, vec, vec.auxData, vec.auxData);
}

template<StrideType Layout>
inline SxVector<int>
SxVecCompute<UseIterator>::absSqrInPlace (SxVecRef<int,Layout> &&x)
{
   return absSqr(x);
}

// Absolute value
template<StrideType Layout>
SxVector<float> SxVecCompute<UseIterator>::abs (const SxVecRef<float,Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(float) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<float,Layout> it = vec.begin ();
      it += iStart;
      float* resPtr = (float*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         float x = *it;
         resPtr[i] = ::fabsf(x) ;
      }
   } else
#endif
   {
      SxVecConstIt<float,Layout> it = vec.begin ();
      float* resPtr = (float*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         float x = *it;
         resPtr[i] = ::fabsf(x) ;
      }
   }
   return SxVector<float>(mem, vec, vec.auxData, vec.auxData);
}

template<StrideType Layout>
inline SxVector<float>
SxVecCompute<UseIterator>::absInPlace (SxVecRef<float,Layout> &&x)
{
   return abs(x);
}

// Absolute value squared
template<StrideType Layout>
SxVector<float> SxVecCompute<UseIterator>::absSqr (const SxVecRef<float,Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(float) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<float,Layout> it = vec.begin ();
      it += iStart;
      float* resPtr = (float*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         float x = *it;
         resPtr[i] = x*x ;
      }
   } else
#endif
   {
      SxVecConstIt<float,Layout> it = vec.begin ();
      float* resPtr = (float*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         float x = *it;
         resPtr[i] = x*x ;
      }
   }
   return SxVector<float>(mem, vec, vec.auxData, vec.auxData);
}

template<StrideType Layout>
inline SxVector<float>
SxVecCompute<UseIterator>::absSqrInPlace (SxVecRef<float,Layout> &&x)
{
   return absSqr(x);
}

// Absolute value
template<StrideType Layout>
SxVector<double> SxVecCompute<UseIterator>::abs (const SxVecRef<double,Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(double) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<double,Layout> it = vec.begin ();
      it += iStart;
      double* resPtr = (double*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         double x = *it;
         resPtr[i] = ::fabs(x) ;
      }
   } else
#endif
   {
      SxVecConstIt<double,Layout> it = vec.begin ();
      double* resPtr = (double*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         double x = *it;
         resPtr[i] = ::fabs(x) ;
      }
   }
   return SxVector<double>(mem, vec, vec.auxData, vec.auxData);
}

template<StrideType Layout>
inline SxVector<double>
SxVecCompute<UseIterator>::absInPlace (SxVecRef<double,Layout> &&x)
{
   return abs(x);
}

// Absolute value squared
template<StrideType Layout>
SxVector<double> SxVecCompute<UseIterator>::absSqr (const SxVecRef<double,Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(double) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<double,Layout> it = vec.begin ();
      it += iStart;
      double* resPtr = (double*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         double x = *it;
         resPtr[i] = x*x ;
      }
   } else
#endif
   {
      SxVecConstIt<double,Layout> it = vec.begin ();
      double* resPtr = (double*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         double x = *it;
         resPtr[i] = x*x ;
      }
   }
   return SxVector<double>(mem, vec, vec.auxData, vec.auxData);
}

template<StrideType Layout>
inline SxVector<double>
SxVecCompute<UseIterator>::absSqrInPlace (SxVecRef<double,Layout> &&x)
{
   return absSqr(x);
}
// --- UseIterator_abs_real

// --- The next 68 lines were generated from snippets/SxVecCompute.hpp snippet UseIterator_abs_complex
template<class T, StrideType Layout>
SxVector<T>
SxVecCompute<UseIterator>::abs (const SxVecRef<SxComplex<T>, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<SxComplex<T>,Layout> it = vec.begin ();
      it += iStart;
      T* resPtr = (T*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         const SxComplex<T> &x = *it;
         resPtr[i] = sqrt(x.absSqr ()) ;
      }
   } else
#endif
   {
      SxVecConstIt<SxComplex<T>,Layout> it = vec.begin ();
      T* resPtr = (T*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         const SxComplex<T> &x = *it;
         resPtr[i] = sqrt(x.absSqr ()) ;
      }
   }
   return SxVector<T>(mem, vec, vec.auxData, vec.auxData);
}

template<class T, StrideType Layout>
SxVector<T>
SxVecCompute<UseIterator>::absSqr (const SxVecRef<SxComplex<T>, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<SxComplex<T>,Layout> it = vec.begin ();
      it += iStart;
      T* resPtr = (T*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         const SxComplex<T> &x = *it;
         resPtr[i] = x.absSqr () ;
      }
   } else
#endif
   {
      SxVecConstIt<SxComplex<T>,Layout> it = vec.begin ();
      T* resPtr = (T*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         const SxComplex<T> &x = *it;
         resPtr[i] = x.absSqr () ;
      }
   }
   return SxVector<T>(mem, vec, vec.auxData, vec.auxData);
}
// --- UseIterator_abs_complex

template<class T>
inline T SxVecCompute<Simd>::sum (const T* elements, ssize_t n)
{
   SX_CHECK (n>0, n);
   T res(0);
#ifdef USE_OPENMP
# pragma omp parallel for reduction(+:res) if (n > sxChunkSize)
#endif
   for (ssize_t i = 0; i < n; i++)
      res += elements[i];
   SX_CHECK_NUM (res);
   return res;
}

template<class T>
T SxVecCompute<Simd>::product (const T* elements, ssize_t n)
{
   SX_CHECK (n>0, n);
   T res(1);
#ifdef USE_OPENMP
# pragma omp parallel for reduction(*:res) if (n > sxChunkSize)
#endif
   for (ssize_t i = 0; i < n; i++)
      res *= elements[i];
   SX_CHECK_NUM (res);
   return res;
}

template<class T>
T SxVecCompute<Simd>::normSqr (const T* elements, ssize_t n)
{
   SX_CHECK (n>0, n);
   T res(1);
#ifdef USE_OPENMP
# pragma omp parallel for reduction(+:res) if (n > sxChunkSize)
#endif
   for (ssize_t i = 0; i < n; i++)
      res *= elements[i] * elements[i];
   SX_CHECK_NUM (res);
   return res;
}

// --- The next 36 lines were generated from snippets/SxVecCompute.hpp snippet SimdCol_reduce
/// Add all elements
template<class T>
T SxVecCompute<SimdCol>::sum (const SxVecRef<T,SubMatrix> &x)
{
   T res = 0.;
   ssize_t cStride = x.getColStride ();
   for (ssize_t ic = 0; ic < x.getNCols (); ic++)
      res += SxVecCompute<Simd>::sum (x.elements + ic * cStride, x.getNRows ());
   SX_CHECK_NUM(res);
   return res;
}

/// Multiply all elements
template<class T>
T SxVecCompute<SimdCol>::product (const SxVecRef<T,SubMatrix> &x)
{
   T res = 1.;
   ssize_t cStride = x.getColStride ();
   for (ssize_t ic = 0; ic < x.getNCols (); ic++)
      res *= SxVecCompute<Simd>::product (x.elements + ic * cStride, x.getNRows ());
   SX_CHECK_NUM(res);
   return res;
}

/// Add abs. squares of all elements
template<class T>
typename SxTypeMapper<T>::TReal SxVecCompute<SimdCol>::normSqr (const SxVecRef<T,SubMatrix> &x)
{
   typename SxTypeMapper<T>::TReal res = 0.;
   ssize_t cStride = x.getColStride ();
#  ifdef USE_OPENMP
#  pragma omp parallel for \
     if (x.getSize () > sxChunkSize && x.getNCols () >= omp_get_max_threads ())
#  endif
   for (ssize_t ic = 0; ic < x.getNCols (); ic++)
      res += SxVecCompute<Simd>::normSqr (x.elements + ic * cStride, x.getNRows ());
   SX_CHECK_NUM(res);
   return res;
}
// --- SimdCol_reduce

// --- The next 33 lines were generated from snippets/SxVecCompute.hpp snippet UseIterator_reduce
/// Add all elements
template<class T, StrideType Layout>
T SxVecCompute<UseIterator>::sum (const SxVecRef<T,Layout> &x)
{
   T res = 0.;
   for (SxVecConstIt<T,Layout> it = x.begin (); it != x.end (); ++it)
      res += (*it);
   SX_CHECK_NUM(res);
   return res;
}

/// Multiply all elements
template<class T, StrideType Layout>
T SxVecCompute<UseIterator>::product (const SxVecRef<T,Layout> &x)
{
   T res = 1.;
   for (SxVecConstIt<T,Layout> it = x.begin (); it != x.end (); ++it)
      res *= (*it);
   SX_CHECK_NUM(res);
   return res;
}

/// Add abs. squares of all elements
template<class T, StrideType Layout>
typename SxTypeMapper<T>::TReal SxVecCompute<UseIterator>::normSqr (const SxVecRef<T,Layout> &x)
{
   typename SxTypeMapper<T>::TReal res = 0.;
   for (SxVecConstIt<T,Layout> it = x.begin (); it != x.end (); ++it)
      res += absSqr(*it);
   SX_CHECK_NUM(res);
   return res;
}
// --- UseIterator_reduce

template<class T, StrideType Layout>
SxVecRef<T, Layout>
SxVecCompute<Simd>::real (const SxVecRef<SxComplex<T>, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

   T* resPtr = (T*)mem->mem;
   // NOTE: it might be advantageous to do explicit AVX/SSE implementation
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         resPtr[i] = vec.elements[i].re;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         resPtr[i] = vec.elements[i].re;
      }
   return SxVecRef<T,Layout>(mem, vec, vec.auxData, vec.auxData);
}

template<class T, StrideType Layout>
SxVecRef<T, Layout>
SxVecCompute<Simd>::imag (const SxVecRef<SxComplex<T>, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * n);

   T* resPtr = (T*)mem->mem;
   // NOTE: it might be advantageous to do explicit AVX/SSE implementation
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         resPtr[i] = vec.elements[i].im;
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         resPtr[i] = vec.elements[i].im;
      }
   return SxVecRef<T,Layout>(mem, vec, vec.auxData, vec.auxData);
}

template<class T, StrideType Layout>
typename SxVecResult<SxComplex<T>, Layout>::VecType
SxVecCompute<Simd>::conj (const SxVecRef<SxComplex<T>, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(SxComplex<T>) * n);

   SxComplex<T>* resPtr = (SxComplex<T>*)mem->mem;
   // NOTE: it might be advantageous to do explicit AVX/SSE implementation
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         resPtr[i] = vec.elements[i].conj ();
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         resPtr[i] = vec.elements[i].conj ();
      }
   return typename SxVecResult<SxComplex<T>,Layout>::VecType(mem, vec, vec.auxData, vec.auxData);
}

template<class T, StrideType Layout>
typename SxVecResult<SxComplex<T>, Layout>::VecType
SxVecCompute<Simd>::conjInPlace (SxVecRef<SxComplex<T>, Layout> &&vec)
{
   SX_CHECK (vec.allocMem);
   if (!vec.isOnlyRef ()) return conj((const SxVecRef<SxComplex<T>,Layout>&)vec);
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();

   // NOTE: it might be advantageous to do explicit AVX/SSE implementation
#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel for simd
      for (ssize_t i = 0; i < n; ++i) {
         vec.elements[i] = vec.elements[i].conj ();
      }
   else
#     pragma omp simd
#endif
      for (ssize_t i = 0; i < n; ++i) {
         vec.elements[i] = vec.elements[i].conj ();
      }
   return std::move(vec);
}

template<class T, StrideType Layout>
SxVector<SxComplex<T> >
SxVecCompute<UseIterator>::conj (const SxVecRef<SxComplex<T>, Layout> &vec)
{
   SX_CHECK (vec.getSize () > 0);
   ssize_t n = vec.getSize ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(SxComplex<T>) * n);

#ifdef USE_OPENMP
   if (n > sxChunkSize)
#     pragma omp parallel
   {
      ssize_t nChunk = sx_omp_chunk (n);
      ssize_t iStart = sx_omp_start (nChunk, n);
      SxVecConstIt<SxComplex<T>,Layout> it = vec.begin ();
      it += iStart;
      SxComplex<T>* resPtr = (SxComplex<T>*)mem->mem + iStart;
      for (ssize_t i = 0 ; i < nChunk ; ++i, ++it) {
         resPtr[i] = (*it).conj ();
      }
   } else
#endif
   {
      SxVecConstIt<SxComplex<T>,Layout> it = vec.begin ();
      SxComplex<T>* resPtr = (SxComplex<T>*)mem->mem;
      for (ssize_t i = 0 ; i < n; ++i, ++it) {
         resPtr[i] = (*it).conj ();
      }
   }
   return SxVector<SxComplex<T> >(mem, vec, vec.auxData, vec.auxData);
}

// --- The next 168 lines were generated from snippets/SxVecCompute.hpp snippet Transpose
// transpose N x M (on src) to M x N (on dest)
// strideS and strideD are the column strides for src/dest respectively
template<class T>
void SxVecCompute<Simd>::transpose (ssize_t N, ssize_t M,
                                    const T* src, ssize_t strideS,
                                    T* dest, ssize_t strideD)
{
   // split into blocks of max 2048 x 2048 (= 4 million elements)
   for (ssize_t jj = 0; jj < M; jj+= 2048)  {
      ssize_t mm = (M-jj) > 2048 ? 2048 : (M-jj);
      for (ssize_t ii = 0; ii < N; ii+=2048)  {
         ssize_t nn = (N-ii) > 2048 ? 2048 : (N-ii);
         // --- use cache-oblivious transpose on these blocks
         SxCacheOblivious task(uint16_t(nn), uint16_t(mm), 256);
         while (task.splitProblem ())  {
            ssize_t offX = ii + task.offX;
            ssize_t offY = jj + task.offY;
            const T* sPtr = src  + offX + offY * strideS;
                  T* dPtr = dest + offY + offX * strideD;
            for (uint16_t j = 0; j < task.nY; j++)
               for (uint16_t i = 0; i < task.nX; i++)
                  dPtr[j + i * strideD] = sPtr[i + j * strideS];
         }
      }
   }
}

// transpose for Compact or SubMatrix
template<class T, StrideType Layout>
SxVector<T > SxVecCompute<Simd>::transpose (const SxVecRef<T, Layout> &vec)
{
   ssize_t n = vec.getSize ();
   SxVector<T > res(n);
   res.reshape (vec.nCols, vec.nRows);
   transpose (vec.nRows, vec.nCols, vec.elements, vec.getColStride (), res.elements, vec.nCols);
   res.auxData = vec.auxData;
   res.auxData.basisPtr = NULL; // basis no longer makes sense
   return res;
}

// transpose N x M (on src) to M x N (on dest)
// rStride and cStride are the row and column strides for src,
// strideD is the column stride for the (compact) result respectively
template<class T>
void SxVecCompute<SimdIt>::transpose (ssize_t N, ssize_t M,
                                      const T* src, ssize_t rStride, ssize_t cStride,
                                      T* dest, ssize_t strideD)
{
   // split into blocks of max 2048 x 2048 (= 4 million elements)
   for (ssize_t jj = 0; jj < M; jj+= 2048)  {
      ssize_t mm = (M-jj) > 2048 ? 2048 : (M-jj);
      for (ssize_t ii = 0; ii < N; ii+=2048)  {
         ssize_t nn = (N-ii) > 2048 ? 2048 : (N-ii);
         // --- use cache-oblivious transpose on these blocks
         SxCacheOblivious task(uint16_t(nn), uint16_t(mm), 256);
         while (task.splitProblem ())  {
            ssize_t offX = ii + task.offX;
            ssize_t offY = jj + task.offY;
            // --- small case
            const T* sPtr = src  + rStride * offX + offY * cStride;
                  T* dPtr = dest + offY + offX * strideD;
            for (uint16_t j = 0; j < task.nY; j++)
               for (uint16_t i = 0; i < task.nX; i++)
                  dPtr[j + i * strideD] = sPtr[rStride * i + j * cStride];
         }
      }
   }
}

// transpose for Strided and GeneralMatrix
template<class T, StrideType Layout>
SxVector<T > SxVecCompute<SimdIt>::transpose (const SxVecRef<T, Layout> &vec)
{
   ssize_t n = vec.getSize ();
   SxVector<T > res(n);
   res.reshape (vec.nCols, vec.nRows);
   transpose (vec.nRows, vec.nCols,
              vec.elements, vec.getRowStride (), vec.getColStride (),
              res.elements, vec.nCols);
   res.auxData = vec.auxData;
   res.auxData.basisPtr = NULL; // basis no longer makes sense
   return res;
}

// adjoint N x M (on src) to M x N (on dest)
// strideS and strideD are the column strides for src/dest respectively
template<class T>
void SxVecCompute<Simd>::adjoint (ssize_t N, ssize_t M,
                                    const SxComplex<T>* src, ssize_t strideS,
                                    SxComplex<T>* dest, ssize_t strideD)
{
   // split into blocks of max 2048 x 2048 (= 4 million elements)
   for (ssize_t jj = 0; jj < M; jj+= 2048)  {
      ssize_t mm = (M-jj) > 2048 ? 2048 : (M-jj);
      for (ssize_t ii = 0; ii < N; ii+=2048)  {
         ssize_t nn = (N-ii) > 2048 ? 2048 : (N-ii);
         // --- use cache-oblivious adjoint on these blocks
         SxCacheOblivious task(uint16_t(nn), uint16_t(mm), 256);
         while (task.splitProblem ())  {
            ssize_t offX = ii + task.offX;
            ssize_t offY = jj + task.offY;
            const SxComplex<T>* sPtr = src  + offX + offY * strideS;
                  SxComplex<T>* dPtr = dest + offY + offX * strideD;
            for (uint16_t j = 0; j < task.nY; j++)
               for (uint16_t i = 0; i < task.nX; i++)
                  dPtr[j + i * strideD] = sPtr[i + j * strideS].conj ();
         }
      }
   }
}

// adjoint for Compact or SubMatrix
template<class T, StrideType Layout>
SxVector<SxComplex<T> > SxVecCompute<Simd>::adjoint (const SxVecRef<SxComplex<T>, Layout> &vec)
{
   ssize_t n = vec.getSize ();
   SxVector<SxComplex<T> > res(n);
   res.reshape (vec.nCols, vec.nRows);
   adjoint (vec.nRows, vec.nCols, vec.elements, vec.getColStride (), res.elements, vec.nCols);
   res.auxData = vec.auxData;
   res.auxData.basisPtr = NULL; // basis no longer makes sense
   return res;
}

// adjoint N x M (on src) to M x N (on dest)
// rStride and cStride are the row and column strides for src,
// strideD is the column stride for the (compact) result respectively
template<class T>
void SxVecCompute<SimdIt>::adjoint (ssize_t N, ssize_t M,
                                      const SxComplex<T>* src, ssize_t rStride, ssize_t cStride,
                                      SxComplex<T>* dest, ssize_t strideD)
{
   // split into blocks of max 2048 x 2048 (= 4 million elements)
   for (ssize_t jj = 0; jj < M; jj+= 2048)  {
      ssize_t mm = (M-jj) > 2048 ? 2048 : (M-jj);
      for (ssize_t ii = 0; ii < N; ii+=2048)  {
         ssize_t nn = (N-ii) > 2048 ? 2048 : (N-ii);
         // --- use cache-oblivious adjoint on these blocks
         SxCacheOblivious task(uint16_t(nn), uint16_t(mm), 256);
         while (task.splitProblem ())  {
            ssize_t offX = ii + task.offX;
            ssize_t offY = jj + task.offY;
            // --- small case
            const SxComplex<T>* sPtr = src  + rStride * offX + offY * cStride;
                  SxComplex<T>* dPtr = dest + offY + offX * strideD;
            for (uint16_t j = 0; j < task.nY; j++)
               for (uint16_t i = 0; i < task.nX; i++)
                  dPtr[j + i * strideD] = sPtr[rStride * i + j * cStride].conj ();
         }
      }
   }
}

// adjoint for Strided and GeneralMatrix
template<class T, StrideType Layout>
SxVector<SxComplex<T> > SxVecCompute<SimdIt>::adjoint (const SxVecRef<SxComplex<T>, Layout> &vec)
{
   ssize_t n = vec.getSize ();
   SxVector<SxComplex<T> > res(n);
   res.reshape (vec.nCols, vec.nRows);
   adjoint (vec.nRows, vec.nCols,
              vec.elements, vec.getRowStride (), vec.getColStride (),
              res.elements, vec.nCols);
   res.auxData = vec.auxData;
   res.auxData.basisPtr = NULL; // basis no longer makes sense
   return res;
}
// --- TransposeInterface

template<class T, StrideType Layout1, StrideType Layout2>
SxVector<T>
SxVecCompute<Simd>::matmult (const SxVecRef<T,Layout1> &a,
                             const SxVecRef<T,Layout2> &b)
{
   SX_CHECK (a.getNRows () > 0);
   SX_CHECK (b.getNCols () > 0);
   SX_CHECK (int(a.getNRows ()) == a.getNRows ());
   SX_CHECK (int(a.getNCols ()) == a.getNCols ());
   SX_CHECK (int(b.getNCols ()) == b.getNCols ());
   SX_CHECK (int(a.getColStride ()) == a.getColStride ());
   SX_CHECK (int(b.getColStride ()) == b.getColStride ());
   ssize_t nRow = a.getNRows (), nCol = b.getNCols ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * nRow * nCol);
   // call BLAS {s,d,c,z}gemm wrapper
   ::matmult ((T*)mem->mem, 0., a.elements, b.elements,
              (int)a.getNRows (), (int)a.getNCols (), (int)b.getNCols (),
              (int)a.getColStride (), (int)b.getColStride (), (int)a.getNRows ());
   return SxVector<T> (mem, SxVecLayout<Compact> (nRow, nCol),
                       a.auxData, a.auxData); /* keep auxData from a ! */
}

template<class T, StrideType Layout1, StrideType Layout2>
SxVector<T> SxVecCompute<SimdCol>::matmult (const SxVecRef<T,Layout1> &a,
                                            const SxVecRef<T,Layout2> &b)
{
   SX_CHECK (a.getNRows () > 0);
   SX_CHECK (b.getNCols () > 0);
   ssize_t nRow = a.getNRows (), nCol = b.getNCols ();
   SX_CHECK (int(nRow) == nRow);
   SX_CHECK (int(a.getColStride ()) == a.getColStride ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * nRow * nCol);
   T *work = new T[1024];
   for (int jj = 0; jj < b.getNCols (); jj += 1024)  {
      ssize_t nj = b.getNCols () - jj;
      if (nj > 1024) nj = 1024;
      for (int ii = 0; ii < b.getNRows (); ii += 1024)  {
         ssize_t ni = b.getNRows () - ii;
         if (ni > 1024) ni = 1024;
         SxCacheOblivious task (uint16_t(ni), uint16_t(nj), 1024);
         while (task.splitProblem ())  {
            ssize_t offi = ii + task.offX;
            ssize_t offj = jj + task.offY;
            int wStride = task.nX;
            ssize_t rStride = b.getRowStride ();
            // collect a small submatrix of b
            for (uint16_t j = 0; j < task.nY; j++)  {
               T *src = b.elements + offi * rStride + (offj + j) * b.getColStride ();
               for (uint16_t i = 0; i < task.nX; i++)
                  work[i + wStride * j] = src[i * rStride];
            }
            // multiply this submatrix (from left) with a into result
            // via BLAS {s,d,c,z}gemm wrapper
            ::matmult ((T*)mem->mem + nRow * offj, (offi == 0) ? 0. : 1.,
                       a.elements + offi * a.getColStride (), work,
                       (int)nRow, task.nX, task.nY,
                       (int)a.getColStride (), wStride, (int)nRow);

         }
      }
   }
   delete [] work;
   return SxVector<T> (mem, SxVecLayout<Compact> (nRow, nCol),
                       a.auxData, a.auxData); /* keep auxData from a ! */
}

template<class T, StrideType Layout1, StrideType Layout2>
SxVector<T> SxVecCompute<SimdIt>::matmult (const SxVecRef<T,Layout1> &a,
                                           const SxVecRef<T,Layout2> &b)
{
   SX_CHECK (a.getNRows () > 0);
   SX_CHECK (b.getNCols () > 0);
   SX_CHECK ((int)a.getNRows () == a.getNRows (), a.getNRows ());
   SX_CHECK ((int)b.getNCols () == b.getNCols (), b.getNCols ());
   SX_CHECK ((int)b.getColStride () == b.getColStride (), b.getColStride ());
   ssize_t nRow = a.getNRows (), nCol = b.getNCols ();
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * nRow * nCol);
   T *work = new T[1024];
   for (int jj = 0; jj < a.getNCols (); jj += 1024)  {
      ssize_t nj = a.getNCols () - jj;
      if (nj > 1024) nj = 1024;
      for (int ii = 0; ii < a.getNRows (); ii += 1024)  {
         ssize_t ni = a.getNRows () - ii;
         if (ni > 1024) ni = 1024;
         SxCacheOblivious task (uint16_t(ni), uint16_t(nj), 1024);
         while (task.splitProblem ())  {
            ssize_t offi = ii + task.offX;
            ssize_t offj = jj + task.offY;
            int wStride = task.nX;
            ssize_t rStride = a.getRowStride ();
            // collect a small submatrix of b
            for (uint16_t j = 0; j < task.nY; j++)  {
               T *src = a.elements + offi * rStride + (offj + j) * a.getColStride ();
               for (uint16_t i = 0; i < task.nX; i++)
                  work[i + wStride * j] = src[i * rStride];
            }
            // multiply this submatrix (from left) with a into result
            // via BLAS {s,d,c,z}gemm wrapper
            ::matmult ((T*)mem->mem + offi, (offj == 0) ? 0. : 1.,
                       work, b.elements + offj,
                       task.nX, task.nY, (int)nCol,
                       wStride, (int)b.getColStride (), (int)nRow);

         }
      }
   }
   delete [] work;
   return SxVector<T> (mem, SxVecLayout<Compact> (nRow, nCol),
                       a.auxData, a.auxData); /* keep auxData from a ! */
}

/// Matrix-matrix multiplication (first one is packed)
template<class T, StrideType Layout2>
SxVector<T>
SxVecCompute<SimdExpand>::matmult (const SxVecRef<T,PackedSymMatrix> &a,
                                   const SxVecRef<T,Layout2> &b)
{
   SX_CHECK (a.getNCols () == b.getNRows (), a.getNCols (), b.getNRows ());
#ifndef NDEBUG
   if (b.getNCols () > 1)  {
      cout << "Warning: inefficient implementation of PackedSymMatrix - matrix multiply!"
           << endl;
      // TODO: due a clever one by cache-oblivious approach + unpacking + gemm / hemm
   }
#endif
   ssize_t nRow = a.getNRows (), nCol = b.getNCols ();
   SX_CHECK ((int)nRow == nRow, nRow);
   SX_CHECK ((int)b.getRowStride () == b.getRowStride (), b.getRowStride ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * nRow * nCol);
   for (ssize_t iCol = 0; iCol < nCol; iCol++)
      symmatmult ((int)nRow, (T*)mem->mem + iCol * nRow, a.elements, UpperRight,
                  b.elements + iCol*b.getColStride (), (int)b.getRowStride ());

   return SxVector<T> (mem, SxVecLayout<Compact> (nRow, nCol),
                       b.auxData, b.auxData); /* keep auxData from b ! */
}

/// Matrix-matrix multiplication (first one is packed)
template<class T, StrideType Layout2>
SxVector<T>
SxVecCompute<ItExpand>::matmult (const SxVecRef<T,PackedSymMatrix> &a,
                                 const SxVecRef<T,Layout2> &b)
{
   SX_CHECK (a.getNCols () == b.getNRows (), a.getNCols (), b.getNRows ());
   SX_CHECK (!isPacked (Layout2));
#ifndef NDEBUG
   if (b.getNCols () > 1)  {
      cout << "Warning: inefficient implementation of PackedSymMatrix - matrix multiply!"
           << endl;
      // TODO: due a clever one by cache-oblivious approach + unpacking + gemm / hemm
   }
#endif
   ssize_t nRow = a.getNRows (), nCol = b.getNCols ();
   SX_CHECK ((int)nRow == nRow, nRow);
   SX_CHECK ((int)b.getRowStride () == b.getRowStride (), b.getRowStride ());
   SxAllocMem *mem = SxAllocMem::create (sizeof(T) * nRow * nCol);
   T* elements = (T*)mem->mem;
   for (ssize_t iCol = 0; iCol < nCol; iCol++)
      symmatmult ((int)nRow, (T*)mem->mem + iCol * nRow, a.elements, UpperRight,
                  b.elements + iCol*b.getColStride (), (int)b.getRowStride ());

   return SxVector<T> (mem, SxVecLayout<Compact> (nRow, nCol),
                       b.auxData, b.auxData); /* keep auxData from b ! */
}

// should be last in header: explicit instantiations
#include<SxVecCompute_instantiate.h>
#endif // _SX_VEC_COMPUTE_HPP
