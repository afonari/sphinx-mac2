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

#ifndef _SX_VEC_COMPUTE_H_
#define _SX_VEC_COMPUTE_H_

#include <SxMath.h>
#include <SxConfig.h>
#ifdef USE_OPENMP
#   include <omp.h>
#endif
#include <type_traits>
#include <SxBlasLib.h>

template<class, StrideType> class SxVecRef;
template<class T> class SxVector;

/// Possible vector-vector algorithms
enum SxVecAlgorithm {
   /// Use CPU vectorization (for compact storage, no type mixing)
   Simd,
   /// Use CPU vectorization, but loop over columns (one vector is SubMatrix)
   SimdCol,
   /// Use CPU vectorization, but one vector is non-contiguous
   SimdIt,
   /// Use CPU vectorization: SubMatrix and non-contiguous
   SimdItCol,  // unused
   /// Use CPU vectorization, but expand packed case
   SimdExpand, // not implemented
   /// Use Iterators (non-contiguous vectors, or type mixing)
   UseIterator,
   /// Use Iterators and expand packed cases
   ItExpand /* not implemented */ };
template<SxVecAlgorithm> class SxVecCompute;

// --- select suitable algorithm for vector-vector operations
SxVecAlgorithm constexpr algoSelect (bool sameType, StrideType layout1,
                                                    StrideType layout2)
{
   if (layout1 == layout2)  {
      if (isContiguous(layout1))
         return sameType ? Simd : UseIterator;
      else if (layout1 == SubMatrix)
         return sameType ? SimdCol : UseIterator;
      else
         return UseIterator;
   }
   // sort layouts
   if (layout1 > layout2)  {
      StrideType l2 = layout2;
      layout2 = layout1;
      layout1 = l2;
   }
   if (isPacked(layout2))  {
      // --- packed layouts combined with unpacked ones...
      if (layout1 == Compact || layout1 == SubMatrix)
         return sameType ? SimdExpand : ItExpand;
      else
         return ItExpand;
   } else if (!sameType)  {
      // --- any type mixing? => work with iterators
      return UseIterator;
   } else if (layout2 == SubMatrix)  {
      if (layout1 == Compact)
         return SimdCol;
      else // layout1 == Strided
         //return SimdItCol;
         return UseIterator;
   } else if (layout1 == Compact)  {
      // (layout2 == Strided || layout2 = GeneralMatrix)
      return SimdIt;
   } else if (layout1 == SubMatrix && layout2 == GeneralMatrix)  {
      //return SimdItCol;
      return UseIterator;
   } else {
      return UseIterator;
   }
}

// --- select suitable algorithm for vector-vector operations
SxVecAlgorithm constexpr algoInPlace (StrideType layout1, StrideType layout2)
{
   if (layout1 == layout2)  {
      if (isContiguous(layout1))
         return Simd;
      else
         return UseIterator;
   }
   // sort layouts
   if (layout1 > layout2)  {
      StrideType l2 = layout2;
      layout2 = layout1;
      layout1 = l2;
   }
   if (isPacked(layout2))  {
      // --- packed layouts combined with unpacked ones...
      if (layout1 == Compact || layout1 == SubMatrix)
         return SimdExpand;
      else
         return ItExpand;
   }
   return UseIterator;
}

// --- select suitable algorithm for vector-scalar operations
SxVecAlgorithm constexpr algoScalar (StrideType layout)
{
   if (isContiguous(layout))
      return Simd;
   else
      return UseIterator;
}

// --- select suitable algorithm for matrix transpose operations
SxVecAlgorithm constexpr algoTrans (StrideType layout)
{
   return (layout == Compact || layout == SubMatrix) ? Simd : SimdIt;
}

// --- select suitable algorithm for matrix multiplication
SxVecAlgorithm constexpr algoMat (StrideType layout1, StrideType layout2)
{
   // Compact and SubMatrix fit to BLAS's "general matrix" layout
   bool blasMat1 = (layout1 == Compact || layout1 == SubMatrix);
   bool blasMat2 = (layout2 == Compact || layout2 == SubMatrix);
   if (blasMat1 && blasMat2) return Simd;
   bool strided2 = (layout2 == Strided || layout2 == GeneralMatrix);
   if (blasMat1 && strided2) return SimdCol;
   bool strided1 = (layout1 == Strided || layout1 == GeneralMatrix);
   if (strided1 && blasMat2) return SimdIt;
   if (strided1 && strided2) return UseIterator;
   // -> at this point, one of the two matrices is packed
   // use SimdExpand for packed/"general"
   // and ItExpand for packed/strided
   // and compiler error for packed/packed
   return (blasMat1 || blasMat2) ? SimdExpand : ItExpand;
}

// --- map result type of (some) vector operations:
//     packed formats (PackedSymMatrix, BandedMatrix) stay packed
//     everything else becomes SxVector
template<class T, StrideType Layout>
class SxVecResult {
   public:
      typedef SxVector<T> VecType;
};

template<class T>
class SxVecResult<T, PackedSymMatrix> {
   public:
      typedef SxVecRef<T, PackedSymMatrix> VecType;
};

template<class T>
class SxVecResult<T, BandedMatrix> {
   public:
      typedef SxVecRef<T, BandedMatrix> VecType;
};


/** \brief Container for Vector operations:
     SIMD for contiguous memory layouts: Compact, PackedSymMatrix, BandedMatrix

  */
template<>
class SxVecCompute<Simd>
{
   public:
      // --- The next 60 lines were generated from snippets/SxVecCompute.h snippet Simd_OP
      /// add two vectors
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      add (const SxVecRef<T, Layout> &x, const SxVecRef<T, Layout> &y);

      /// add compute routine
      template<class T, class T2>
      static void add (decltype(T(0) + T2(0)) *resPtr, ssize_t n,
                        const T* x, const T2* y);

      // vector x + y (x, y are contiguous, y is about to die)
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      add (const SxVecRef<T, Layout> &x, SxVecRef<T, Layout> &&y);

      /// subtract two vectors
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      subtract (const SxVecRef<T, Layout> &x, const SxVecRef<T, Layout> &y);

      /// subtract compute routine
      template<class T, class T2>
      static void subtract (decltype(T(0) - T2(0)) *resPtr, ssize_t n,
                        const T* x, const T2* y);

      // vector x - y (x, y are contiguous, y is about to die)
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      subtract (const SxVecRef<T, Layout> &x, SxVecRef<T, Layout> &&y);

      /// multiply two vectors
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      multiply (const SxVecRef<T, Layout> &x, const SxVecRef<T, Layout> &y);

      /// multiply compute routine
      template<class T, class T2>
      static void multiply (decltype(T(0) * T2(0)) *resPtr, ssize_t n,
                        const T* x, const T2* y);

      // vector x * y (x, y are contiguous, y is about to die)
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      multiply (const SxVecRef<T, Layout> &x, SxVecRef<T, Layout> &&y);

      /// divide two vectors
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      divide (const SxVecRef<T, Layout> &x, const SxVecRef<T, Layout> &y);

      /// divide compute routine
      template<class T, class T2>
      static void divide (decltype(T(0) / T2(0)) *resPtr, ssize_t n,
                        const T* x, const T2* y);

      // vector x / y (x, y are contiguous, y is about to die)
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      divide (const SxVecRef<T, Layout> &x, SxVecRef<T, Layout> &&y);
      // --- Simd_OP

      /** \brief Multiply columns of a matrix with the same vector
          @param resPtr destination (n x m)
          @param n      number of rows
          @param m      number of columns of matrix
          @param x      source matrix (n x m)
          @param y      source vector (n x 1)
      */
      template<class T1, class T2>
      static inline
      void multiplyMx1 (decltype(T1(0) * T2(0)) *resPtr,
                        /* ssize_t colStrideRes, */
                        ssize_t n, ssize_t m,
                        const T1* x,
                        ssize_t colStrideX,
                        const T2* y);

      // --- The next 17 lines were generated from snippets/SxVecCompute.h snippet Simd_OPEQ
      /// vector x+= y
      template<class T, StrideType Layout>
      static void addInPlace (      SxVecRef<T, Layout> &x,
                              const SxVecRef<T, Layout> &y);
      /// vector x-= y
      template<class T, StrideType Layout>
      static void subtractInPlace (      SxVecRef<T, Layout> &x,
                              const SxVecRef<T, Layout> &y);
      /// vector x*= y
      template<class T, StrideType Layout>
      static void multiplyInPlace (      SxVecRef<T, Layout> &x,
                              const SxVecRef<T, Layout> &y);
      /// vector x/= y
      template<class T, StrideType Layout>
      static void divideInPlace (      SxVecRef<T, Layout> &x,
                              const SxVecRef<T, Layout> &y);
      // --- Simd_OPEQ

      /// y += a * x
      template<class T, StrideType Layout>
      static inline void axpy(SxVecRef<T,Layout> &y, const T &a,
                              const SxVecRef<T, Layout> &x)
      {
         ::axpy (y.elements, a, x.elements, x.getSize ());
      }

      // --- The next 27 lines were generated from snippets/SxVecCompute.h snippet Simd_OP_scalar
      /// add vector and scalar
      template<class T, StrideType Layout, class T2>
      static typename SxVecResult<decltype(T(0)+T2(0)),Layout>::VecType
      add (const SxVecRef<T, Layout> &x, const T2 &y);

      /// add vector and scalar (in place)
      template<class T, StrideType Layout, class T2>
      static void addInPlace (SxVecRef<T, Layout> &x, const T2 &y);

      /// subtract vector and scalar
      template<class T, StrideType Layout, class T2>
      static typename SxVecResult<decltype(T(0)+T2(0)),Layout>::VecType
      subtract (const SxVecRef<T, Layout> &x, const T2 &y);

      /// subtract vector and scalar (in place)
      template<class T, StrideType Layout, class T2>
      static void subtractInPlace (SxVecRef<T, Layout> &x, const T2 &y);

      /// multiply vector and scalar
      template<class T, StrideType Layout, class T2>
      static typename SxVecResult<decltype(T(0)+T2(0)),Layout>::VecType
      multiply (const SxVecRef<T, Layout> &x, const T2 &y);

      /// multiply vector and scalar (in place)
      template<class T, StrideType Layout, class T2>
      static void multiplyInPlace (SxVecRef<T, Layout> &x, const T2 &y);
      // --- Simd_OP_scalar

      // vector * scalar (vector contiguous)
      template<class T, class T2>
      static void multiply (decltype(T(0) + T2(0)) *resPtr, ssize_t n,
                            const T *x, const T2 &y);

      // --- divide can be done via inverse if T is not any integer type
      //     For integer types, the division must be done for every elementsent
      //     TODO: there exists an algorithm for fast integer division via
      //     a bit-shifted inverse (multiply with 2^N / y; then right-shift
      //     by N bits)
      template<bool NoInverse>
      class Divide;

      /// divide vector and scalar
      template<class T, StrideType Layout,class T2>
      static typename SxVecResult<decltype(T(0)+T2(0)),Layout>::VecType
      divide (const SxVecRef<T, Layout> &x, const T2 &y)
      {
         return Divide<std::is_integral<T2>::value>::divide (x,y);
      }

      /// divide vector and scalar (in place)
      template<class T, StrideType Layout, class T2>
      static void
      divideInPlace (SxVecRef<T, Layout> &x, const T2 &y)
      {
         Divide<std::is_integral<T2>::value>::divideInPlace (x,y);
      }

      // --- The next 10 lines were generated from snippets/SxVecCompute.h snippet Simd_scalar_vec
      /// compute scalar - vector
      template<class T, StrideType Layout, class T2>
      static typename SxVecResult<decltype(T(0)+T2(0)),Layout>::VecType
      subtract (const T2 &y, const SxVecRef<T, Layout> &x);

      /// compute scalar / vector
      template<class T, StrideType Layout, class T2>
      static typename SxVecResult<decltype(T(0)+T2(0)),Layout>::VecType
      divide (const T2 &y, const SxVecRef<T, Layout> &x);
      // --- Simd_scalar_vec

      /// Wrap sqrt to float-specific function
      static inline float sqrt(float x) { return ::sqrtf(x); }
      /// Wrap sqrt to double-specific function
      static inline double sqrt(double x) { return ::sqrt(x); }

      // --- The next 80 lines were generated from snippets/SxVecCompute.h snippet Simd_vecfunc
      /// Unary minus
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      minus (const SxVecRef<T,Layout> &x);

      /// Unary minus (in-place, if feasible)
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      minusInPlace (SxVecRef<T,Layout> &&);

      /// Square
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      sqr (const SxVecRef<T,Layout> &x);

      /// Square (in-place, if feasible)
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      sqrInPlace (SxVecRef<T,Layout> &&);

      /// Cube (x*x*x)
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      cub (const SxVecRef<T,Layout> &x);

      /// Cube (x*x*x) (in-place, if feasible)
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      cubInPlace (SxVecRef<T,Layout> &&);

      /// Exponential
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      exp (const SxVecRef<T,Layout> &x);

      /// Exponential (in-place, if feasible)
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      expInPlace (SxVecRef<T,Layout> &&);

      /// Square root
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      sqrt (const SxVecRef<T,Layout> &x);

      /// Square root (in-place, if feasible)
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      sqrtInPlace (SxVecRef<T,Layout> &&);

      /// Logarithm
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      log (const SxVecRef<T,Layout> &x);

      /// Logarithm (in-place, if feasible)
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      logInPlace (SxVecRef<T,Layout> &&);

      /// Error function
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      erf (const SxVecRef<T,Layout> &x);

      /// Error function (in-place, if feasible)
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      erfInPlace (SxVecRef<T,Layout> &&);

      /// Complementary error function
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      erfc (const SxVecRef<T,Layout> &x);

      /// Complementary error function (in-place, if feasible)
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      erfcInPlace (SxVecRef<T,Layout> &&);
      // --- Simd_vecfunc
      // --- The next 33 lines were generated from snippets/SxVecCompute.h snippet Simd_abs_real
      /// Absolute value
      static SxVector<int> abs (const SxVecRef<int,Compact> &);
      /// Absolute value (in-place, if feasible)
      static SxVector<int> absInPlace (SxVecRef<int,Compact> &&);
      /// Absolute value squared
      static SxVector<int> absSqr (const SxVecRef<int,Compact> &);
      /// Absolute value squared (in-place, if feasible)
      static SxVector<int> absSqrInPlace (SxVecRef<int,Compact> &&);
      /// Absolute value
      static SxVector<long> abs (const SxVecRef<long,Compact> &);
      /// Absolute value (in-place, if feasible)
      static SxVector<long> absInPlace (SxVecRef<long,Compact> &&);
      /// Absolute value squared
      static SxVector<long> absSqr (const SxVecRef<long,Compact> &);
      /// Absolute value squared (in-place, if feasible)
      static SxVector<long> absSqrInPlace (SxVecRef<long,Compact> &&);
      /// Absolute value
      static SxVector<float> abs (const SxVecRef<float,Compact> &);
      /// Absolute value (in-place, if feasible)
      static SxVector<float> absInPlace (SxVecRef<float,Compact> &&);
      /// Absolute value squared
      static SxVector<float> absSqr (const SxVecRef<float,Compact> &);
      /// Absolute value squared (in-place, if feasible)
      static SxVector<float> absSqrInPlace (SxVecRef<float,Compact> &&);
      /// Absolute value
      static SxVector<double> abs (const SxVecRef<double,Compact> &);
      /// Absolute value (in-place, if feasible)
      static SxVector<double> absInPlace (SxVecRef<double,Compact> &&);
      /// Absolute value squared
      static SxVector<double> absSqr (const SxVecRef<double,Compact> &);
      /// Absolute value squared (in-place, if feasible)
      static SxVector<double> absSqrInPlace (SxVecRef<double,Compact> &&);
      // --- Simd_abs_real
      // --- The next 26 lines were generated from snippets/SxVecCompute.h snippet Simd_abs_complex
      /// Absolute value
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      abs (const SxVecRef<SxComplex<T>,Layout> &x);

      /// Absolute value (cannot be done in-place, wraps abs)
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      absInPlace (SxVecRef<SxComplex<T>,Layout> &&x)
      {
         return abs(x);
      }

      /// Absolute value squared
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      absSqr (const SxVecRef<SxComplex<T>,Layout> &x);

      /// Absolute value squared (cannot be done in-place, wraps absSqr)
      template<class T, StrideType Layout>
      static inline typename SxVecResult<T,Layout>::VecType
      absSqrInPlace (SxVecRef<SxComplex<T>,Layout> &&x)
      {
         return absSqr(x);
      }
      // --- Simd_abs_complex

      // --- specialized implementations of reductions (sum/product/normSqr)
      /// Sum elements
      static double sum (const double*, ssize_t n);
      /// Sum elements
      static SxComplex16 sum (const SxComplex16*, ssize_t n);

      /// Sum square of elements
      static double normSqr (const double*, ssize_t n);

      /// Add abs. squares of elements
      template<class T>
      static inline T normSqr (const SxComplex<T> *elements, ssize_t n)
      {
         return normSqr((const T*)elements, 2 * n);
      }

      /// Add abs. squares of elements (generic implementation)
      template<class T>
      static T normSqr (const T *elements, ssize_t n);

      // --- The next 24 lines were generated from ./snippets/SxVecCompute.h snippet Simd_reduce
      /// Add elements
      template<class T>
      static inline T sum (const T*, ssize_t n);

      /// Add all elements
      template<class T, StrideType Layout>
      static inline T
      sum (const SxVecRef<T,Layout> &x)
      {
         return sum(x.elements, x.getSize ());
      }

      /// Multiply elements
      template<class T>
      static inline T product (const T*, ssize_t n);

      /// Multiply all elements
      template<class T, StrideType Layout>
      static inline T
      product (const SxVecRef<T,Layout> &x)
      {
         return product(x.elements, x.getSize ());
      }
      // --- Simd_reduce

      /// Euclidean norm of vector
      template<class T, StrideType Layout>
      static inline typename SxTypeMapper<T>::TReal
      norm (const SxVecRef<T,Layout> &x)
      {
         return norm2 (x.elements, x.getSize ());
      }

      /// Add abs. squares of all elements
      template<class T, StrideType Layout>
      static inline typename SxTypeMapper<T>::TReal
      normSqr (const SxVecRef<T,Layout> &x)
      {
         // TODO: test if this is faster than our own normSqr
         typename SxTypeMapper<T>::TReal nrm = norm (x);
         return nrm * nrm;
      }

      /// Real part for packed formats
      template<class T, StrideType Layout>
      static SxVecRef<T, Layout>
      real (const SxVecRef<SxComplex<T>, Layout> &vec);

      /// Imaginary part for packed formats
      template<class T, StrideType Layout>
      static SxVecRef<T, Layout>
      imag (const SxVecRef<SxComplex<T>, Layout> &vec);

      /// Conjugate (only complex!)
      template<class T, StrideType Layout>
      static typename SxVecResult<SxComplex<T>, Layout>::VecType
      conj (const SxVecRef<SxComplex<T>, Layout> &vec);

      /// Conjugate in place (only complex!)
      template<class T, StrideType Layout>
      static typename SxVecResult<SxComplex<T>, Layout>::VecType
      conjInPlace (SxVecRef<SxComplex<T>, Layout> &&vec);

      // --- The next 20 lines were generated from snippets/SxVecCompute.h snippet Simd_transpose
      /// Matrix transpose
      template<class T>
      static void transpose (ssize_t N, ssize_t M,
                             const T* src, ssize_t strideS,
                             T* dest, ssize_t strideD);

      /// Matrix transpose
      template<class T, StrideType Layout> // Layout = Compact,SubMatrix
      static inline SxVector<T > transpose (const SxVecRef<T, Layout> &vec);

      /// Matrix adjoint
      template<class T>
      static void adjoint (ssize_t N, ssize_t M,
                             const SxComplex<T>* src, ssize_t strideS,
                             SxComplex<T>* dest, ssize_t strideD);

      /// Matrix adjoint
      template<class T, StrideType Layout> // Layout = Compact,SubMatrix
      static inline SxVector<SxComplex<T> > adjoint (const SxVecRef<SxComplex<T>, Layout> &vec);
      // --- Simd_transpose

      /// Matrix "adjoint" for real-type matrix = transpose
      static inline void adjoint (ssize_t N, ssize_t M,
                                  const float* src, ssize_t strideS,
                                  float* dest, ssize_t strideD)
      {
         transpose (N, M, src, strideS, dest, strideD);
      }

      /// Matrix "adjoint" for real-type matrix = transpose
      static inline void adjoint (ssize_t N, ssize_t M,
                                  const double* src, ssize_t strideS,
                                  double* dest, ssize_t strideD)
      {
         transpose (N, M, src, strideS, dest, strideD);
      }

      /// Matrix-matrix multiplication
      template<class T, StrideType Layout1, StrideType Layout2>
      static inline SxVector<T> matmult (const SxVecRef<T,Layout1> &a,
                                         const SxVecRef<T,Layout2> &b);
      /// check whether a matrix is Hermitian
      static bool isHermitian (const SxVector<float> &H);
      static bool isHermitian (const SxVector<double> &H);
      static bool isHermitian (const SxVector<SxComplex8> &H);
      static bool isHermitian (const SxVector<SxComplex16> &H);
};

/** \brief Container for Vector operations:
     SIMD for contiguous column layouts: Compact / SubMatrix

  */
template<>
class SxVecCompute<SimdCol>
{
   public:
      // --- The next 20 lines were generated from snippets/SxVecCompute.h snippet SimdCol_OP
      /// add two vectors
      template<class T, StrideType Layout, StrideType Layout2>
      static SxVector<T>
      add (const SxVecRef<T, Layout> &x, const SxVecRef<T, Layout2> &y);

      /// subtract two vectors
      template<class T, StrideType Layout, StrideType Layout2>
      static SxVector<T>
      subtract (const SxVecRef<T, Layout> &x, const SxVecRef<T, Layout2> &y);

      /// multiply two vectors
      template<class T, StrideType Layout, StrideType Layout2>
      static SxVector<T>
      multiply (const SxVecRef<T, Layout> &x, const SxVecRef<T, Layout2> &y);

      /// divide two vectors
      template<class T, StrideType Layout, StrideType Layout2>
      static SxVector<T>
      divide (const SxVecRef<T, Layout> &x, const SxVecRef<T, Layout2> &y);
      // --- SimdCol_OP

      // --- The next 15 lines were generated from snippets/SxVecCompute.h snippet SimdCol_reduce
      /// Add all elements
      template<class T>
      static inline T
      sum (const SxVecRef<T,SubMatrix> &x);

      /// Multiply all elements
      template<class T>
      static inline T
      product (const SxVecRef<T,SubMatrix> &x);

      /// Add abs. squares of all elements
      template<class T>
      static inline typename SxTypeMapper<T>::TReal
      normSqr (const SxVecRef<T,SubMatrix> &x);
      // --- SimdCol_reduce

      /// Euclidean norm
      template<class T>
      static inline typename SxTypeMapper<T>::TReal
      norm (const SxVecRef<T,SubMatrix> &x)
      {
         return typename SxTypeMapper<T>::TReal(sqrt(x.normSqr ()));
      }

      /// Matrix-matrix multiplication (first matrix is col-compact, 2nd one row-strided)
      template<class T, StrideType Layout1, StrideType Layout2>
      static SxVector<T> matmult (const SxVecRef<T,Layout1> &a,
                                  const SxVecRef<T,Layout2> &b);
};

/** \brief Container for Vector operations:
     One operand is Compact, one has strided layout

  */
template<>
class SxVecCompute<SimdIt>
{
   public:
      // --- The next 80 lines were generated from snippets/SxVecCompute.h snippet SimdIt_OP
      /// add two vectors
      template<class T, StrideType Layout>
      static SxVector<T>
      add (const SxVecRef<T, Compact> &x, const SxVecRef<T, Layout> &y);

      /// add two vectors: compute routine
      template<class T, class T2, StrideType Layout>
      static inline void add (const T* xPtr, SxVecConstIt<T2, Layout> it,
                              ssize_t n, decltype(T() + T2()) *resPtr);

      /// add two vectors
      template<class T, StrideType Layout>
      static SxVector<T>
      add (const SxVecRef<T, Layout> &x, const SxVecRef<T, Compact> &y);

      /// add two vectors: compute routine
      template<class T, class T2, StrideType Layout>
      static inline void add (SxVecConstIt<T, Layout> it, const T2* xPtr,
                              ssize_t n, decltype(T() + T2()) *resPtr);

      /// subtract two vectors
      template<class T, StrideType Layout>
      static SxVector<T>
      subtract (const SxVecRef<T, Compact> &x, const SxVecRef<T, Layout> &y);

      /// subtract two vectors: compute routine
      template<class T, class T2, StrideType Layout>
      static inline void subtract (const T* xPtr, SxVecConstIt<T2, Layout> it,
                              ssize_t n, decltype(T() - T2()) *resPtr);

      /// subtract two vectors
      template<class T, StrideType Layout>
      static SxVector<T>
      subtract (const SxVecRef<T, Layout> &x, const SxVecRef<T, Compact> &y);

      /// subtract two vectors: compute routine
      template<class T, class T2, StrideType Layout>
      static inline void subtract (SxVecConstIt<T, Layout> it, const T2* xPtr,
                              ssize_t n, decltype(T() - T2()) *resPtr);

      /// multiply two vectors
      template<class T, StrideType Layout>
      static SxVector<T>
      multiply (const SxVecRef<T, Compact> &x, const SxVecRef<T, Layout> &y);

      /// multiply two vectors: compute routine
      template<class T, class T2, StrideType Layout>
      static inline void multiply (const T* xPtr, SxVecConstIt<T2, Layout> it,
                              ssize_t n, decltype(T() * T2()) *resPtr);

      /// multiply two vectors
      template<class T, StrideType Layout>
      static SxVector<T>
      multiply (const SxVecRef<T, Layout> &x, const SxVecRef<T, Compact> &y);

      /// multiply two vectors: compute routine
      template<class T, class T2, StrideType Layout>
      static inline void multiply (SxVecConstIt<T, Layout> it, const T2* xPtr,
                              ssize_t n, decltype(T() * T2()) *resPtr);

      /// divide two vectors
      template<class T, StrideType Layout>
      static SxVector<T>
      divide (const SxVecRef<T, Compact> &x, const SxVecRef<T, Layout> &y);

      /// divide two vectors: compute routine
      template<class T, class T2, StrideType Layout>
      static inline void divide (const T* xPtr, SxVecConstIt<T2, Layout> it,
                              ssize_t n, decltype(T() / T2()) *resPtr);

      /// divide two vectors
      template<class T, StrideType Layout>
      static SxVector<T>
      divide (const SxVecRef<T, Layout> &x, const SxVecRef<T, Compact> &y);

      /// divide two vectors: compute routine
      template<class T, class T2, StrideType Layout>
      static inline void divide (SxVecConstIt<T, Layout> it, const T2* xPtr,
                              ssize_t n, decltype(T() / T2()) *resPtr);
      // --- SimdIt_OP

      // --- The next 48 lines were generated from snippets/SxVecCompute.h snippet SimdIt_transpose
      /// Matrix transpose
      template<class T>
      static void transpose (ssize_t N, ssize_t M,
                             const T* src, ssize_t rStride, ssize_t cStride,
                             T* dest, ssize_t strideD);

      /// Matrix transpose
      template<class T, StrideType Layout> // Layout = Strided, GeneralMatrix
      static inline SxVector<T > transpose (const SxVecRef<T, Layout> &vec);

      /// Matrix transpose
      template<class T>
      static SxVector<T > transpose (const SxVecRef<T, BandedMatrix> &vec)
      {
         SX_EXIT; // not implemented
      }

      /// Matrix transpose
      template<class T>
      static SxVector<T> transpose (const SxVecRef<T, PackedSymMatrix> &vec)
      {
         SX_EXIT; // nonsense
      }

      /// Matrix adjoint
      template<class T>
      static void adjoint (ssize_t N, ssize_t M,
                             const SxComplex<T>* src, ssize_t rStride, ssize_t cStride,
                             SxComplex<T>* dest, ssize_t strideD);

      /// Matrix adjoint
      template<class T, StrideType Layout> // Layout = Strided, GeneralMatrix
      static inline SxVector<SxComplex<T> > adjoint (const SxVecRef<SxComplex<T>, Layout> &vec);

      /// Matrix adjoint
      template<class T>
      static SxVector<SxComplex<T> > adjoint (const SxVecRef<SxComplex<T>, BandedMatrix> &vec)
      {
         SX_EXIT; // not implemented
      }

      /// Matrix adjoint
      template<class T>
      static SxVector<SxComplex<T>> adjoint (const SxVecRef<SxComplex<T>, PackedSymMatrix> &vec)
      {
         SX_EXIT; // nonsense
      }
      // --- SimdIt_transpose

      /// Matrix-matrix multiplication (first matrix is row-strided, 2nd one col-compact)
      template<class T, StrideType Layout1, StrideType Layout2>
      static inline SxVector<T> matmult (const SxVecRef<T,Layout1> &a,
                                         const SxVecRef<T,Layout2> &b);
};

/** \brief Container for Vector operations:
     Use Iterators (non-contiguous vectors, or type mixing)

  */
template<>
class SxVecCompute<UseIterator>
{
   public:
      // --- The next 20 lines were generated from snippets/SxVecCompute.h snippet UseIterator_OP
      /// add two vectors
      template<class T, StrideType Layout, class T2, StrideType Layout2>
      static SxVector<decltype(T() + T2())>
      add (const SxVecRef<T, Layout>  &x, const SxVecRef<T2, Layout2> &y);

      /// subtract two vectors
      template<class T, StrideType Layout, class T2, StrideType Layout2>
      static SxVector<decltype(T() + T2())>
      subtract (const SxVecRef<T, Layout>  &x, const SxVecRef<T2, Layout2> &y);

      /// multiply two vectors
      template<class T, StrideType Layout, class T2, StrideType Layout2>
      static SxVector<decltype(T() + T2())>
      multiply (const SxVecRef<T, Layout>  &x, const SxVecRef<T2, Layout2> &y);

      /// divide two vectors
      template<class T, StrideType Layout, class T2, StrideType Layout2>
      static SxVector<decltype(T() + T2())>
      divide (const SxVecRef<T, Layout>  &x, const SxVecRef<T2, Layout2> &y);
      // --- UseIterator_OP

      // --- The next 17 lines were generated from snippets/SxVecCompute.h snippet UseIterator_OPEQ
      /// vector x+= y
      template<class T, StrideType Layout, class T2, StrideType Layout2>
      static void addInPlace (      SxVecRef<T, Layout> &x,
                              const SxVecRef<T2, Layout2> &y);
      /// vector x-= y
      template<class T, StrideType Layout, class T2, StrideType Layout2>
      static void subtractInPlace (      SxVecRef<T, Layout> &x,
                              const SxVecRef<T2, Layout2> &y);
      /// vector x*= y
      template<class T, StrideType Layout, class T2, StrideType Layout2>
      static void multiplyInPlace (      SxVecRef<T, Layout> &x,
                              const SxVecRef<T2, Layout2> &y);
      /// vector x/= y
      template<class T, StrideType Layout, class T2, StrideType Layout2>
      static void divideInPlace (      SxVecRef<T, Layout> &x,
                              const SxVecRef<T2, Layout2> &y);
      // --- UseIterator_OPEQ

      /// y += a * x
      template<class T, StrideType Layout, StrideType Layout2>
      static inline void axpy(SxVecRef<T,Layout> &y, const T &a,
                              const SxVecRef<T, Layout2> &x);
      // y += a * x
      template<class T,StrideType Layout,class Ta,class T2,StrideType Layout2>
      static void axpy(SxVecRef<T,Layout> &y, const Ta &a,
                       const SxVecRef<T2, Layout2> &x);

      // --- The next 27 lines were generated from snippets/SxVecCompute.h snippet UseIterator_OP_scalar
      /// add vector and scalar
      template<class T, StrideType Layout, class T2>
      static SxVector<decltype(T(0)+T2(0))>
      add (const SxVecRef<T, Layout> &x, const T2 &y);

      /// add vector and scalar (in place)
      template<class T, StrideType Layout, class T2>
      static void addInPlace (SxVecRef<T, Layout> &x, const T2 &y);

      /// subtract vector and scalar
      template<class T, StrideType Layout, class T2>
      static SxVector<decltype(T(0)+T2(0))>
      subtract (const SxVecRef<T, Layout> &x, const T2 &y);

      /// subtract vector and scalar (in place)
      template<class T, StrideType Layout, class T2>
      static void subtractInPlace (SxVecRef<T, Layout> &x, const T2 &y);

      /// multiply vector and scalar
      template<class T, StrideType Layout, class T2>
      static SxVector<decltype(T(0)+T2(0))>
      multiply (const SxVecRef<T, Layout> &x, const T2 &y);

      /// multiply vector and scalar (in place)
      template<class T, StrideType Layout, class T2>
      static void multiplyInPlace (SxVecRef<T, Layout> &x, const T2 &y);
      // --- UseIterator_OP_scalar

      // --- divide can be done via inverse if T is not any integer type
      //     For integer types, the division must be done for every elementsent
      //     TODO: there exists an algorithm for fast integer division via
      //     a bit-shifted inverse (multiply with 2^N / y; then right-shift
      //     by N bits)
      template<bool NoInverse>
      class Divide;

      /// divide vector and scalar
      template<class T, StrideType Layout,class T2>
      static SxVector<decltype(T(0)+T2(0))>
      divide (const SxVecRef<T, Layout> &x, const T2 &y)
      {
         return Divide<std::is_integral<T2>::value>::divide (x,y);
      }

      /// divide vector and scalar (in place)
      template<class T, StrideType Layout, class T2>
      static void
      divideInPlace (SxVecRef<T, Layout> &x, const T2 &y)
      {
         Divide<std::is_integral<T2>::value>::divideInPlace (x,y);
      }

      // --- The next 10 lines were generated from snippets/SxVecCompute.h snippet UseIterator_scalar_vec
      /// compute scalar - vector
      template<class T, StrideType Layout, class T2>
      static SxVector<decltype(T(0)+T2(0))>
      subtract (const T2 &y, const SxVecRef<T, Layout> &x);

      /// compute scalar / vector
      template<class T, StrideType Layout, class T2>
      static SxVector<decltype(T(0)+T2(0))>
      divide (const T2 &y, const SxVecRef<T, Layout> &x);
      // --- UseIterator_scalar_vec

      // --- The next 88 lines were generated from snippets/SxVecCompute.h snippet UseIterator_vecfunc
      /// Unary minus
      template<class T, StrideType Layout>
      static inline SxVector<T> minus (const SxVecRef<T,Layout> &x);

      /// Unary minus (cannot be done in-place, wraps minus)
      template<class T, StrideType Layout>
      static inline SxVector<T> minusInPlace (SxVecRef<T,Layout> &&x)
      {
         return minus(x);
      }

      /// Square
      template<class T, StrideType Layout>
      static inline SxVector<T> sqr (const SxVecRef<T,Layout> &x);

      /// Square (cannot be done in-place, wraps sqr)
      template<class T, StrideType Layout>
      static inline SxVector<T> sqrInPlace (SxVecRef<T,Layout> &&x)
      {
         return sqr(x);
      }

      /// Cube (x*x*x)
      template<class T, StrideType Layout>
      static inline SxVector<T> cub (const SxVecRef<T,Layout> &x);

      /// Cube (x*x*x) (cannot be done in-place, wraps cub)
      template<class T, StrideType Layout>
      static inline SxVector<T> cubInPlace (SxVecRef<T,Layout> &&x)
      {
         return cub(x);
      }

      /// Exponential
      template<class T, StrideType Layout>
      static inline SxVector<T> exp (const SxVecRef<T,Layout> &x);

      /// Exponential (cannot be done in-place, wraps exp)
      template<class T, StrideType Layout>
      static inline SxVector<T> expInPlace (SxVecRef<T,Layout> &&x)
      {
         return exp(x);
      }

      /// Square root
      template<class T, StrideType Layout>
      static inline SxVector<T> sqrt (const SxVecRef<T,Layout> &x);

      /// Square root (cannot be done in-place, wraps sqrt)
      template<class T, StrideType Layout>
      static inline SxVector<T> sqrtInPlace (SxVecRef<T,Layout> &&x)
      {
         return sqrt(x);
      }

      /// Logarithm
      template<class T, StrideType Layout>
      static inline SxVector<T> log (const SxVecRef<T,Layout> &x);

      /// Logarithm (cannot be done in-place, wraps log)
      template<class T, StrideType Layout>
      static inline SxVector<T> logInPlace (SxVecRef<T,Layout> &&x)
      {
         return log(x);
      }

      /// Error function
      template<class T, StrideType Layout>
      static inline SxVector<T> erf (const SxVecRef<T,Layout> &x);

      /// Error function (cannot be done in-place, wraps erf)
      template<class T, StrideType Layout>
      static inline SxVector<T> erfInPlace (SxVecRef<T,Layout> &&x)
      {
         return erf(x);
      }

      /// Complementary error function
      template<class T, StrideType Layout>
      static inline SxVector<T> erfc (const SxVecRef<T,Layout> &x);

      /// Complementary error function (cannot be done in-place, wraps erfc)
      template<class T, StrideType Layout>
      static inline SxVector<T> erfcInPlace (SxVecRef<T,Layout> &&x)
      {
         return erfc(x);
      }
      // --- UseIterator_vecfunc
      // --- The next 64 lines were generated from snippets/SxVecCompute.h snippet UseIterator_abs_real
      /// Absolute value
      template<StrideType Layout>
      static SxVector<int> abs (const SxVecRef<int,Layout> &);

      /// Absolute value (cannot be done in-place, wraps abs)
      template<StrideType Layout>
      static inline SxVector<int> absInPlace (SxVecRef<int,Layout> &&x);

      /// Absolute value squared
      template<StrideType Layout>
      static SxVector<int> absSqr (const SxVecRef<int,Layout> &);

      /// Absolute value squared (cannot be done in-place, wraps absSqr)
      template<StrideType Layout>
      static inline SxVector<int> absSqrInPlace (SxVecRef<int,Layout> &&x);

      /// Absolute value
      template<StrideType Layout>
      static SxVector<long> abs (const SxVecRef<long,Layout> &);

      /// Absolute value (cannot be done in-place, wraps abs)
      template<StrideType Layout>
      static inline SxVector<long> absInPlace (SxVecRef<long,Layout> &&x);

      /// Absolute value squared
      template<StrideType Layout>
      static SxVector<long> absSqr (const SxVecRef<long,Layout> &);

      /// Absolute value squared (cannot be done in-place, wraps absSqr)
      template<StrideType Layout>
      static inline SxVector<long> absSqrInPlace (SxVecRef<long,Layout> &&x);

      /// Absolute value
      template<StrideType Layout>
      static SxVector<float> abs (const SxVecRef<float,Layout> &);

      /// Absolute value (cannot be done in-place, wraps abs)
      template<StrideType Layout>
      static inline SxVector<float> absInPlace (SxVecRef<float,Layout> &&x);

      /// Absolute value squared
      template<StrideType Layout>
      static SxVector<float> absSqr (const SxVecRef<float,Layout> &);

      /// Absolute value squared (cannot be done in-place, wraps absSqr)
      template<StrideType Layout>
      static inline SxVector<float> absSqrInPlace (SxVecRef<float,Layout> &&x);

      /// Absolute value
      template<StrideType Layout>
      static SxVector<double> abs (const SxVecRef<double,Layout> &);

      /// Absolute value (cannot be done in-place, wraps abs)
      template<StrideType Layout>
      static inline SxVector<double> absInPlace (SxVecRef<double,Layout> &&x);

      /// Absolute value squared
      template<StrideType Layout>
      static SxVector<double> absSqr (const SxVecRef<double,Layout> &);

      /// Absolute value squared (cannot be done in-place, wraps absSqr)
      template<StrideType Layout>
      static inline SxVector<double> absSqrInPlace (SxVecRef<double,Layout> &&x);
      // --- UseIterator_abs_real
      // --- The next 22 lines were generated from snippets/SxVecCompute.h snippet UseIterator_abs_complex
      /// Absolute value
      template<class T, StrideType Layout>
      static inline SxVector<T> abs (const SxVecRef<SxComplex<T>,Layout> &x);

      /// Absolute value (cannot be done in-place, wraps abs)
      template<class T, StrideType Layout>
      static inline SxVector<T> absInPlace (SxVecRef<SxComplex<T>,Layout> &&x)
      {
         return abs(x);
      }

      /// Absolute value squared
      template<class T, StrideType Layout>
      static inline SxVector<T> absSqr (const SxVecRef<SxComplex<T>,Layout> &x);

      /// Absolute value squared (cannot be done in-place, wraps absSqr)
      template<class T, StrideType Layout>
      static inline SxVector<T> absSqrInPlace (SxVecRef<SxComplex<T>,Layout> &&x)
      {
         return absSqr(x);
      }
      // --- UseIterator_abs_complex

      /// Wrap absSqr for SxComplex16
      static inline double absSqr (const SxComplex16& x) { return x.absSqr (); }
      /// Wrap absSqr for SxComplex16
      static inline float absSqr (const SxComplex8& x) { return x.absSqr (); }
      /// Wrap absSqr for real numbers
      template<class T> static inline T absSqr (const T x) { return x*x; }

      /// Power function
      template<StrideType Layout>
      static SxVector<double> pow (const SxVecRef<double, Layout> &vec, double a);

      // --- The next 15 lines were generated from snippets/SxVecCompute.h snippet UseIterator_reduce
      /// Add all elements
      template<class T, StrideType Layout>
      static inline T
      sum (const SxVecRef<T,Layout> &x);

      /// Multiply all elements
      template<class T, StrideType Layout>
      static inline T
      product (const SxVecRef<T,Layout> &x);

      /// Add abs. squares of all elements
      template<class T, StrideType Layout>
      static inline typename SxTypeMapper<T>::TReal
      normSqr (const SxVecRef<T,Layout> &x);
      // --- UseIterator_reduce

      /// Euclidean norm
      template<class T>
      static inline typename SxTypeMapper<T>::TReal
      norm (const SxVecRef<T,SubMatrix> &x)
      {
         return typename SxTypeMapper<T>::TReal(sqrt(x.normSqr ()));
      }

      /// Euclidean norm (strided vector) via BLAS
      template<class T>
      static inline typename SxTypeMapper<T>::TReal
      norm (const SxVecRef<T,Strided> &x)
      {
         SX_CHECK ((int)x.getRowStride () == x.getRowStride (),
                   x.getRowStride ());
         return norm2 (x.elements, x.getSize (), (int)x.getRowStride ());
      }

      /// Add abs. squares of all elements (strided vector) via BLAS
      template<class T>
      static inline typename SxTypeMapper<T>::TReal
      normSqr (const SxVecRef<T,Strided> &x)
      {
         typename SxTypeMapper<T>::TReal nrm = norm(x);
         return nrm * nrm;
      }

      /// Conjugate (only complex!)
      template<class T, StrideType Layout>
      static SxVector<SxComplex<T> >
      conj (const SxVecRef<SxComplex<T>, Layout> &vec);

      /// Conjugate (only complex!) - cannot be done in place...
      template<class T, StrideType Layout>
      static SxVector<SxComplex<T> >
      conjInPlace (SxVecRef<SxComplex<T>, Layout> &&vec)
      {
         return conj (vec);
      }

      /// Matrix-matrix multiplication
      template<class T, StrideType Layout1, StrideType Layout2>
      static inline SxVector<T> matmult (const SxVecRef<T,Layout1> & /* a */,
                                         const SxVecRef<T,Layout2> & /* b */)
      {
         // not implemented - probably inefficient
         SX_EXIT;
      }
};

/** \brief Container for Vector operations:
     One operand is packed, one is a contiguous vector
  */
template<>
class SxVecCompute<SimdExpand>
{
   public:
      // --- The next 17 lines were generated from snippets/SxVecCompute.h snippet Packed_matmult
      /// Matrix-matrix multiplication (first one is packed)
      template<class T, StrideType Layout2>
      static SxVector<T> matmult (const SxVecRef<T,PackedSymMatrix> &a,
                                  const SxVecRef<T,Layout2> &b);
      /// Matrix-matrix multiplication (second one is packed)
      template<class T, StrideType Layout2>
      static SxVector<T> matmult (const SxVecRef<T,Layout2> &a,
                                  const SxVecRef<T,PackedSymMatrix> &b);
      /// Matrix-matrix multiplication (first one is packed)
      template<class T, StrideType Layout2>
      static SxVector<T> matmult (const SxVecRef<T,BandedMatrix> &a,
                                  const SxVecRef<T,Layout2> &b);
      /// Matrix-matrix multiplication (second one is packed)
      template<class T, StrideType Layout2>
      static SxVector<T> matmult (const SxVecRef<T,Layout2> &a,
                                  const SxVecRef<T,BandedMatrix> &b);
      // --- Packed_matmult
};

/** \brief Container for Vector operations:
     One operand is packed, one is a strided vector
  */
template<>
class SxVecCompute<ItExpand>
{
   public:
      // --- The next 17 lines were generated from snippets/SxVecCompute.h snippet Packed_matmult
      /// Matrix-matrix multiplication (first one is packed)
      template<class T, StrideType Layout2>
      static SxVector<T> matmult (const SxVecRef<T,PackedSymMatrix> &a,
                                  const SxVecRef<T,Layout2> &b);
      /// Matrix-matrix multiplication (second one is packed)
      template<class T, StrideType Layout2>
      static SxVector<T> matmult (const SxVecRef<T,Layout2> &a,
                                  const SxVecRef<T,PackedSymMatrix> &b);
      /// Matrix-matrix multiplication (first one is packed)
      template<class T, StrideType Layout2>
      static SxVector<T> matmult (const SxVecRef<T,BandedMatrix> &a,
                                  const SxVecRef<T,Layout2> &b);
      /// Matrix-matrix multiplication (second one is packed)
      template<class T, StrideType Layout2>
      static SxVector<T> matmult (const SxVecRef<T,Layout2> &a,
                                  const SxVecRef<T,BandedMatrix> &b);
      // --- Packed_matmult
};

template<>
class SxVecCompute<Simd>::Divide<false> {
   public:
      // vector / scalar (vector contiguous, same non-integer type)
      template<class T, StrideType Layout, class T2>
      static typename SxVecResult<decltype(T(0)+T2(0)), Layout>::VecType
      divide (const SxVecRef<T, Layout> &x, const T2 &y)
      {
         // multiply with inverse
         SX_CHECK_DIV(y);
         return multiply (x, typename SxTypeMapper<T2>::TReal(1) / y);
      }

      /// divide vector and scalar (in place)
      template<class T, StrideType Layout, class T2>
      static void
      divideInPlace (SxVecRef<T, Layout> &x, const T2 &y)
      {
         // multiply with inverse
         SX_CHECK_DIV(y);
         multiplyInPlace (x, typename SxTypeMapper<T2>::TReal(1) / y);
      }

};

template<>
class SxVecCompute<Simd>::Divide<true> {
   public:
      /// divide vector and scalar
      template<class T, StrideType Layout, class T2>
      static typename SxVecResult<decltype(T(0)+T2(0)), Layout>::VecType
      divide (const SxVecRef<T, Layout> &x, const T2 &y);

      /// divide vector and scalar (in place)
      template<class T, StrideType Layout, class T2>
      static void
      divideInPlace (SxVecRef<T, Layout> &x, const T2 &y);
};

template<>
class SxVecCompute<UseIterator>::Divide<false> {
   public:
      // vector / scalar (any vector, non-integer divisor)
      template<class T, StrideType Layout, class T2>
      static SxVector<decltype(T(0)+T2(0))>
      divide (const SxVecRef<T, Layout> &x, const T2 &y)
      {
         // multiply with inverse
         SX_CHECK_DIV(y);
         return multiply (x, typename SxTypeMapper<T2>::TReal(1) / y);
      }

      /// divide vector and scalar (in place)
      template<class T, StrideType Layout, class T2>
      static void
      divideInPlace (SxVecRef<T, Layout> &x, const T2 &y)
      {
         // multiply with inverse
         SX_CHECK_DIV(y);
         multiplyInPlace (x, typename SxTypeMapper<T2>::TReal(1) / y);
      }

};

template<>
class SxVecCompute<UseIterator>::Divide<true> {
   public:
      /// divide vector and scalar
      template<class T, StrideType Layout, class T2>
      static SxVector<decltype(T(0)+T2(0))>
      divide (const SxVecRef<T, Layout> &x, const T2 &y);

      /// divide vector and scalar (in place)
      template<class T, StrideType Layout, class T2>
      static void
      divideInPlace (SxVecRef<T, Layout> &x, const T2 &y);
};

#endif /* _SX_VEC_COMPUTE_H_ */
