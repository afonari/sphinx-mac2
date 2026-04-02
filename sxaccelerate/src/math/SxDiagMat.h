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

#ifndef _SX_DIAG_MAT_H_
#define _SX_DIAG_MAT_H_

#include <SxDiagMat.h>

/** \brief Diagonal matrix

    \author C. Freysoldt, freysoldt@mpie.de */
template <class T = double>
class SX_EXPORT_MATH SxDiagMat
{
   public:
      /// Diagonal elements
      SxVector<T> diagVals;

      /// Constructor
      template<class T2, StrideType Layout>
      SxDiagMat (SxVector<T> &&in)
      : diagVals (std::move (in))
      {
         SX_CHECK (diagVals.getNCols () == 1, diagVals.getNCols ());
      }

      /// Constructor
      template<class T2, StrideType Layout>
      SxDiagMat (const SxVecRef<T2, Layout> &in)
      : diagVals (in)
      {
         SX_CHECK (diagVals.getNCols () == 1, diagVals.getNCols ());
      }

      /// Empty constructor
      SxDiagMat () {}

      /// Matrix multiply diag ^ A
      template<class T2, StrideType Layout>
      inline SxVector<decltype(T(0) * T2(0))>
      operator^ (const SxVecRef<T2, Layout> &A) const;

      /// Matrix multiply diag ^ A (A is dying)
      template<class T2>
      inline SxVector<decltype(T(0) * T2(0))>
      operator^ (SxVecRef<T2> &&A) const
      {
         /* constexpr */ if (std::is_same<T2, decltype(T(0) * T2(0))>::value) {
            if (A.isOnlyRef ())   {
               SX_CHECK (A.getNRows () == diagVals.getSize (),
                         A.getNRows (), diagVals.getSize ());
                  for (ssize_t i = 0; i < A.getNCols (); i++)
                     A.colRef(i) *= diagVals;
                  return std::move (A);
            }
         }
         return operator^ ((const SxVecRef<T2>&)A);
      }
};

template<class T>
template<class T2, StrideType Layout>
inline SxVector<decltype(T(0) * T2(0))>
SxDiagMat<T>::operator^ (const SxVecRef<T2, Layout> &A) const
{
   SX_CHECK (!isPacked(Layout)); // not yet implemented
   SX_CHECK (A.getNRows () == diagVals.getSize (),
             A.getNRows (), diagVals.getSize ());
   ssize_t N = A.getNRows (), M = A.getNCols ();
   SxVector<decltype(T(0) * T2(0))> res(N, M);
   res.auxData = A.auxData;
   if (Layout == Compact || Layout == SubMatrix)  {
      SxVecCompute<Simd>::multiplyMx1 (res.elements, N, M,
                                       A.elements, A.getColStride (),
                                       diagVals.elements);
   } else {
      SxVecConstIt<T2, Layout> colBegin = A.begin ();
      for (ssize_t i = 0; i < M; ++i)  {
         SxVecCompute<SimdIt>::multiply (colBegin, diagVals.elements, N,
                                         res.elements + i * N);
         colBegin += N;
      }
   }
   return res;
}

template<class T, class T2, StrideType Layout>
inline SxVector<decltype(T(0) * T2(0))>
operator^(const SxVecRef<T2, Layout> &A, const SxDiagMat<T> &D)
{
   SX_CHECK (!isPacked(Layout)); // not yet implemented
   SX_CHECK (A.getNCols () == D.diagVals.getSize (),
             A.getNCols (), D.diagVals.getSize ());
   SxVector<decltype(T(0) * T2(0))> res(A.getNRows (), A.getNCols ());
   res.auxData = A.auxData;
   if (Layout == Compact || Layout == SubMatrix)  {
      for (ssize_t i = 0; i < A.getNCols (); ++i)  {
         SxVecCompute<Simd>::multiply (res.elements + i * res.getColStride (),
                                       res.getNRows (),
                                       A.elements + i * A.getColStride(),
                                       D.diagVals(i));
      }
   } else {
      SxVecConstIt<T2, Layout> aIt = A.begin ();
      SxVecIt<decltype(T(0) * T2(0)),Compact> resIt = res.begin ();
      for (ssize_t i = 0; i < A.getNCols (); ++i)  {
         T d = D.diagVals(i);
         for (ssize_t j = 0; j < A.getNRows (); ++j)
            *resIt++ = *aIt++ * d;
      }
   }
   return res;
}

template<class T, class T2>
inline SxVector<decltype(T(0) * T2(0))>
operator^(SxVecRef<T2> &&A, const SxDiagMat<T> &D)
{
   /* constexpr */ if (std::is_same<T2, decltype(T(0) * T2(0))>::value)  {
      if (A.isOnlyRef ())  {
         SX_CHECK (A.getNCols () == D.diagVals.getSize (),
                   A.getNCols (), D.diagVals.getSize ());
         for (ssize_t i = 0; i < A.getNCols (); ++i)
            A.colRef(i) *= D.diagVals(i);
         return std::move (A);
      }
   }
   return static_cast<const SxVecRef<T2>&>(A) ^ D;
}

#endif /* _SX_DIAG_MAT_H_ */
