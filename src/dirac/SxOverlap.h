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

#ifndef _SX_OVERLAP_H_
#define _SX_OVERLAP_H_

#include <SxDiracLib.h>
#include <SxVector.h>
#include <SxPtr.h>
#include <SxTimer.h>

/** Determines what method in orthogonalization/orthonormalization
    will be applied.
 */
enum SxOrthoMethod { GramSchmidt, Loewdin};

/** \brief Generic overlap operator interface

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DIRAC SxOverlapBase
{
   public:
      /// Virtual destructor
      virtual ~SxOverlapBase () = default;

      /// Norm square
      virtual double normSqr (const SxVecRef<SxComplex16> &psi) const = 0;

      /// Scalar product
      virtual SxComplex16 dot (const SxVecRef<SxComplex16> &x,
                               const SxVecRef<SxComplex16> &y) const = 0;

      /// Apply overlap operator
      virtual SxVecRef<SxComplex16>
      apply (const SxVecRef<SxComplex16> &psi) const = 0;

      /// Set states x orthogonal to states y
      virtual void setOrthogonal (SxVecRef<SxComplex16> *xPtr,
                                  const SxVecRef<SxComplex16> &y) const = 0;

      /// Orthonormalize states x
      virtual void orthonormalize (SxVecRef<SxComplex16> *xPtr,
                                   const SxOrthoMethod how = GramSchmidt) const = 0;

      /// Normalize states x
      virtual void normalize (SxVecRef<SxComplex16> *xPtr) const;

      /// Get overlap matrix
      virtual SxVector<SxComplex16>
      getMatrix (const SxVecRef<SxComplex16> &x,
                 const SxVecRef<SxComplex16> &y) const
      {
         SX_CHECK (x.getSize () > 0);
         SX_CHECK (y.getSize () > 0);
         SX_CHECK (x.getBasisPtr () == y.getBasisPtr ());
         if (x.getNCols () > y.getNCols ())
            return x.overlap (apply (y));
         else
            return apply (x).overlap (y);
      }

};
/** \brief Generic overlap operator

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DIRAC SxOverlap
{
   protected:
      SxConstPtr<SxOverlapBase> S;
   public:
      SxOverlap (const SxConstPtr<SxOverlapBase> &in)
         : S(in)
      { /* empty */ }

      /// Set states x orthogonal to states y
      void setOrthogonal (SxVecRef<SxComplex16> *xPtr,
                          const SxVecRef<SxComplex16> &y) const
      {
         S->setOrthogonal (xPtr, y);
      }

      /// Orthonormalize states x
      void orthonormalize (SxVecRef<SxComplex16> *xPtr,
                   const SxOrthoMethod how = GramSchmidt) const
      {
         S->orthonormalize (xPtr, how);
      }

      /// Normalize states x
      void normalize (SxVecRef<SxComplex16> *xPtr) const
      {
         S->normalize (xPtr);
      }

      SxVector<SxComplex16>
      getMatrix (const SxVecRef<SxComplex16> &x,
                 const SxVecRef<SxComplex16> &y) const
      {
         return S->getMatrix (x, y);
      }

   protected:
      /// Overlap operator expression
      class SxOvlpExpr  {
         protected:
            /// Type of expression
            const enum { SPsi, PsiS } type;

            /// Overlap operator reference
            const SxOverlapBase &S;

            /// psi reference
            const SxVecRef<SxComplex16> &psi;
         public:
            SxOvlpExpr (const SxOverlapBase &SIn,
                        const SxVecRef<SxComplex16> &psiIn)
               : type(SPsi), S(SIn), psi(psiIn)
            { /* empty */ }

            SxOvlpExpr (const SxVecRef<SxComplex16> &psiIn,
                        const SxOverlapBase &SIn)
               : type(PsiS), S(SIn), psi(psiIn)
            { /* empty */ }

            /// Perform S | psi
            operator SxVecRef<SxComplex16> () const
            {
               SX_CHECK (type == SPsi);
               return S.apply (psi);
            }

            /// Perform (psi | S | psi2)
            SxComplex16 operator| (const SxVecRef<SxComplex16> &psi2)
            {
               SX_CHECK (type == PsiS);
               if (&psi == &psi2)
                  return S.normSqr (psi);
               else
                  return S.dot (psi, psi2);
            }
      };

   public:
      /// Get S|psi expression
      SxOvlpExpr operator| (const SxVecRef<SxComplex16> &psi)
      {
         return SxOvlpExpr (*S, psi);
      }

      /// Apply S
      SxVecRef<SxComplex16> operator* (const SxVecRef<SxComplex16> &psi) const
      {
         return S->apply (psi);
      }

      friend SxOvlpExpr operator| (const SxVecRef<SxComplex16> &,
                                   const SxOverlap &S);
};

/// Get psi|S expression
inline
SxOverlap::SxOvlpExpr operator| (const SxVecRef<SxComplex16> &psi,
                                 const SxOverlap &S)
{
   return SxOverlap::SxOvlpExpr (psi, *S.S);
}

namespace Timer {
   enum OrthoTimer {
      ortho,
      loewdinS
   };
}

SX_REGISTER_TIMERS (Timer::OrthoTimer)
{
   using namespace Timer;
   regTimer (ortho,      "Orthogonalization");
   regTimer (loewdinS,   "Loewdin S");
}


#endif /* _SX_OVERLAP_H_ */
