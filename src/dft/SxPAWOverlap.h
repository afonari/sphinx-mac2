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

#ifndef _SX_PAW_OVERLAP_H_
#define _SX_PAW_OVERLAP_H_

#include <SxDFT.h>
#include <SxOverlap.h>
#include <SxPartialWaveBasis.h>

/** \brief Projector augmented wave overlap operator

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxPAWOverlap : public SxOverlapBase
{
   protected:
      /// PW projectors
      SxConstPtr<SxPartialWaveBasis> pBasis;
   public:
      /** \brief Constructor
        */
      SxPAWOverlap (const SxConstPtr<SxPartialWaveBasis> &pIn)
         : pBasis(pIn)
      { /* empty */ }

      /// Norm square
      virtual double normSqr (const SxVecRef<SxComplex16> &psi) const;

      /// Scalar product
      virtual SxComplex16 dot (const SxVecRef<SxComplex16> &x,
                               const SxVecRef<SxComplex16> &y) const;

   protected:
      /// Partial wave contribution to scalar product
      SxComplex16 dotPartial (const SxVecRef<SxComplex16> &px,
                              const SxVecRef<SxComplex16> &py) const;

   public:
      /// Apply overlap operator
      virtual SxVecRef<SxComplex16>
      apply (const SxVecRef<SxComplex16> &psi) const;

      /// Set states x orthogonal to states y
      virtual void setOrthogonal (SxVecRef<SxComplex16> *xPtr,
                                  const SxVecRef<SxComplex16> &y) const;

      /// Get PAW correction to overlap matrix from projections <p|psi>
      SxVector<PrecCoeffG>
      getDeltaS (const SxVecRef<PrecCoeffG> &px,
                 const SxVecRef<PrecCoeffG> &py) const;

      /// Get PAW overlap matrix elements
      virtual SxVector<SxComplex16>
      getMatrix (const SxVecRef<SxComplex16> &x,
                 const SxVecRef<SxComplex16> &y) const;

      /// Orthonormalize states x
      virtual void orthonormalize (SxVecRef<SxComplex16> *xPtr,
                                   const SxOrthoMethod how=GramSchmidt) const;
      /// Normalize states x
      virtual void normalize (SxVecRef<SxComplex16> *xPtr) const;

};

#endif /* _SX_PAW_OVERLAP_H_ */

