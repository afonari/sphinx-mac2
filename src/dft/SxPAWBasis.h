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

#ifndef _SX_PAW_BASIS_H_
#define _SX_PAW_BASIS_H_

#include <SxDFT.h>
#include <SxGBasis.h>
#include <SxBasisAug.h>
#include <SxPartialWaveBasis.h>

/** \brief Full PAW basis

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxPAWBasis
: public SxBasisAug
{
   public:
      /// Pointer to G+k basis
      SxConstPtr<SxGBasis> gBasis;

      /// Pointer to partial wave basis
      SxConstPtr<SxPartialWaveBasis> pBasis;

      /// Constructor
      SxPAWBasis (const SxConstPtr<SxGBasis> &gIn,
                  const SxConstPtr<SxPartialWaveBasis> &pIn)
         : gBasis (gIn),
           pBasis (pIn)
      {
         SX_CHECK (gIn);
         SX_CHECK (pIn);
      }  

      // --- SxBasis stuff
      typedef SxComplex16 CoeffType;

      virtual ~SxPAWBasis () {
         // empty
      }

      // --- overloaded virtual SxBasis functions ----------
      virtual ssize_t getNElements () const
      {
         return gBasis->getNElements () + pBasis->getNElements ();
      }

      /** \brief Print debug information about the basis */
      virtual void print () const;

      /** Very simple description of basis a la "|R>" or "|G>"
        */
      virtual SxString getType () const  {
         return "|G+k,p>";
      }

      virtual 
      SxComplex16 scalarProduct (const SxVecRef<SxComplex16> &x,
                                 const SxVecRef<SxComplex16> &y) const;

      REGISTER_PROJECTOR (SxGBasis, toPWBasis);
      REGISTER_PROJECTOR (SxRBasis, toRBasis);
      REGISTER_PROJECTOR (SxAOBasis, toAO);

      //PROJECT_TOEXTRA_VIA_FROM (SxBasisAug, float);
      //PROJECT_TOEXTRA_VIA_FROM (SxBasisAug, double);
      //PROJECT_TOEXTRA_VIA_FROM (SxBasisAug, SxComplex8);
      PROJECT_TOEXTRA_VIA_FROM (SxBasisAug, SxComplex16);

      /** \brief Project to G-basis */
      SxVector<PrecCoeffG> toPWBasis (const SxGBasis *,
                                      const SxVecRef<PrecCoeffG> &) const;
      /** \brief Project to R-basis */
      SxVector<PrecCoeffR> toRBasis (const SxRBasis *,
                                     const SxVecRef<PrecCoeffG> &) const;
      /** \brief Project to R-basis */
      SxVector<PrecCoeffG> toAO     (const SxAOBasis *,
                                     const SxVecRef<PrecCoeffG> &) const;
      /** \brief Projector from \b PAW space to some anonymous basis

          This is the versatile interface: we try to dynamically determine
          the basis to project to.

          \sa \ref page_dirac */
      virtual SxVecRef<SxComplex16>
      projectTo (const SxBasis *, const PsiRef &) const;

      /// Project from SxGBasis
      virtual PsiRef projectFrom (const SxGBasis *gBasis,
                                  const PsiRef   &psi    ) const;

      /// Identity
      virtual PsiRef projectFrom (const SxPAWBasis *me,
                                  const PsiRef     &psi    ) const
      {
         SX_CHECK (me == this);
         return const_cast<PsiRef&>(psi);
      }

      /// Reference to G-basis related coefficients
      virtual SxVecRef<SxComplex16, SubMatrix>
      toBasis (const SxGBasis *gIn, const SxVecRef<SxComplex16> &psi) const;

      /// Reference to partial wave basis related coefficients
      virtual SxVecRef<SxComplex16, SubMatrix>
      toBasis (const SxPartialWaveBasis *pIn,
               const SxVecRef<SxComplex16> &psi) const;
};

#endif /* _SX_PAW_BASIS_H_ */
