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

#ifndef _SX_PARTIAL_WAVE_BASIS_H_
#define _SX_PARTIAL_WAVE_BASIS_H_

#include <SxDFT.h>
#include <SxBasis.h>
#include <SxPAWPot.h>
#include <SxAOBasis.h>
#include <SxBasisAug.h>

class SxPAWBasis;

/** \brief Partial wave basis (PAW method)

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxPartialWaveBasis : public SxBasisAug
{
   public:
      typedef SxComplex16 CoeffType;
   protected:
      /// Pointer to PAW potentials (access via getPotPtr or getPot)
      SxConstPtr<SxPAWPot> potPtr;

      /// Number of elements
      ssize_t nElements;

   public:
      /// Get pointer to PAW potential
      SxConstPtr<SxPAWPot> getPotPtr () const
      {
         return potPtr;
      }

      /// Return reference to PAW potential
      const SxPAWPot& getPot () const
      {
         SX_CHECK (potPtr);
         return *potPtr;
      }

      /// Pointer to projector basis
      SxConstPtr<SxAOBasis> projectors;

      /// Constructor
      SxPartialWaveBasis (const SxConstPtr<SxPAWPot> &pot,
                          const SxAtomicStructure &str);

      virtual ~SxPartialWaveBasis ()
      {
         // empty
      }

      // --- overloaded virtual SxBasis functions ----------

      /// \brief Number of sampling points of the corresponding Dirac vector.
      virtual ssize_t getNElements () const
      {
         return nElements;
      }

      /** \brief Print debug information about the basis */
      virtual void print () const;

      /** Very simple description of basis a la "|R>" or "|G>"
        */
      virtual SxString getType () const  {
         return "|p>";
      }

      /** \brief Project to G-basis */
      PsiG toPWBasis (const SxGBasis *G, const PsiRef &pCoeff) const
      {
         SX_CHECK (projectors);
         return projectors->toPWBasis (G, pCoeff);
      }

      REGISTER_PROJECTOR (SxGBasis, toPWBasis);

      /// Project from SxPAWBasis
      virtual PsiRef projectFrom (const SxPAWBasis *pawIn, const PsiRef &psi) const;
      /// Project from SxGBasis
      virtual PsiRef projectFrom (const SxGBasis *gBasis, const PsiRef &psi) const;

      /** \brief Create a projector basis
        @param gk G+k basis
        @param pawPot partial wave basis
        */
      static SxPtr<SxAOBasis> createProjBasis (const SxGkBasis &gk,
                                               const SxPAWPot  &pawPot);

      /** \brief Create a projector basis
        @param gk G+k basis
        */
      inline void createProjBasis (const SxGkBasis &gk)
      {
         SX_CHECK (potPtr);
         projectors = createProjBasis (gk, *potPtr);
      }
};

#endif /* _SX_PARTIAL_WAVE_BASIS_H_ */
