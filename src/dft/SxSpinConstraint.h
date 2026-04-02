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

#ifndef _SX_SPINCONSTRAINT_H_
#define _SX_SPINCONSTRAINT_H_

#include <SxSpinConstraint.h>
#include <SxFermi.h>
#include <SxPAWSet.h>
#include <SxPAWRho.h>
#include <SxRho.h>
#include<SxHamiltonian.h>
#include <SxPAWHamiltonian.h>
#include <SxPartialWaveBasis.h>
#include <SxAtomicStructure.h>
#include <SxRBasis.h>
#include <SxGrid.h>


class SX_EXPORT_DFT SxSpinConstraint
{
   public:

      SxPtr<SxPAWPot> pawPotPtr;
      SxVector<double> nuA;
      SxVector<double> targetSpin;
      /// whether or not a specific atom is spin-constrained (:iTlAtom)
      SxArray<bool> constrainedAtom;
      SxSpinConstraint( const SxPtr<SxPAWPot> &pawPotPtrIn);
      /// Read constraints from symbol table
      void read (const SxSymbolTable *topLvl,
                 const SxAtomicStructure &structure);
      
      void setNuA (double value);
      /** \brief Set the constraints
          @param targetSpinIn - the desired target spins
          @param str - structure for symmetry checking
        */
      void setConstraints ( const SxVecRef<double> &targetSpinIn,
                           const SxAtomicStructure &str);
      /// Check compatibility of spin constraint with symmetry
      void checkSym (const SxAtomicStructure &str) const;
      double calcKappaOpt (const SxVector<double> &MSpinIn, const SxVector<double> &MSpinPlusIn, double kappa);
      /// Compute new nu Lagrangian parameters for given waves
      SxVector<double> computeNu (const SxPtr<SxPAWSet> &wavesPtr, SxFermi &fermi, const Real8 ekt, double epsilon = 1e-20);

   private:
      /// Get partial volume in state-space for one atom
      SxVector<SxComplex16>
      getOmegaNN (const PsiRef &psi, ssize_t ipBasis, int is) const;
      /// Get partial volume for each state for all atoms from PAW projections
      SxArray<SxVector<SxComplex16> >
      getOmegaN (const SxVecRef<SxComplex16> &proj) const;
      /// Get spins for current nuA with subspace diagonalization
      SxVector<double> getSpinsDiag ( const Eps &eps0In, SxFermi &fermiIn, const double &ektIn , SxPtr<SxPAWSet> wavesPtr = SxPtr<SxPAWSet>() );
      /// Pointer to the currently used waves
      const SxPAWSet *wvPtr;

};
#endif  /*_SX_SPINCONSTRAINT_H_*/
