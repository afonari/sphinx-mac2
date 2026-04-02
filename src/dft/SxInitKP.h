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

#ifndef _SX_INIT_KP_H_
#define _SX_INIT_KP_H_

#include <SxDFT.h>
#include <SxHamInit.h>

/** \brief Initializer class for k.p Hamiltonians

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxInitKP : public SxHamInit
{
   public:
      enum KpType { KdotPgeneral, StrainField };
   protected:
      KpType type;
      const SxRBasis *rPtr;
      int nComp;
   public:
      SxInitKP (KpType typeIn, const SxRBasis* rPtrIn, const SxSymbolTable *top = NULL)
         : type(typeIn), rPtr (rPtrIn)
      {
         if (type == KdotPgeneral)  {
            SX_CHECK (top);
            nComp = SxHamInit::getNComp (top);
         } else if (type == StrainField)  {
            nComp = 3;
         }
      }

      /// Randomize density
      virtual void randomRho () {}
      /// Atomic density
      virtual void atomicRho (const SxSymbolTable *,
                              const SxAtomicStructure &) {}
      /// Density from file
      virtual void readRho (const SxString &) {}
      /// Density from wave functions
      virtual void rhoFromWaves (const SxPtr<SxPWSet> &,
                                 const Focc &) {}
      /// Add extra charge to density
      virtual void addCharge (const SxMeshR &) {}
      /// Initialize waves
      virtual SxPtr<SxPWSet> setupWaves (const SxSymbolTable *wavesGroup,
                                         const SxPtr<SxGkBasis> &gkPtr,
                                         int nStates, int nSpin)
      {
         return setupPW (wavesGroup, gkPtr, nStates, nSpin, nComp);
      }
      /// Randomize waves
      virtual void randomize (const SxPtr<SxPWSet> &wavesPtr);
      /// Set up the Hamiltonian
      virtual SxPtr<SxHamiltonian> setupHam (const SxPtr<SxPWSet> &,
                                          const SxPtr<SxSpeciesData> &,
                                          const SxAtomicStructure &);

      /// Virtual destructor
      virtual ~SxInitKP ();
};

#endif /* _SX_INIT_KP_H_ */
