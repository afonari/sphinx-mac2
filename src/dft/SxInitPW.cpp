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

#include <SxInitPW.h>
#include <SxPWHamiltonian.h>

SxInitPW::SxInitPW (const SxRBasis &R, int nSpin, double nElectrons,
                    const SxPtr<SxPseudoPot> &potPtrIn)
   : potPtr(potPtrIn),
     rho(R, nSpin, nElectrons)
{
   // empty
}

void SxInitPW::randomRho ()
{
   rho.randomize ();
}

void SxInitPW::atomicRho (const SxSymbolTable *rhoGroup,
                          const SxAtomicStructure &structure)
{
   SxPseudoPot &psPot = dynamic_cast<SxPseudoPot&> (*potPtr);
   int nSpin = rho.getNSpin ();
   
   int iSpecies, nSpecies = psPot.nSpecies;
   int iSpin;
   int l, lMax;
   double spinMoment, s = 1., sumFocc;
   // --- create local LCAO occupations
   SxArray<SxArray<SxVector<double> > >  foccAtomicRho(nSpecies);
   // --- resize local LCAO occupations
   for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      lMax = psPot.lMax(iSpecies);
      foccAtomicRho(iSpecies).resize(lMax+1);
      for (l=0; l <= lMax; l++)  {
         foccAtomicRho(iSpecies)(l).resize(nSpin);
      }
   }
   if (rhoGroup && rhoGroup->contains("spinMoment"))  {
      double nElectrons = rho.nElectrons;
      spinMoment = rhoGroup->get("spinMoment")->toReal();
      s          = (nElectrons - spinMoment)/(nElectrons + spinMoment);
   }
   
   for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      lMax = psPot.lMax(iSpecies);
      if (fabs(s-1.) > 1e-12 && (nSpin == 2))  {
         for (l=0; l <= lMax; l++)  {
            sumFocc = psPot.foccAtomicRho(iSpecies)(l).sum();
            foccAtomicRho(iSpecies)(l)( SPIN_UP ) = sumFocc * 1./(1.+s);
            foccAtomicRho(iSpecies)(l)(SPIN_DOWN) = sumFocc *  s/(1.+s);
         }      
      }  else  {
         if (psPot.foccAtomicRho(iSpecies)(0).getSize () == nSpin)  {
            foccAtomicRho(iSpecies) = psPot.foccAtomicRho(iSpecies);
         }  else  {
            for (l=0; l <= lMax; l++)  {
               sumFocc = psPot.foccAtomicRho(iSpecies)(l).sum();
               for (iSpin=0; iSpin < nSpin; iSpin++)  {
                  foccAtomicRho(iSpecies)(l)(iSpin) = 
                     sumFocc / (Real8)nSpin;
               }
            }
         }
      }
   }
   
   // --- copy local occupations to PP object 
   psPot.foccAtomicRho = foccAtomicRho;
   
   // --- calculate atomic charge density
   rho.atomicChargeDensity (structure, psPot);
   
}

void SxInitPW::readRho (const SxString &file)
{
   bool Hirshfeld = false; // to be tested
   rho.readRho (file);
   const SxGBasis &G = rho.rBasisPtr->getGBasis ();
   const SxAtomicStructure &structure = G.getTau ();
   SxAtomicStructure tauFile;
   try {
      SxBinIO io(file, SxBinIO::BINARY_READ_ONLY);
      if (io.contains ("tau")) tauFile.read (io);
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   if (tauFile.getNAtoms () > 0)  {
      if (tauFile.getNSpecies () == structure.getNSpecies ()
          && tauFile.getNAtoms () == tauFile.getNAtoms ()
          && (  tauFile  .atomInfo->nAtoms
              - structure.atomInfo->nAtoms).sqr ().sum () == 0)
      {
         tauFile.replaceInfo (structure.atomInfo);
         if ((tauFile - structure).absSqr ().sum ()
               > 1e-6 * structure.getNAtoms ())
         {
            // --- structure has changed: move atomic density
            cout << SX_SEPARATOR;
            cout << "| WARNING: initial density from displaced atoms."
                 << endl << SX_SEPARATOR;
            const_cast<SxGBasis&>(G).changeTau (tauFile);
            const SxRBasis &R = *rho.rBasisPtr;
            // --- get shape of ideal atomic densities in G space
            //     and overlap of all those densities
            SxArray<SxMeshG> atomRho(potPtr->getNSpecies ());
            SxMeshG allAtom(0., G);
            SX_LOOP(is)  {
               atomRho(is) = potPtr->getAtomRhoG (G, (int)is);
               allAtom += atomRho(is) * G.structureFactors(is);
            }
            SxVector<double> allAtomR = R.symmetrize (R | allAtom);
            // switch to the deformation density: subtract ideal atoms
            SX_LOOP(iSpin)
               rho(iSpin).plus_assign_ax(-1./rho.getNSpin (), allAtomR);

            if (Hirshfeld) rho.displaceHirshfeld (structure, atomRho, &allAtomR);

            // --- now change the structure
            const_cast<SxGBasis&>(G).changeTau (structure);

            // switch back from deformation density to full density
            // by adding the atomic contributions
            allAtom.set (0.);
            SX_LOOP(is) allAtom += atomRho(is) * G.structureFactors(is);
            allAtomR = R | allAtom;
            SX_LOOP(iSpin)
               rho(iSpin).plus_assign_ax(+1./rho.getNSpin (), allAtomR);
         }
      } else {
         cout << "WARNING: initial density from different structure"
            << endl;
      }
   }
}

void SxInitPW::addCharge (const SxMeshR &chargeR)
{
   int nSpin = rho.getNSpin ();
   double totalExtra = chargeR.sum () * rho.rBasisPtr->dOmega;
   rho.nElectrons += totalExtra;
   rho.normalizeRho ();
   rho.nElectrons -= totalExtra;
   for (int iSpin = 0.; iSpin < nSpin; ++iSpin)
      rho(iSpin).plus_assign_ax(-1./nSpin, chargeR);
}


SxPtr<SxPWSet> SxInitPW::setupWaves (const SxSymbolTable *wavesGroup,
                                     const SxPtr<SxGkBasis> &gkPtr,
                                     int nStates, int nSpin)
{
   return setupPW (wavesGroup, gkPtr, nStates, nSpin);
}

void SxInitPW::randomize (const SxPtr<SxPWSet> &wavesPtr)
{
   SxPtr<SxPW> (wavesPtr)->randomize ();
}

void SxInitPW::rhoFromWaves (const SxPtr<SxPWSet> &wavesPtr, const Focc &focc)
{
   rho.computeRho (focc, *SxPtr<SxPW> (wavesPtr));
}

SxPtr<SxHamiltonian>
SxInitPW::setupHam (const SxPtr<SxPWSet> &wavesPtr,
                    const SxPtr<SxSpeciesData> &,
                    const SxAtomicStructure &structure)
{
   SX_CHECK (dynamic_cast<SxPseudoPot*> (potPtr.getPtr ()));
   SxPseudoPot &psPot = dynamic_cast<SxPseudoPot&> (*potPtr);
   SxPW &waves = *SxPtr<SxPW>(wavesPtr);
   return SxPtr<SxPWHamiltonian>::create (waves, rho.rhoR, psPot, structure);
}
