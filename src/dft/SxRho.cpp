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

#include <SxRho.h>
#include <SxGkBasis.h>
#include <SxRadBasis.h>
#include <SxRadialBasis.h>
#include <SxQuantumNumbers.h>
#include <SxConstants.h>
#include <SxYlm.h>
#include <SxLoopMPI.h>  // LoopMPI
#include <SxPAWSet.h>

SxRho::SxRho ()
{
   rBasisPtr = NULL;
   nElectrons = 0.;
}

SxRho::SxRho (const RhoR &rhoRIn)
{
   rhoR       = rhoRIn;
   rBasisPtr  = dynamic_cast<const SxRBasis *>(rhoRIn(0).getBasisPtr());
   SX_CHECK (rBasisPtr);
   nElectrons = getNorm ();
}
SxRho::SxRho (const SxMeshR &rhoIn)
{
   rhoR.resize (1);
   rhoR(0) = rhoIn;
   rBasisPtr  = dynamic_cast<const SxRBasis *>(rhoIn.getBasisPtr());
   SX_CHECK (rBasisPtr);
   nElectrons = getNorm ();
}

SxRho::SxRho (const SxRBasis &R, int nSpin, Real8 nElectrons_)
{
   SX_CHECK (nSpin == 1 || nSpin == 2, nSpin);
   rBasisPtr = &R;
   nElectrons = nElectrons_;
   rhoR.resize (nSpin);
   for (int iSpin=0; iSpin < nSpin; ++iSpin) 
      rhoR(iSpin).resize (R.getMeshSize());
}

SxRho::SxRho (SxBinIO &io, const SxRBasis* rBasisPtrIn)
   : rBasisPtr(rBasisPtrIn)
{
   SX_CHECK (io.mode == SxBinIO::BINARY_READ_ONLY);
   int nSpin;
   try  {
      nSpin = io.getDimension ("nMeshes");
      SX_CHECK(nSpin > 0, nSpin);

      if (rBasisPtr)  {
         rhoR.resize (nSpin);
         readRho (io);
      } else {

         // --- read meshes
         SxMatrix3<double> cell;
         SxVector3<int>    dim;
         rhoR = io.readMesh (&cell, &dim);

         // --- determine nElectrons
         nElectrons = 0.;
         for (int iSpin = 0; iSpin < nSpin; iSpin++)
            nElectrons += rhoR(iSpin).sum ();
         
         double dOmega = fabs(cell.determinant ()) / double(dim.product ());

         nElectrons *= dOmega;

      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

}

void SxRho::operator= (const SxRho &in)
{
   rBasisPtr  = in.rBasisPtr;
   rhoR       = in.rhoR;
   nElectrons = in.nElectrons;
}

SxMeshR &SxRho::operator() (int iSpin)
{
   return rhoR(iSpin);
}

const SxMeshR &SxRho::operator() (int iSpin) const
{
   return rhoR(iSpin);
}

const SxBasis *SxRho::getBasisPtr () const
{
   return rhoR(0).getBasisPtr();
}

const SxRBasis *SxRho::getRBasisPtr () const
{
   return rBasisPtr;
}


RhoR &SxRho::computeRho (const Focc &focc, const SxPWSet &psiSet)
{
   SX_CHECK (rBasisPtr);
   SX_CLOCK (Timer::Rho);

   // --- basis sets
   const SxGkBasis &gk = psiSet.getGkBasis ();  // |G+k>
   const SxRBasis  &R  = *rBasisPtr;            // |R>

   int i, nStates, iSpin, ik;
   int nSpin = psiSet.getNSpin ();

   // --- initialize rhoR
   for (iSpin=0; iSpin < nSpin; iSpin++)  {
      rhoR(iSpin).resize (R.getMeshSize());
      rhoR(iSpin).set (0.);
      rhoR(iSpin).setBasis (&R);
   }
   bool blockStates =   dynamic_cast<const SxPW*>(&psiSet)
                     || dynamic_cast<const SxPAWSet*>(&psiSet);

   // --- compute density
   double f;
   SX_NEW_LOOP (psiSet);
   for (ik=0; ik < gk.nk; ik++)  {
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         SX_MPI_LEVEL("waves-k");
         if (SxLoopMPI::myWork(ik)) {  // SxParallelHierarchy
            nStates = minimum (psiSet.getNStates (ik), focc.getNStates(ik));
            if (blockStates && gk(ik).hasMixedFFT ())  {
               SX_CHECK (psiSet.getNStates (ik) <= focc.getNStates(ik),
                         psiSet.getNStates (ik), focc.getNStates(ik));

               gk(ik).addToRho (psiSet(iSpin, ik),
                                gk.weights(ik) * focc(iSpin,ik),
                                 &rhoR(iSpin));
            } else {
               for (i=0; i < nStates; i++)  {
                  f = gk.weights(ik) * focc(i,iSpin,ik);
                  if (fabs(f) > 1e-12)  {
                     /*
                  rhoR(iSpin) += (gk.weights(ik) * focc(i,iSpin,ik))
                      * (R | psiSet(i,iSpin,ik) ).absSqr();
                      */
                     PsiR psiR = R | psiSet(i,iSpin,ik);
                     rhoR(iSpin).plus_assign_ax (f, psiR.absSqr());

                  }
               }
            }
         }  // LoopMPI
      }
   }

   for (iSpin=0; iSpin < nSpin; iSpin++)  {
      SX_MPI_SOURCE("waves-k", TaskGroupMaster);
      SX_MPI_TARGET(TopLevel, TaskGroupAll);
      SxLoopMPI::sum (rhoR(iSpin));
      rhoR(iSpin) = R.symmetrize (rhoR(iSpin));
   }

   checkNorm (rhoR);

   return rhoR;
}


RhoR SxRho::computeRho (const Focc &/*focc*/,
                        const SxPW &/*psiSet0*/, const SxPW &/*psiSet1*/)
{
   RhoR rhoR;
   SX_EXIT;  // to be implemented
   return rhoR;
}

// --- *** UGLY QUICK&DIRTY COPY FROM REL-1-0 ***
//     *** TOBE REWRITEN IN A CLEAR WAY LATER *** {
RhoR &SxRho::atomicChargeDensity (const SxAtomicStructure &str, 
                                  const SxPseudoPot       &pot)
{
   SX_CHECK (rBasisPtr);
   
   // --- basis sets
   const SxRBasis &R = *rBasisPtr;                  // |R>
   const SxGBasis &G = R.getGBasis ();              // |G>
   SX_CHECK (str.cell.volume > 0.);

   SxMeshG rhoG (0., G);

   int nSpin = getNSpin ();

   // --- print-out
   cout << SX_SEPARATOR;
   cout << "|  Initial occupations\n";
   cout << SX_SEPARATOR;
   for (int iSpecies=0; iSpecies < pot.nSpecies; iSpecies++)  {
      cout << "| " << pot.elementName(iSpecies) << ":  ";
      for (int l=0; l <= pot.lMax(iSpecies); l++)  {
         switch (l)  {
            case 0 : sxprintf ("  s="); break;
            case 1 : sxprintf ("  p="); break;
            case 2 : sxprintf ("  d="); break;
            case 3 : sxprintf ("  f="); break;           
            default: sxprintf ("  %d=", l);
         }
         for (int iSpin=0; iSpin < nSpin; iSpin++)  {
            if (iSpin == 1) sxprintf ("/");
            sxprintf ("%2.2fe", pot.foccAtomicRho(iSpecies)(l)(iSpin));
         }
      }
      sxprintf ("\n");
   }	
   cout << SX_SEPARATOR;

   for (int iSpin=0; iSpin < nSpin; iSpin++)  {
      for (int iSpecies=0; iSpecies < pot.nSpecies; iSpecies++)  {
         rhoG += G.structureFactors(iSpecies)
               * pot.getAtomRhoG(G, iSpecies, iSpin);
      }
      rhoR(iSpin) = R.symmetrize (R | rhoG);
      rhoG.set(0.);
   }

   checkNorm ();
   normalizeRho ();

   return rhoR;
}
// --- *** UGLY QUICK&DIRTY COPY FROM REL-1-0 ***
//     *** TOBE REWRITEN IN A CLEAR WAY LATER *** }



//------------------------------------------------------------------------------
// Check for correct norm of charge density
//------------------------------------------------------------------------------
void SxRho::checkNorm ()
{
   checkNorm (rhoR);
}


Real8 SxRho::getNorm () const
{
   Real8 norm = 0.;
   for (int iSpin=0; iSpin < rhoR.getSize (); iSpin++)  {
      norm += rhoR(iSpin).sum() * rBasisPtr->dOmega;
   }
   return norm;
}


void SxRho::normalizeRho ()
{
   normalizeRho (&rhoR, getNorm ());
}


void SxRho::randomize ()
{
   for (int iSpin=0; iSpin < rhoR.getSize(); ++iSpin)  {
      rhoR(iSpin).randomize ();
      rhoR(iSpin).setBasis (rBasisPtr);
   }
   normalizeRho ();
}

void SxRho::checkNorm (RhoR &rhoIn) const
{
   if (nElectrons > 1e-10)  {
      Real8 norm = getNorm ();
      if (fabs (norm - nElectrons) > 0.001)  {
         cout << "\nThe norm of charge density is getting lost. " << endl;
         cout << "Norm = " << norm;
         cout << " (should be " << nElectrons << ")" << endl;
         normalizeRho (&rhoIn, norm);
      }
   }  else  {
      // sxprintf ("Skipping rho norm check.\n");
   }
}

void SxRho::normalizeRho (RhoR *rhoIn, Real8 rhoInNorm) const
{
   for (int iSpin=0; iSpin < rhoIn->getSize(); iSpin++)
      (*rhoIn)(iSpin) *= nElectrons/rhoInNorm;
}

//SxRho &SxRho::mixRhoLinearly (const RhoR &rhoIn, Real8 mixFactor)
//{
//   SX_CHECK (mixFactor >= 0. && mixFactor <= 1., mixFactor);
//
//   for (int iSpin=0; iSpin < control->nSpin; iSpin++)  {
//      rhoR(iSpin) = rhoR(iSpin)  *       mixFactor
//                  + rhoIn(iSpin) * (1. - mixFactor);
//   }
//
//   return *this;
//}



void SxRho::readRho (const SxBinIO &io)
{
   SX_CHECK (rBasisPtr);
   const SxRBasis &R  = *rBasisPtr;

   int meshSize = R.getMeshSize ();
   const SxVector3<int> mesh = R.getMesh ();

   int iSpin, nSpin = getNSpin ();

   try  {
      SxMatrix3<double> cell;
      SxVector3<int>    dim;
      rhoR = io.readMesh (&cell, &dim);
      if (nSpin <= 0) nSpin = int(rhoR.getSize ());
      for (iSpin=0; iSpin < nSpin; iSpin++)  rhoR(iSpin).setBasis (&R);


      if (rhoR.getSize() != nSpin)  {
         cout << SX_SEPARATOR;
         sxprintf ("| Charge density file does not contain the correct number "
                 "of spin channels!\n");
         sxprintf ("| The file contains %d spin components whereas %d are "
                 "expected.\n", int(rhoR.getSize()), nSpin);
         cout << SX_SEPARATOR;
         SX_EXIT;
      }

      if (mesh(0) != dim(0) || mesh(1) != dim(1) || mesh(2) != dim(2))
      {
         cout << SX_SEPARATOR;
         sxprintf ("| The mesh size %dx%dx%d of the input density is wrong.\n",
                  dim(0), dim(1), dim(2));
         sxprintf ("| (should be %dx%dx%d)\n", mesh(0), mesh(1), mesh(2));
         cout << SX_SEPARATOR;
         SX_EXIT;
      }

      if (rhoR(0).getSize() != meshSize)  {
         cout << SX_SEPARATOR;
         sxprintf ("| The mesh size of the input density file is wrong!\n");
         sxprintf ("| The input density has %d elements (%d expected).\n", 
                 (int)rhoR(0).getSize(), meshSize);
         cout << SX_SEPARATOR;
         SX_EXIT;
      }

      if ( (cell - R.cell).absSqr().sum() > 1e-5 )  {
         cout << SX_SEPARATOR;
         sxprintf ("| The geometry of the input density is wrong!\n");
         sxprintf ("| Input density: [[%g,%g,%g],[%g,%g,%g],[%g,%g,%g]]\n",
                 cell(0,0), cell(0,1), cell(0,2),
                 cell(1,0), cell(1,1), cell(1,2),
                 cell(2,0), cell(2,1), cell(2,2));
         sxprintf ("| should be:     [[%g,%g,%g],[%g,%g,%g],[%g,%g,%g]]\n",
                 R.cell(0,0), R.cell(0,1), R.cell(0,2),
                 R.cell(1,0), R.cell(1,1), R.cell(1,2),
                 R.cell(2,0), R.cell(2,1), R.cell(2,2));
         cout << SX_SEPARATOR;
         SX_EXIT;
      }

   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}


void SxRho::writeRho (SxBinIO &io) const
{
   SX_CHECK (rBasisPtr);
   const SxAtomicStructure *strPtr = rBasisPtr->hasGBasis ()
                                   ? rBasisPtr->getGBasis ().structPtr
                                   : NULL;
   try  {
      io.writeMesh (rhoR, rBasisPtr->cell, rBasisPtr->getMesh ());
      if (strPtr) strPtr->write (io);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

// --- SxDensity interface
void SxRho::operator= (const SxDensity &in)
{
   const SxRho &inRef = in.getRef<SxRho> ();
   rhoR.resize (inRef.rhoR.getSize ());
   for (int iSpin = 0; iSpin < rhoR.getSize (); ++iSpin)
      rhoR(iSpin) = inRef.rhoR(iSpin);
   if (!rBasisPtr) rBasisPtr = inRef.rBasisPtr;
}

void SxRho::operator+= (const SxDensity &x)
{
   const SxRho &xRef = x.getRef<SxRho> ();
   SX_CHECK (rhoR.getSize () == xRef.rhoR.getSize (),
             rhoR.getSize (), xRef.rhoR.getSize ());
   for (int iSpin = 0; iSpin < rhoR.getSize (); ++iSpin)
      rhoR(iSpin) += xRef.rhoR(iSpin);
}

void SxRho::operator-= (const SxDensity &x)
{
   const SxRho &xRef = x.getRef<SxRho> ();
   SX_CHECK (rhoR.getSize () == xRef.rhoR.getSize (),
             rhoR.getSize (), xRef.rhoR.getSize ());
   for (int iSpin = 0; iSpin < rhoR.getSize (); ++iSpin)
      rhoR(iSpin) -= xRef.rhoR(iSpin);
}

void SxRho::plus_assign_ax (double a, const SxDensity &x)
{
   const SxRho &xRef = x.getRef<SxRho> ();
   SX_CHECK (rhoR.getSize () == xRef.rhoR.getSize (),
             rhoR.getSize (), xRef.rhoR.getSize ());
   for (int iSpin = 0; iSpin < rhoR.getSize (); ++iSpin)
      rhoR(iSpin).plus_assign_ax(a, xRef.rhoR(iSpin));
}

void SxRho::plus_assign_aspin (double a, const SxDensity &x)
{
   const SxRho &xRef = x.getRef<SxRho> ();
   SX_CHECK (getNSpin () == 2, getNSpin ());
   SX_CHECK (xRef.getNSpin () == 1, xRef.getNSpin ());
   rhoR(0).plus_assign_ax( 0.5 * a, xRef.rhoR(0));
   rhoR(1).plus_assign_ax(-0.5 * a, xRef.rhoR(0));
}

/// Scalar product
double SxRho::operator| (const SxDensity &x) const
{
   const SxRho &xRef = x.getRef<SxRho> ();
   SX_CHECK (rhoR.getSize () == xRef.rhoR.getSize (),
             rhoR.getSize (), xRef.rhoR.getSize ());
   double res = 0.;
   //const SxGBasis &G = rBasisPtr->getGBasis ();
   for (int iSpin = 0; iSpin < rhoR.getSize (); ++iSpin)  {
      res += dot((*this)(iSpin), xRef(iSpin));
      // SxVector<PrecCoeffG> xG = G | xRef(iSpin), myG = G | ((*this)(iSpin));
      // res += 10. * dot (myG, (1. - G.g2/(2. + G.g2) ) * xG).re; 
   }
   return res / (double)rhoR(0).getSize ();
}

double SxRho::normSqr () const
{
   double res = 0.;
   for (int iSpin = 0; iSpin < rhoR.getSize (); ++iSpin)
      res += rhoR(iSpin).normSqr ();
   return res / (double)rhoR(0).getSize ();
}

SxDensity SxRho::operator- (const SxDensity &x) const
{
   const SxRho &xRef = x.getRef<SxRho> ();
   SX_CHECK (rhoR.getSize () == xRef.rhoR.getSize (),
             rhoR.getSize (), xRef.rhoR.getSize ());
   SxPtr<SxRho> res = SxPtr<SxRho>::create ();
   if (rBasisPtr)
      res->rBasisPtr = rBasisPtr;
   else
      res->rBasisPtr = xRef.rBasisPtr;
   int nSpin = getNSpin ();
   res->rhoR.resize (nSpin);
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)
      res->rhoR(iSpin) = rhoR(iSpin) - xRef.rhoR(iSpin);
   return SxDensity(res);
}

SxDensity SxRho::getCopy () const
{
   SxPtr<SxRho> res = SxPtr<SxRho>::create ();
   res->rBasisPtr = rBasisPtr;
   res->nElectrons = nElectrons;
   int nSpin = getNSpin ();
   res->rhoR.resize (nSpin);
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)
      res->rhoR(iSpin).copy(rhoR(iSpin));
   return SxDensity(res);
}

SxDensity SxRho::spin () const
{
   SX_CHECK (rhoR.getSize () == 2, rhoR.getSize ());
   SxPtr<SxRho> res = SxPtr<SxRho>::create ();
   res->rBasisPtr = rBasisPtr;
   res->nElectrons = -1.;
   res->rhoR.resize (1);
   res->rhoR(0) = rhoR(0) - rhoR(1);
   return SxDensity(res);
}

void SxRho::renormalize ()
{
   /*
   int nSpin = rhoR.getSize ();
   double norm = 0.;
   for (int iSpin=0; iSpin < nSpin; iSpin++)
      norm += rhoR(iSpin).sum ();
   norm *= dOmega;
   if (fabs (norm - nElectrons) > 1e-8 * nElectrons)  {
      double adjustRho = (nElectrons - norm) / (dOmega * nSpin);
      for (iSpin=0; iSpin < nSpin; iSpin++)
         rhoR(iSpin) += adjustRho;
   }
   */
   normalizeRho ();
}

void SxRho::syncMPI ()
{
#ifdef USE_LOOPMPI
   for (int iSpin = 0; iSpin < rhoR.getSize (); ++iSpin)  {
      if (rBasisPtr)
         rhoR(iSpin).resize (rBasisPtr->getNElements ());
      SxLoopMPI::bcast (rhoR(iSpin), 0);
      SX_VALIDATE_VECTOR(rhoR(iSpin));
   }
#endif
}

void SxRho::displaceHirshfeld (const SxAtomicStructure &toStr,
                               const SxArray<SxMeshG > &atomRhoG,
                               SxVector<double> *allAtomRptr)
{
   SX_CHECK (atomRhoG.getSize () > 0);
   SX_CHECK (rhoR.getSize () > 0);
   const SxGBasis &G = atomRhoG(0).getBasis<SxGBasis> ();
   const SxRBasis &R = rhoR(0).getBasis<SxRBasis> ();
   const SxAtomicStructure structure = G.getTau ();
   SX_CHECK (atomRhoG.getSize () == structure.getNSpecies (),
             atomRhoG.getSize (), structure.getNSpecies ());
   cout << "Displacing deformation density via Hirshfeld..." << endl;

   // --- get overlap of all atomic densities (if not available)
   bool haveAllAtom = allAtomRptr;
   if (!haveAllAtom)  {
      allAtomRptr = new SxVector<double> ();
      SxMeshG allAtom(0., G);
      SX_LOOP(is)  {
         allAtom += atomRhoG(is) * G.structureFactors(is);
      }
      *allAtomRptr = R | allAtom;
   }
   SxVector<double> &allAtomR = *allAtomRptr;

   // inverse of all-atom density, forced to be > 1e-4
   SX_LOOP(ir) allAtomR(ir) = 1. / sqrt(1e-8 + sqr(allAtomR(ir)));

   // --- idea: displace the Hirshfeld decomposition of the
   //     deformation density
   int nSpin = getNSpin ();
   SxVector<SxComplex16> dRhoG(G.ng, nSpin);
   dRhoG.set (0.);
   dRhoG.setBasis (G);

   // --- now extract contribution of each atom to the deformation
   //     density weighted by w(r) = thisAtom(r) / allAtoms(r)
   SX_LOOP2(is,ia)  {
      Coord dx = toStr(is, ia) - structure(is, ia);
      //if (dx.norm () < 1e-6) continue;
      SxMeshG myAtomG = atomRhoG(is) * G.getPhaseFactors ((int)is, (int)ia);
      // decompose density like the overlap of atomic densities
      SxVector<double> weight = (R | myAtomG) * allAtomR;
      // --- get the contribution and shift it in G-space

      for (int iSpin = 0.; iSpin < nSpin; ++iSpin)  {
         SxVector<double> myAtomR = weight * rhoR(iSpin);
         dRhoG.colRef(iSpin) += (G | myAtomR)
                              * (G.getPhaseFactors (dx) - 1.);
      }
   }
   // and add back the shifted contributions
   for (int iSpin = 0.; iSpin < nSpin; ++iSpin)  {
      rhoR(iSpin) += (R | dRhoG.colRef (iSpin)).real ();
   }
   if (!haveAllAtom) delete allAtomRptr;
}


SxArray<double>
SxRho::getHirshfeldEffVolumes (const SxArray<SxVector<double> > &atomRhoG,
                               const SxPAWPot &pawPot,
                               double gCut,
                               SxVector<double> *allAtomRptr)
{
   SX_CHECK (atomRhoG.getSize () > 0);
   SX_CHECK (rhoR.getSize () > 0);
   const SxGBasis &G = atomRhoG(0).getBasis<SxGBasis> ();
   const SxRBasis &R = rhoR(0).getBasis<SxRBasis> ();
   const SxAtomicStructure structure = G.getTau ();

   int nSpecies = structure.getNSpecies ();
   int nAtoms = structure.getNAtoms ();
   int nSpin = getNSpin ();

   SX_CHECK (atomRhoG.getSize () == nSpecies,
             atomRhoG.getSize (), nSpecies);

   SxArray<SxVector<double> > rhoAtomFree(nSpecies),
                                rhoAtomFreeR3(nSpecies),
                                shapeG(nSpecies),
                                shapeR3G(nSpecies),
                                rhoAtomG(nSpecies),
                                rhoAtomG3(nSpecies);

   SxArray<double> r3Free(nSpecies), effectiveVolumes(nAtoms);
   effectiveVolumes.set (0.);

   const SxRadBasis &rad  = pawPot.getRadBasis ();
   double rMax = rad.getRMax ();
   rMax = 20.;

   double gMax = sqrt(gCut)*1.05;

   int nrRad = 1000; // ad hoc
   SxRadialBasis radR(0.,rMax,nrRad,true,SxRadialBasis::Linear);
   int ngRad = 1000; // ad hoc
   SxRadialBasis radG(0., gMax, ngRad,false, SxRadialBasis::Linear);
   const SxVector<double> &rVals = radR.getRadFunc ();

   for (int is = 0; is < nSpecies; is++)  {
      // calculate and store free atomic "volumes"
      rhoAtomFree(is) = radR | pawPot.rhoInit(is);
      rhoAtomFreeR3(is) = rhoAtomFree(is) * rVals.cub ();

      shapeG(is) = radG | rhoAtomFree(is);
      shapeR3G(is) = radG | rhoAtomFreeR3(is);

      rhoAtomG(is) = G | shapeG(is);
      rhoAtomG3(is) = G | shapeR3G(is);
      SxVector<double> rhoAtomR3 = R | rhoAtomG3(is);

      SxVector<double> unity(rhoAtomR3.getSize());
      unity.set (1.);

      r3Free(is) = dot(unity, rhoAtomR3);
   }

   // --- get overlap of all atomic densities (if not available)
   bool haveAllAtom = allAtomRptr;
   if (!haveAllAtom)  {
      SxMeshG allAtom(0., G);
      SX_LOOP(is)  {
         allAtom += rhoAtomG(is) * G.structureFactors(is);
      }
      allAtomRptr = new SxVector<double> ();
      *allAtomRptr = R | allAtom;
   }
   SxVector<double> &allAtomR = *allAtomRptr;

   // all-atom density, forced to be >= 1e-4
   SX_LOOP(ir) allAtomR(ir) = sqrt(1e-10 + sqr(allAtomR(ir)));

   // --- now extract contribution of each atom to the deformation
   //     density weighted by w(r) = thisAtom(r) / allAtoms(r)
   for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
      int is = structure.getISpecies(iAtom);
      int ia = iAtom - structure.atomInfo->offset(is);

      SxMeshG myAtomG = rhoAtomG3(is) * G.getPhaseFactors ((int)is, (int)ia);
      SxVector<double> myAtomR = (R | myAtomG);

      // Also force this atom's density to be >= 1e-4
      //SX_LOOP(ir) myAtomR(ir) = sqrt(1e-10 + sqr(myAtomR(ir)));
      SxVector<double> weight = myAtomR / allAtomR;
      double r3Real = 0.;
      for (int iSpin = 0.; iSpin < nSpin; ++iSpin)  {
         r3Real += dot(weight, rhoR(iSpin));// * sqrt(FOUR_PI);
      }
      effectiveVolumes(iAtom) = r3Real / r3Free(is);
   }
   if (!haveAllAtom) delete allAtomRptr;
   return effectiveVolumes;
}

SxArray<double>
SxRho::getHirshfeldEffCharges (const SxArray<SxVector<double> > &atomRhoG,
                               const SxPAWPot &pawPot,
                               double gCut,
                               SxVector<double> *allAtomRptr)
{
   SX_CHECK (atomRhoG.getSize () > 0);
   SX_CHECK (rhoR.getSize () > 0);
   const SxGBasis &G = atomRhoG(0).getBasis<SxGBasis> ();
   const SxRBasis &R = rhoR(0).getBasis<SxRBasis> ();
   const SxAtomicStructure structure = G.getTau ();

   int nSpecies = structure.getNSpecies ();
   int nAtoms = structure.getNAtoms ();
   int nSpin = getNSpin ();

   SX_CHECK (atomRhoG.getSize () == nSpecies,
             atomRhoG.getSize (), nSpecies);

   SxArray<SxVector<double> > rhoAtomFree(nSpecies),
                                rhoAtomFreeR3(nSpecies),
                                shapeG(nSpecies),
                                shapeR3G(nSpecies),
                                rhoAtomG(nSpecies),
                                rhoAtomG3(nSpecies);

   SxArray<double> r3Free(nSpecies), effectiveCharges(nAtoms);
   effectiveCharges.set (0.);

   const SxRadBasis &rad  = pawPot.getRadBasis ();
   double rMax = rad.getRMax ();

   double gMax = sqrt(gCut)*1.05;

   int nrRad = 1000; // ad hoc
   SxRadialBasis radR(0.,rMax,nrRad, true, SxRadialBasis::Linear);
   int ngRad = 1000; // ad hoc
   SxRadialBasis radG(0., gMax, ngRad, false, SxRadialBasis::Linear);
   const SxVector<double> &rVals = radR.getRadFunc ();

   for (int is = 0; is < nSpecies; is++)  {
      // calculate and store free atomic "volumes"
      rhoAtomFree(is) = radR | pawPot.rhoInit(is);
      rhoAtomFreeR3(is) = rhoAtomFree(is) * rVals.cub ();

      shapeG(is) = radG | rhoAtomFree(is);
      shapeR3G(is) = radG | rhoAtomFreeR3(is);

      rhoAtomG(is) = G | shapeG(is);
      rhoAtomG3(is) = G | shapeR3G(is);
      SxVector<double> rhoAtomR3 = R | rhoAtomG3(is);

      SxVector<double> unity(rhoAtomR3.getSize());
      unity.set (1.);

      r3Free(is) = dot(unity, rhoAtomR3);
   }

   // --- get overlap of all atomic densities (if not available)
   bool haveAllAtom = allAtomRptr;
   if (!haveAllAtom)  {
      allAtomRptr = new SxVector<double> ();
      SxMeshG allAtom(0., G);
      SX_LOOP(is)  {
         allAtom += atomRhoG(is) * G.structureFactors(is);
      }
      *allAtomRptr = R | allAtom;
   }
   SxVector<double> &allAtomR = *allAtomRptr;

   // all-atom density, forced to be >= 1e-4
   SX_LOOP(ir) allAtomR(ir) = sqrt(1e-8 + sqr(allAtomR(ir)));

   // --- now extract contribution of each atom to the deformation
   //     density weighted by w(r) = thisAtom(r) / allAtoms(r)
   for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
      int is = structure.getISpecies(iAtom);
      int ia = iAtom - structure.atomInfo->offset(is);

      SxMeshG myAtomG = rhoAtomG3(is) * G.getPhaseFactors ((int)is, (int)ia);
      SxVector<double> myAtomR = (R | myAtomG);

      // Also force this atom's density to be >= 1e-4
      SX_LOOP(ir) myAtomR(ir) = sqrt(1e-8 + sqr(myAtomR(ir)));
      SxVector<double> weight = myAtomR / allAtomR;
      double r3Real = 0.;
      for (int iSpin = 0.; iSpin < nSpin; ++iSpin)  {
         r3Real += dot(weight, rhoR(iSpin));
      }
      effectiveCharges(iAtom) = pawPot.valenceCharge(is)*r3Real/r3Free(is) - pawPot.valenceCharge(is);
   }
   if (!haveAllAtom) delete allAtomRptr;
   return effectiveCharges;
}
