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

#include <SxSpinConstraint.h>
#include <SxFermi.h>
#include <SxPAWSet.h>
#include <SxAtomicStructure.h>
#include <SxFileParser.h>
#include <SxSimpleParser.h>
#include <SxEigensystem.h>

namespace Timer {
   enum SpinConstraintTimers { OptimizeNu, NuIntro , PAWInt, OmegaInit,
                               SpinsDiag, NuDiagH, OmegaRot, OmegaN, NuRest};
}

SX_REGISTER_TIMERS (Timer::SpinConstraintTimers)
{
   using namespace Timer;
   regTimer (OptimizeNu, "optimize nu");
   regTimer (NuIntro   , "opt. nu startup");
   regTimer (PAWInt    , "spin constr. PAW int");
   regTimer (OmegaInit , "setup omega_nn");
   regTimer (SpinsDiag , "spinsDiag total");
   regTimer (NuDiagH   , "spinsDiag H diag");
   regTimer (OmegaRot  , "spinsDiag rot. omega");
   regTimer (OmegaN    , "spinsDiag omega_n");
   regTimer (NuRest    , "spinsDiag rest");
}

SxSpinConstraint::SxSpinConstraint ( const SxPtr<SxPAWPot> &pawPotPtrIn)
   : wvPtr (NULL)
{
   pawPotPtr = pawPotPtrIn;

}

void SxSpinConstraint::read (const SxSymbolTable *topLvl,
                             const SxAtomicStructure &structure)
{
   // --- create variables for spin constraint
   targetSpin.resize (structure.getNAtoms());
   targetSpin.set (0.);
   constrainedAtom.resize (structure.getNAtoms ());
   constrainedAtom.set (false);
   // --- read the target spins from input.sx
   SYMBOLPARSE(topLvl)  {
      FOREACH_SYMBOLGROUP("spinConstraint")  {
         SxString constraintFile = SYMBOLGET("file") || "";
         if (constraintFile.getSize () > 0)  {
            SxFileParser fp(constraintFile);
            fp.topic ("spin constraints"); // get nicer errors
            SX_LOOP(ia)  {
               fp.skipWhite ();
               if (fp.reads("X") || fp.reads("x"))  {
                  constrainedAtom(ia) = false;
                  fp.nextLine ();
               } else {
                  constrainedAtom(ia) = true;
                  fp >> targetSpin(ia);
               }
            }
         } else {
            const SxArray<SxString> &labels = structure.getLabels ();

            SxString labelConstraint = SYMBOLGET("label");
            double   spinConst       = SYMBOLGET("constraint");

            for (int iAtom = 0; iAtom < structure.getNAtoms(); iAtom++) {
               if ( labelConstraint == labels(iAtom) ) {
                  targetSpin(iAtom) = spinConst;
                  constrainedAtom(iAtom) = true;
               }
            }
         }
      }
   }
   cout << SX_SEPARATOR;
   cout << "| Spin constraints" << endl;
   cout << SX_SEPARATOR;
   SX_LOOP(ia)  {
      sxprintf("| Atom %3d ",int(ia+1));
      if (constrainedAtom(ia))
         cout << "spin=" << targetSpin(ia);
      else
        cout << "not constrained.";
      cout << endl;
   }
   cout << SX_SEPARATOR;
   checkSym (structure);
}

namespace {
   class SxConstrainedInfo {
      public:
         int iAtom, is;
         ssize_t ipStart;
   };

   SxArray<SxConstrainedInfo>
   getInfoAtom(const SxArray<bool> &constrainedAtom,
               const SxAtomicStructure &structure,
               const SxPAWPot &pot)
   {
      int nAtomsConstr = 0;
      SX_LOOP(iAtom)
         if (constrainedAtom(iAtom)) nAtomsConstr++;

      SxArray<SxConstrainedInfo> infoAtom(nAtomsConstr);
      int iAtomC = 0;
      ssize_t ipStart = 0;
      SX_LOOP (iAtom) {
         int is = structure.getISpecies((int)iAtom);
         int npl = pot.getNProj(is);
         if (constrainedAtom(iAtom)) {
            infoAtom(iAtomC).is = is;
            infoAtom(iAtomC).iAtom = (int)iAtom;
            infoAtom(iAtomC).ipStart = ipStart;
            iAtomC++;
         }
         ipStart += npl;
      }
      return infoAtom;
   }

}

SxVector<SxComplex16>
SxSpinConstraint::getOmegaNN (const PsiRef &psi,
                              ssize_t ipBasis,
                              int is) const
{
   const SxPAWPot &pawPot = *pawPotPtr;
   int nStates = (int)psi.getNCols ();
   int npt = pawPot.getNProjType (is);
   int npl = pawPot.getNProj(is);
   int iSpin = psi.auxData.iSpin;
   const SxVector<double> &omegaIJ = pawPot.omegaPAW(is);
   // --- get compact copy of projectors for this atom
   //SxVector<PrecCoeffG> projAtom
   SxVecRef<PrecCoeffG,SubMatrix> projAtom
      = psi.getRef<SubMatrix> (ipBasis, npl, nStates);

   // --- apply PAW volume operator on projAtom (and transpose)
   SxVector<PrecCoeffG> projOmega (nStates, npl);
   projOmega.set (0.);
   for (int ipt = 0, ipl0 = 0; ipt < npt; ipt++) {
      int li = pawPot.lPhi(is)(ipt);
      int nm = 2 * li + 1;
      for (int jpt = 0, jpl0 = 0; jpt < npt; jpt++) {
         int lj = pawPot.lPhi(is)(jpt);
         if (li == lj) {
            double integral = omegaIJ(ipt, jpt);
            if (iSpin == 1) integral = -integral;
            for (int im = 0; im < nm; ++im)  {
               for (int iState = 0; iState < nStates; ++iState)
                  projOmega(iState, ipl0 + im) += integral * projAtom(jpl0 + im, iState).conj ();
            }
         }
         jpl0 += 2 * lj + 1;
      }
      ipl0 += nm;
   }
   // get matrix elements
   return projOmega ^ projAtom;
}

SxArray<SxVector<SxComplex16> >
SxSpinConstraint::getOmegaN (const SxVecRef<SxComplex16> &proj) const
{
   SX_CLOCK (Timer::OmegaN);
   SX_CHECK (pawPotPtr);
   SX_CHECK (wvPtr);
   int nAtom = (int)targetSpin.getSize ();
   int iSpin = proj.auxData.iSpin;
   ssize_t nStates = proj.getNCols ();
   const SxPAWPot &pawPot = *pawPotPtr;
   const SxArray<SxVector<double> > &omegaIJ = pawPot.omegaPAW;
   const SxAtomicStructure &structure = wvPtr->getGkBasisPtr ()->getTau ();

   SxArray<SxVector<SxComplex16> > result(nAtom);
   SX_LOOP(iAtom)
      if (constrainedAtom(iAtom)) result(iAtom).resize (nStates);

#ifdef USE_OPENMP
#  pragma omp parallel
#endif
   {
      // --- loop over atoms
      ssize_t ipBasis = 0;
      for (int iAtom = 0; iAtom < nAtom; iAtom++) {
         int is = structure.getISpecies(iAtom);
         int npt = pawPot.getNProjType (is);
         int npl = pawPot.getNProj(is);
         if (!constrainedAtom(iAtom)) {
            ipBasis += npl;
            continue;
         }
         // --- PAW volume operator expectation value for each state
#ifdef USE_OPENMP
#        pragma omp for nowait
#endif
         for (ssize_t iState = 0; iState < nStates; ++iState)  {
            SxComplex16 omega = 0.;
            for (int ipt = 0, ipl0 = 0; ipt < npt; ipt++) {
               int li = pawPot.lPhi(is)(ipt);
               int nm = 2 * li + 1;
               for (int jpt = 0, jpl0 = 0; jpt < npt; jpt++) {
                  int lj = pawPot.lPhi(is)(jpt);
                  if (li == lj) {
                     double integral = omegaIJ(is)(ipt, jpt);
                     if (iSpin == 1) integral = -integral;
                     SxComplex16 mSum = 0.;
                     for (int im = 0; im < nm; ++im)
                        mSum += proj(ipBasis + jpl0 + im, iState).conj ()
                              * proj(ipBasis + ipl0 + im, iState);
                     omega += mSum * integral;
                  }
                  jpl0 += 2 * lj + 1;
               }
               ipl0 += nm;
            }
            result(iAtom)(iState) = omega;
         }
         // next atom
         ipBasis += npl;
      }
   }
   return result;
}

void SxSpinConstraint::setNuA(double value)
{
   SX_CHECK (targetSpin.getSize () > 0);
   nuA.resize (targetSpin.getSize ());
   nuA.set (value);
}
void SxSpinConstraint::setConstraints(const SxVecRef<double> &targetSpinIn,
                                      const SxAtomicStructure &str)
{
   SX_CHECK(targetSpinIn.getSize () == str.getNAtoms (),
            targetSpinIn.getSize (), str.getNAtoms ());
   targetSpin.resize (str.getNAtoms ());
   SX_LOOP(ia)
      targetSpin(ia) = constrainedAtom(ia) ? targetSpinIn(ia) : 0.;
   checkSym (str);
}

void SxSpinConstraint::checkSym (const SxAtomicStructure &str) const
{
   if (!str.cell.symGroupPtr) return;
   const SxSymGroup &syms =  *str.cell.symGroupPtr;
   if (syms.getNSymmorphic () == 1) return;

   SxGrid grid(str, 10);
   bool fatal = false;
   SxArray2<bool> fatalIJ;
   for (int iSym = 0; iSym < syms.getNSymmorphic (); ++iSym)  {
      SxAtomInfo::ConstPtr map = str.match (grid,
                                            syms.getSymmorphic (iSym) ^ str);
      if (!map)  {
         cout << "Symmetry inconsistency!" << endl;
         cout << "Please check the symmetries of your structure." << endl;
         cout << "If you are not close to the numerical noise limit ("
              << str.epsEqual << ")," << endl
              << "contact the developers with your structure file" << endl;
         SX_EXIT;
      }
      SX_LOOP(ia)  {
         int ja = (int)map->parentMap(ia);
         if (fabs(targetSpin(ia) - targetSpin(ja)) > 1e-3
             || (constrainedAtom(ia)? 0 : 1) != (constrainedAtom(ja)? 0 : 1))  {
            if (!fatal)  {
               cout << endl << SX_SEPARATOR;
               cout << "ERROR: Symmetry inconsistent spin constraints:" << endl;
               cout << SX_SEPARATOR;
               int nAtoms = str.getNAtoms ();
               fatalIJ.reformat (nAtoms, nAtoms);
               fatalIJ.set (false);
            }
            if (!fatalIJ(ia,ja))  {
               cout << "Atom " << (ia + 1) << " @ " << str(ia) << " and" << endl
                    << "atom " << (ja + 1) << " @ " << str(ja) << " are symmetry equivalent"
                    << endl;
               cout << "Target spin " << (ia + 1) << ": " << targetSpin(ia);
               if (!constrainedAtom(ia)) cout << " (not constrained)";
               cout << endl;
               cout << "Target spin " << (ja + 1) << ": " << targetSpin(ja);
               if (!constrainedAtom(ja)) cout << " (not constrained)";
               cout << endl;
               fatalIJ(ia, ja) = fatalIJ(ja, ia) = true;
            }
            fatal = true;
         }
      }
   }
   if (fatal)  {
      cout << "=================" << endl;
      cout << "--- WHAT TO DO?" << endl;
      cout << "=================" << endl;
      cout << "* Fix the target spins to a symmetry-compatible set"
           << endl;
      cout << "* use labels to reduce the symmetry" << endl;
      cout << "* specify symmetry in the input file" << endl;
      cout << "  To switch off all symmetries, insert" << endl
           << "     symmetry { operator { S = [[1,0,0],[0,1,0],[0,0,1]]; } }"
           << endl << "  into the structure {} group" << endl;
      SX_QUIT;
   }
   cout << "Spin constraints: symmetry check OK" << endl;
}

double SxSpinConstraint::calcKappaOpt (const SxVector<double> &MSpinIn, const SxVector<double> &MSpinPlusIn, const double kappa)
{
   if ((MSpinIn - targetSpin).normSqr () / double(MSpinIn.getSize ()) > 1e-2)  {
      cout << "MSpinIN   and MSpinPlusIN   = " << MSpinIn << "  " << MSpinPlusIn << endl;
      cout << "targetSpin  =  " << targetSpin << endl;
   }
   double sumK = 0, sumK2 = 0; 
   sumK = dot( (MSpinIn - targetSpin) , (MSpinPlusIn - MSpinIn) );
   sumK2 = (MSpinPlusIn - MSpinIn).normSqr();
   return - (sumK * kappa)/(sumK2);
}

SxVector<double>
SxSpinConstraint::computeNu (const SxPtr<SxPAWSet> &wavesPtr, SxFermi &fermi,
                             const Real8 ekt, double epsilon)
{
   SX_CLOCK(Timer::OptimizeNu);
   SX_START_TIMER (Timer::NuIntro);
   SxPAWSet &waves = *wavesPtr;
   double kappa = 0.01, kappaOpt = 0, meanError;
   int nSpin = fermi.getNSpin ();
   int nk = fermi.getNk ();
   int nStates = fermi.getNStates ();
   Eps eps0(nStates,nSpin,nk);
   ssize_t nAtom = nuA.getSize();

   wvPtr = &waves;
   for (int ik = 0; ik < nk; ik++) {
      SX_MPI_LEVEL("waves-k");
      if (SxLoopMPI::myWork(ik)) {
         for (int iSpin = 0; iSpin < nSpin; iSpin++)
            eps0(iSpin,ik) <<= fermi.eps(iSpin,ik);
      }
   }
   SxVector<double> MSpin(nAtom), MSpinPlus(nAtom);
   // --- iteration loop for spin constraint
   int nStep = 100;
   SxVector<double> searchOld, deltaSpinOld;
   double meanOld = 0.;
   SX_STOP_TIMER(Timer::NuIntro);
   for (int iStep = 0; iStep < nStep; iStep++) {
      //cout << "nuA=" << nuA << endl;
      MSpin = getSpinsDiag(eps0, fermi, ekt);
      //cout << "MSpin=" << MSpin << endl;
      SxVector<double> deltaSpin = MSpin - targetSpin;
      SX_LOOP(ia) if (!constrainedAtom(ia)) deltaSpin(ia) = 0.;
      SxVector<double> searchNu = deltaSpin;

      meanError = deltaSpin.normSqr() * 1./double(nAtom);
      if ( meanError < epsilon ) break;

      if (iStep > 0) {
         // Fletcher-Reeves
         double gamma = meanError / meanOld;
         // Polak-Ribiere
         //double gamma = (meanError - dot(deltaSpinOld,deltaSpin)/double(nAtom)) / meanOld;
         // Hestenes-Stiefel
         //double gamma = dot(deltaSpinOld-deltaSpin,deltaSpin)/dot(deltaSpinOld - deltaSpin,searchOld);
         cout << "gamma=" << gamma << " R_M^2=" << meanError << endl;
         searchNu += gamma * searchOld;
      } else {
         cout << "gamma=X R_M^2=" << meanError << endl;
      }

      nuA += kappa * searchNu;
      MSpinPlus = getSpinsDiag(eps0, fermi, ekt);

      // --- calculate optimal kappa with test Spin^(i+1)
      kappaOpt = calcKappaOpt(MSpin, MSpinPlus, kappa);
      cout << "kappa opt=" << kappaOpt << endl;
      // --- update nuA for next iteration
      if (kappaOpt * searchNu.maxval () > 3. * ekt) {
         // note: factor 3 * ekt is not tested
         kappaOpt = 3. * ekt / searchNu.maxval ();
         cout << "Restricting step size" << endl;
      }
      nuA += (kappaOpt - kappa) * searchNu;
      //cout << "delta Nu=" << nuA << endl;
      // check for abort criterion
      searchOld = std::move (searchNu);
      deltaSpinOld = std::move (deltaSpin);
      meanOld = meanError;

      // --- adjust trial step
      double g = 1.5 * fabs(kappaOpt) / kappa;
      if (g > 2.) g = 2.;
      if (g < 0.5) g = 0.5;
      kappa *= pow(g,0.7);
      cout << "kappa=" << kappa << endl;
   } // iStep

   // rotate waves
   getSpinsDiag (eps0, fermi, ekt, wavesPtr);
   wvPtr = NULL;
   return nuA;
}

SxVector<double> SxSpinConstraint::getSpinsDiag ( const Eps &eps0In, SxFermi &fermiIn, const double &ektIn , SxPtr<SxPAWSet> wavesPtr ) 
{
   SX_CHECK (wvPtr);
   SX_CLOCK (Timer::SpinsDiag);
   int nSpin = fermiIn.getNSpin ();
   int nk = fermiIn.getNk ();
   int nStates = fermiIn.getNStates ();
   ssize_t nAtom = nuA.getSize();
   SxArray<Eps> omega(nAtom);
   SX_MPI_LEVEL("waves-k");

   // --- selective resize of omega
   SX_LOOP(iAtom)  {
      // only constrained atoms
      if (constrainedAtom(iAtom))  {
         omega(iAtom).bundle.resize (nk);
         for (int ik = 0; ik < nk; ik++) {
            // only ik's for this MPI task
            if (SxLoopMPI::myWork(ik))  {
               omega(iAtom).bundle(ik).resize (nSpin);
               SX_LOOP(iSpin) omega(iAtom).bundle(ik)(iSpin).resize (nStates);
            }
         }
      }
   }

   // --- get condensed list of constrained atoms for openmp parallelization
   SxArray<SxConstrainedInfo> infoAtom
      = getInfoAtom (constrainedAtom, wvPtr->getGkBasisPtr ()->getTau (),
                     *pawPotPtr);
   ssize_t nAtomsConstr = infoAtom.getSize ();

   SxSymEigensystem<SxComplex16> eig;
   for( int iSpin = 0; iSpin < nSpin; iSpin++) {
      for (int ik = 0; ik < nk; ik++) {
         if (!SxLoopMPI::myWork(ik)) continue;

         // --- set up new Hamiltonian matrix for current deltaNu
         SxVector<SxComplex16> Hnew(nStates, nStates);
         // spin constraint contribution
         const SxVecRef<SxComplex16> &psi = (*wvPtr)(iSpin, ik);
         ssize_t ipBasis = psi.getBasis<SxPAWBasis> ().gBasis->ng;
         SX_START_TIMER(Timer::OmegaInit);
#ifdef USE_OPENMP
         SxArray<SxVector<SxComplex16> > extraH(omp_get_max_threads ());
#        pragma omp parallel
         {
            SX_ALLOC_CACHE;
            int nThreads = omp_get_num_threads ();
            int iThread = omp_get_thread_num ();
            // --- get contribution to Hamiltonian from constraints
            //     separately for threads -> extraH
            if (iThread > 0)  {
               extraH(iThread).reformat (nStates, nStates);
               extraH(iThread).set (0.);
            } else {
               Hnew.set (0.);
               // diagonal elements
               SX_LOOP(n) Hnew(n,n) += eps0In(n,iSpin,ik);
            }
            SxVector<SxComplex16> &myH
               = (iThread == 0) ? Hnew : extraH(iThread);
#           pragma omp for
            for (ssize_t iAtomC = 0; iAtomC < nAtomsConstr; ++iAtomC) {
               int iAtom = infoAtom(iAtomC).iAtom;
               SxVector<SxComplex16> omegaNN
                  = getOmegaNN (psi, infoAtom(iAtomC).ipStart + ipBasis,
                                infoAtom(iAtomC).is);

               // --- now add to Hamiltonian
               myH.plus_assign_ax (nuA(iAtom), omegaNN);
            }
            // --- now sum over thread contributions
            if (nThreads > 1)  {
               ssize_t nij = (Hnew.getSize () - 1) / nThreads + 1;
               while (nij & 3) nij++;
               ssize_t ijFrom = iThread * nij;
               ssize_t ijTo = min(ijFrom + nij, Hnew.getSize ());
               for (int jThread = 1; jThread < nThreads; ++jThread)  {
                  const SxVector<SxComplex16> &Hj = extraH(jThread);
#                 pragma omp simd
                  for (ssize_t ij = ijFrom; ij < ijTo; ++ij)
                     Hnew(ij) += Hj(ij);
               }
            }
         }
         extraH.resize (0);
#else
         Hnew.set (0.);
         // diagonal elements
         SX_LOOP(n) Hnew(n,n) += eps0In(n,iSpin,ik);
         for (ssize_t iAtomC = 0; iAtomC < nAtomsConstr; ++iAtomC) {
            int iAtom = infoAtom(iAtomC).iAtom;
            ssize_t ipTot0 = infoAtom(iAtomC).ipStart + ipBasis;
            SxVector<SxComplex16> omegaNN
               = getOmegaNN (psi, ipTot0, infoAtom(iAtomC).is);
            Hnew.plus_assign_ax(nuA(iAtom), omegaNN);
         }
#endif
         SX_STOP_TIMER(Timer::OmegaInit);

         // --- diagonalize
         { 
            SX_CLOCK (Timer::NuDiagH);
            eig.compute (std::move(Hnew), All, true);
            SX_LOOP(in) fermiIn.eps(in,iSpin,ik) = eig.vals(in);
         }

         // --- rotate omegann according to new eigenfunctions
         {
            SX_CLOCK(Timer::OmegaRot);
            ssize_t ng = psi.getBasis<SxPAWBasis> ().gBasis->ng;
            ssize_t nProj = psi.getBasis<SxPAWBasis> ().pBasis->getNElements ();
            SxVecRef<SxComplex16, SubMatrix> proj
               = psi.getRef<SubMatrix> (ng, nProj, nStates);
            proj.auxData.iSpin = (char)iSpin;
            SxVector<SxComplex16> projRot = proj ^ eig.vecs;
            SxArray<SxVector<SxComplex16> > omegaN = getOmegaN(projRot);
            SX_LOOP(iAtom)
               if (constrainedAtom(iAtom))
                  omega(iAtom)(iSpin, ik) = std::move (omegaN(iAtom));
         }

         if (wavesPtr) (*wavesPtr)(iSpin, ik).rotate (eig.vecs);
      }//ik
   }//iSpin
   SX_CLOCK(Timer::NuRest);
   // tell every ik after MPI for fermiDistribution
   fermiIn.eps.synMPI();
   // --- get the new fermi occupations for current deltaNu
   fermiIn.fermiDistribution(ektIn);

   SxVector<double> MSpin(nAtom);
   for (int iAtom = 0; iAtom < nAtom; iAtom++) {
      MSpin(iAtom) = 0;
      if (!constrainedAtom(iAtom)) continue;
      for (int iSpin = 0; iSpin < nSpin; iSpin++) {
         for (int ik = 0; ik < nk; ik++) {
            if (SxLoopMPI::myWork(ik)) {
               for (int in = 0; in < fermiIn.getNStates(ik); ++in)  {
                   MSpin(iAtom) += omega(iAtom)(in,iSpin,ik) * fermiIn.focc(in,iSpin,ik) * fermiIn.kpPtr->weights(ik);
               }//in
            }//LoopMPI
         }//ik
      }//iSpin
   }//iAtom
   SX_MPI_SOURCE ("waves-k", TaskGroupMaster);
   SX_MPI_TARGET (TopLevel, TaskGroupAll);
   SxLoopMPI::sum(MSpin);
   //SX_LOOP(iAtom) cout << "MSpin(" << iAtom << ") = " << setprecision(8) << MSpin(iAtom) << endl;

   return MSpin;
}
