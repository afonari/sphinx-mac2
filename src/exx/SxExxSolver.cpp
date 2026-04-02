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

#include <SxExxSolver.h>
#include <SxPW.h>
#include <SxMuPW.h>
#include <SxRhoMixer.h>
#include <SxParser.h>
#include <SxConstants.h>
#include <SxFile.h>
#include <SxRhoMixerG.h> // This is Blazejs G-space version of the mixer.
//#include <SxPulay.h>
#include <SxChi0.h>
#include <SxPWHamiltonian.h>
#include <SxFockGk.h>
#include <SxSymFockGk.h>
#include <SxFileIO.h>
#include <SxTextIO.h>
#include <SxIO.h>

SxExxSolver::SxExxSolver (const SxAtomicStructure &structureIn,
                          const SxSymbolTable *table)
         : SxHamSolver(structureIn, table)
{
   SX_CHECK(wavesPtr);
   // --- initialize n123inv table of |G+k> basises
   SxXC::XCFunctional xc = SxXC::getXCFunctional (table->topLevel ());
   if (   xc == SxXC::EXX
       || xc == SxXC::EXX_LDA
       || xc == SxXC::EXX_PBE
       || xc == SxXC::EXX_PBE_WB
       || xc == SxXC::Slater
       || xc == SxXC::KLI)
      wavesPtr->getGkBasis().setupN123inv ();
}

bool SxExxSolver::isExx (const SxSymbolTable *cmd)
{
   if (cmd == NULL) return false;
   SxString name = cmd->getName();
   return (   name == "SH" 
           || name == "STEXX"
           || name == "gsEXX"
           || name == "orbDepPot");
}

bool SxExxSolver::isRegistered (const SxSymbolTable *cmd) const
{
   if (isExx (cmd)) return true;
   return (SxHamSolver::isRegistered (cmd));
}

void SxExxSolver::execute (const SxSymbolTable *cmd, bool calc)
{
   SX_CHECK (cmd);

   SxString str = cmd->getName();
   bool     storeWaves, storeRho;

   rhoMixing = foccMixing = deltaT = dPsiConv = dEnergy = dEps = dRelEps
             = ekt = 0.;
   maxSteps  = 1;

   try  {
      storeWaves =  (cmd->contains("noWavesStorage"))
                 ? !(cmd->get("noWavesStorage")->toAttribute())
                 :   true;
      storeRho   =  (cmd->contains("noRhoStorage"))
                 ? !(cmd->get("noRhoStorage")->toAttribute())
                 :   true;
   }  catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   if (hamPtr)  {
      SX_CHECK (dynamic_cast<SxPWHamiltonian*>(hamPtr.getPtr ()));
      SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian&>(*hamPtr);
      if (!dynamic_cast<SxExx*>(H.xcPtr.getPtr ()))  {
         // --- replace standard XC with SxExx
         SxPtr<SxExx> exxPtr = SxPtr<SxExx>::create (*H.xcPtr);
         H.xcPtr = exxPtr;
         exxPtr->read (cmd, G);
      }
   }

   if (str == "STEXX")          stEXX (cmd, calc);
   else if (str == "gsEXX")     gsEXX (cmd, calc);
   else if (str == "orbDepPot") odp (cmd, calc);
   else                         SxHamSolver::execute (cmd, calc);

   // --- write data to disc
   if (calc)  {
      energy = hamPtr->getEnergy ();
      writeData (storeWaves, storeRho);
      printTiming (/* restart= */true);
   }
}

// ----------------------- Sternheimer -----------------------
PsiG
SxExxSolver::sternheimerStDesc_write (const PsiRef &h1psi0,
                                      const PsiRef &psi1Start,
                                      bool write)
{
   SX_CHECK (wavesPtr && dynamic_cast<SxPW*>(wavesPtr.getPtr ()));
   SxPW &waves = *SxPtr<SxPW>(wavesPtr);
#ifndef NDEBUG
   SxGkBasis &Gk = waves.getGkBasis();
#endif
   SX_CHECK (dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ()));
   SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian &>(*hamPtr);
   SX_CLOCK      (Timer::Sternheimer);
   SX_CHECK (h1psi0.getSize() == psi1Start.getSize(),
             h1psi0.getSize(), psi1Start.getSize());
   SX_CHECK (h1psi0.auxData.iSpin == psi1Start.auxData.iSpin,
             (int)h1psi0.auxData.iSpin, (int)psi1Start.auxData.iSpin);
   SX_CHECK (h1psi0.auxData.ik == psi1Start.auxData.ik,
             h1psi0.auxData.ik, psi1Start.auxData.ik);
   SX_CHECK (h1psi0.getBasisPtr() == &(Gk(h1psi0.auxData.ik)),
             (size_t)h1psi0.getBasisPtr(),
             (size_t)&(Gk(h1psi0.auxData.ik)));
   SX_CHECK (h1psi0.getBasisPtr() == psi1Start.getBasisPtr(),
             (size_t)h1psi0.getBasisPtr(),
             (size_t)psi1Start.getBasisPtr());
   int     i      = h1psi0.auxData.iOther[0];
   int     iSpin  = h1psi0.auxData.iSpin;
   int     ik     = h1psi0.auxData.ik;

   double  eps0   = fermi.eps(i,iSpin,ik);
   PsiG psi1 = psi1Start; // copy

   double   lambda = 1., lambdaMin, s0;
   double   p, q, pTrial, qTrial, f0, f1;
   double   r2, eps2;

   SxString lhsFile, eFile, rFile;
   if (write)  {
      lhsFile = SxString("lhs")+"-i-"+i+"-ik-"+ik+".dat";
      eFile = SxString("eFunc")+"-i-"+i+"-ik-"+ik+".dat";
      rFile = SxString("res")+"-i-"+i+"-ik-"+ik+".dat";

      // --- clean up possible previous files produced here
      SX_OUTPUT_FILE (lhsFile);
      SX_OUTPUT_FILE (eFile);
      SX_OUTPUT_FILE (rFile);

      // --- check Sternheimer equation
      if (checkSternheimerEq)  {
         sxprintf ("Start:  ");
         checkSternheimer (psi1, h1psi0);
         sxprintf ("\n");
      }  else  {
         checkSternheimer (psi1, h1psi0, false);
      }
   }

   // --- convergence criterion
   if (epsResidue > 0.)  eps2 = epsResidue * epsResidue;
   else                  eps2 = -1.;

   // --- optimised steepest descent scheme
   int it;
   for (it = 0; it < maxSternSteps; it++)  {

      // --- gradient at point |psi1>
      PsiG h0eps0psi1 = H * psi1 - eps0 * psi1;
      PsiG gOrtho = h0eps0psi1 + h1psi0;

      // --- Only search in concave part of Hilbert space, i.e., among the
      //     occupied states!
      waves.setOrthogonal (&gOrtho, nValStates, iSpin, ik,SxPW::DONT_NORMALIZE);

      // --- convergence?
      r2 = gOrtho.normSqr ();
      if (r2 < eps2 && it >= minSternSteps)  break;

      // --- negative slope and search direction at |psi1>
      s0        = - sqrt (r2);
      PsiG searchDir = gOrtho / s0;  // normalised to length 1

      // --- a further point |psi_i^1>^(trial) for line minimization
      PsiG psi1Trial = psi1 + lambda * searchDir;

      // --- 1/2 < phi^1_i | (H^0 - eps^0_i) | phi^1_i > =: p
      p      = 0.5 * dot(psi1, h0eps0psi1).re;
      pTrial = 0.5 * dot(psi1Trial, H * psi1Trial - eps0 * psi1Trial).re;

      // --- Re {< phi_i^1 | H1 | phi_i^0 >} =: q
      q         =   dot(psi1, h1psi0).re;
      qTrial    =   dot(psi1Trial, h1psi0).re;

      // --- functional: values at |psi1> and |psi1Trial>
      f0        =   p + q;
      f1        =   pTrial + qTrial;

      double lhsVal = -1.;
      if (write)  {
         // --- variables to be written to output files
         // |LHS| of Sthmr. eq.
         lhsVal = checkSternheimer (psi1, h1psi0, false);
      }

      // --- minimum of line minimisation
      lambdaMin =   0.5 *  s0 * lambda * lambda / (s0 * lambda + f0 - f1);
      psi1.plus_assign_ax (lambdaMin, searchDir);

      // --- project |psi1> onto concave part of Hilbert space
      waves.setOrthogonal (&psi1, nValStates, iSpin, ik, SxPW::DONT_NORMALIZE);

      if (write)  {
         // --- append values to convergence files
         SxFileIO::appendToFile
         (   SxString(it)                + "\t"  // iteration
           + SxString(lhsVal, "%15.12f") + "\n"  // value
         , lhsFile);
         SxFileIO::appendToFile
         (   SxString(it)                + "\t"  // iteration
           + SxString(sqrt (r2), "%15.12f")   + "\n"  // value
         , rFile);
         SxFileIO::appendToFile
         (   SxString(it)                + "\t"  // iteration
           + SxString(f0, "%15.12f")     + "\n"  // value
         , eFile);
      }
   }

   if (write) {
      if (checkSternheimerEq)  {
         sxprintf ("    ->  ");
         checkSternheimer (psi1, h1psi0);
         sxprintf ("   (after %d iter.)\n", it);
      }
   } else {
      checkSternheimer (psi1, h1psi0);
      cout << endl;  cout.flush ();
   }

   if (it >= maxSternSteps-1)
      sxprintf ("WARNING: maximum number of steps reached!\n");

   return psi1;
}


PsiG
SxExxSolver::sternheimerConjGrad_write (const PsiRef &h1psi0,
                                        const PsiRef &psi1Start,
                                        bool write)
{
   SX_CHECK (dynamic_cast<SxPW *>(wavesPtr.getPtr ()));
   SxPW &waves = dynamic_cast<SxPW &>(*wavesPtr);
#ifndef NDEBUG
   SxGkBasis &Gk = waves.getGkBasis();
#endif
   SX_CHECK (dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ()));
   SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian &>(*hamPtr);
   SX_CLOCK      (Timer::Sternheimer);
   SX_CHECK (h1psi0.getSize() == psi1Start.getSize(),
             h1psi0.getSize(), psi1Start.getSize());
   SX_CHECK (h1psi0.auxData.iSpin == psi1Start.auxData.iSpin,
             (int)h1psi0.auxData.iSpin, (int)psi1Start.auxData.iSpin);
   SX_CHECK (h1psi0.auxData.ik == psi1Start.auxData.ik,
             h1psi0.auxData.ik, psi1Start.auxData.ik);
   SX_CHECK (h1psi0.getBasisPtr() == &(Gk(h1psi0.auxData.ik)),
             (size_t)h1psi0.getBasisPtr(),
             (size_t)&(Gk(h1psi0.auxData.ik)));
   SX_CHECK (h1psi0.getBasisPtr() == psi1Start.getBasisPtr(),
             (size_t)h1psi0.getBasisPtr(),
             (size_t)psi1Start.getBasisPtr());

   int     i      = h1psi0.auxData.iOther[0];
   int     iSpin  = h1psi0.auxData.iSpin;
   int     ik     = h1psi0.auxData.ik;

   double  eps0   = fermi.eps(i,iSpin,ik);
   PsiG psi1 = psi1Start; // copy

   double  alpha, beta;
   double  r2, r2old, eps2;
   int     it;

   SxString lhsFile, eFile, rFile;
   if (write)  {
      lhsFile = SxString("lhs")+"-i-"+i+"-ik-"+ik+".dat";
      eFile = SxString("eFunc")+"-i-"+i+"-ik-"+ik+".dat";
      rFile = SxString("res")+"-i-"+i+"-ik-"+ik+".dat";

      // --- clean up possible previous files produced here
      SX_OUTPUT_FILE (lhsFile);
      SX_OUTPUT_FILE (eFile);
      SX_OUTPUT_FILE (rFile);
   }

   // --- check Sternheimer equation
   double lhsVal = -1.;
   if (checkSternheimerEq)  {
      sxprintf ("Start:  ");
      lhsVal = checkSternheimer (psi1, h1psi0);
      sxprintf ("\n");
   }  else if (write) {
      lhsVal = checkSternheimer (psi1, h1psi0, false);
   }

   // --- convergence criterion
   if (epsResidue > 0.)  eps2 = epsResidue * epsResidue;
   else                  eps2 = -1.;

   // --- (negative) gradient
   PsiG h0eps0psi1 = H * psi1 - eps0 * psi1;
   PsiG r = -h1psi0 - h0eps0psi1;

   // --- project gradient onto "convex" part of Hilbert space
   waves.setOrthogonal (&r, nValStates, iSpin, ik, SxPW::DONT_NORMALIZE, 0.l);
   PsiG searchDir = r;
   r2 = r2old = r.normSqr ();

   if (write)  {
      // --- value of the Sternheimerian energy functional
      double eVal = sternFunctional (psi1, h1psi0);

      // --- append values to convergence files
      SxFileIO::appendToFile
      (   SxString(0)                 + "\t"  // iteration
        + SxString(lhsVal, "%15.12f") + "\n"  // value
      , lhsFile);
      SxFileIO::appendToFile
      (   SxString(0)                 + "\t"  // iteration
        + SxString(sqrt (r2), "%15.12f")   + "\n"  // value
      , rFile);
      SxFileIO::appendToFile
      (   SxString(0)                 + "\t"  // iteration
        + SxString(eVal, "%15.12f")   + "\n"  // value
      , eFile);
   }

   // --- Conjugate(d) Gradient scheme
   for (it = 1; it <= maxSternSteps; it++)  {

      // --- convergence?
      if (r2 < eps2 && it >= minSternSteps)  break;

      PsiG h0eps0searchDir = H * searchDir - eps0 * searchDir;

      alpha        = r2 / dot(searchDir, h0eps0searchDir).re;

      // --- improved psi1 in "convex" part of Hilbert space
      psi1.plus_assign_ax (alpha, searchDir);
      waves.setOrthogonal (&psi1, nValStates, iSpin, ik,
                           SxPW::DONT_NORMALIZE, 0.l);  // TODO goes to end?

      // --- new gradient in "convex" part of Hilbert space
      r           -= alpha * h0eps0searchDir;
      waves.setOrthogonal (&r, nValStates, iSpin, ik,
                           SxPW::DONT_NORMALIZE, 0.l);
      r2 = r.normSqr ();

      if (write)  {
         // --- variables to be written to output files
         lhsVal = checkSternheimer (psi1, h1psi0, false);  // |LHS| of Sthmr. eq.
         double eVal = sternFunctional (psi1, h1psi0);

         // --- append values to convergence files
         SxFileIO::appendToFile
         (   SxString(it)                + "\t"  // iteration
           + SxString(lhsVal, "%15.12f") + "\n"  // value
         , lhsFile);
         SxFileIO::appendToFile
         (   SxString(it)                  + "\t"  // iteration
           + SxString(sqrt(r2), "%15.12f") + "\n"  // value
         , rFile);
         SxFileIO::appendToFile
         (   SxString(it)                + "\t"  // iteration
           + SxString(eVal, "%15.12f")   + "\n"  // value
         , eFile);
      }

      // --- get new search direction
      beta         = r2 / r2old;
      searchDir    = r + beta * searchDir;

      r2old = r2;
   }

   if (checkSternheimerEq)  {
      sxprintf ("    ->  ");
      checkSternheimer (psi1, h1psi0);
      sxprintf ("   (after %d iter.)\n", it-1);
   }

   if (it >= maxSternSteps)
      sxprintf ("WARNING: maximum number of steps reached!\n");

   return psi1;
}

PsiG
SxExxSolver::sternheimerPrecCG_write (const PsiRef &h1psi0,
                                      const PsiRef &psi1Start,
                                      bool write)
{
   //SVN_HEAD;
   SX_CHECK (dynamic_cast<SxPW *>(wavesPtr.getPtr ()));
   SxPW &waves = dynamic_cast<SxPW &>(*wavesPtr);
#ifndef NDEBUG
   SxGkBasis &Gk = waves.getGkBasis();
#endif
   SX_CHECK (dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ()));
   SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian &>(*hamPtr);
   SX_CLOCK      (Timer::Sternheimer);
   SX_CHECK (h1psi0.getSize() == psi1Start.getSize(),
             h1psi0.getSize(), psi1Start.getSize());
   SX_CHECK (h1psi0.auxData.iSpin == psi1Start.auxData.iSpin,
             (int)h1psi0.auxData.iSpin, (int)psi1Start.auxData.iSpin);
   SX_CHECK (h1psi0.auxData.ik == psi1Start.auxData.ik,
             h1psi0.auxData.ik, psi1Start.auxData.ik);
   SX_CHECK (h1psi0.getBasisPtr() == &(Gk(h1psi0.auxData.ik)),
             (size_t)h1psi0.getBasisPtr(),
             (size_t)&(Gk(h1psi0.auxData.ik)));
   SX_CHECK (h1psi0.getBasisPtr() == psi1Start.getBasisPtr(),
             (size_t)h1psi0.getBasisPtr(),
             (size_t)psi1Start.getBasisPtr());

   int         i      = h1psi0.auxData.iOther[0];
   int         iSpin  = h1psi0.auxData.iSpin;
   int         ik     = h1psi0.auxData.ik;
              
   double      eps0   = fermi.eps(i,iSpin,ik);
   PsiG psi1 = psi1Start; // copy

   SxComplex16 alpha, beta;
   PrecCoeffG  tr;
   double      r2, eps2;
   int         it;

   SxString lhsFile, eFile, rFile;
   if (write)  {
      lhsFile = SxString("lhs")+"-i-"+i+"-ik-"+ik+".dat";
      eFile = SxString("eFunc")+"-i-"+i+"-ik-"+ik+".dat";
      rFile = SxString("res")+"-i-"+i+"-ik-"+ik+".dat";

      // --- clean up possible previous files produced here
      SX_OUTPUT_FILE (lhsFile);
      SX_OUTPUT_FILE (eFile);
      SX_OUTPUT_FILE (rFile);
   }

   // --- check Sternheimer equation
   double lhsVal = -1.;
   if (checkSternheimerEq)  {
      sxprintf ("Start:  ");
      lhsVal = checkSternheimer (psi1, h1psi0);
      sxprintf ("\n");
   }  else if (write) {
      lhsVal = checkSternheimer (psi1, h1psi0, false);
   }

   // --- convergence criterion
   if (epsResidue > 0.)  eps2 = epsResidue * epsResidue;
   else                  eps2 = -1.;

   // --- (negative) gradient
   PsiG h0eps0psi1 = H * psi1 - eps0 * psi1;
   PsiG r = -h1psi0 - h0eps0psi1;

   // --- project gradient onto "convex" part of Hilbert space
   waves.setOrthogonal (&r, nValStates, iSpin, ik, SxPW::DONT_NORMALIZE, 0.l);
   r2 = r.normSqr ();
   if (write)  {
      // --- value of the Sternheimerian energy functional
      double eVal = sternFunctional (psi1, h1psi0);

      // --- append values to convergence files
      SxFileIO::appendToFile
      (   SxString(0)                 + "\t"  // iteration
        + SxString(lhsVal, "%15.12f") + "\n"  // value
      , lhsFile);
      SxFileIO::appendToFile
      (   SxString(0)                 + "\t"  // iteration
        + SxString(sqrt(r2), "%15.12f")   + "\n"  // value
      , rFile);
      SxFileIO::appendToFile
      (   SxString(0)                 + "\t"  // iteration
        + SxString(eVal, "%15.12f")   + "\n"  // value
      , eFile);
   }
 
   // --- setup preconditioner
   PsiRef psiI = waves(i,iSpin,ik);
   SxVector<double> preconditioner;
   if (precondType == 0)  {
      preconditioner = H.preconditioner (psiI, SxPWHamiltonian::Payne);
   }  else  {
      preconditioner = H.preconditioner (psiI, SxPWHamiltonian::Arias);
   }

   // --- get search direction
   PsiG K = preconditioner * r;
   PsiG searchDir = K;  // = d
   waves.setOrthogonal (&searchDir, nValStates, iSpin, ik,
                        SxPW::DONT_NORMALIZE, 0.l);

   tr = dot(r,K);  // = r M^-1 r

   // --- Preconditioned Conjugate(d) Gradient scheme
   for (it = 1; it <= maxSternSteps; it++)  {

      // --- convergence?
      if (r2 < eps2 && it >= minSternSteps)  break;

      // = A * d
      PsiG h0eps0searchDir  = H * searchDir - eps0 * searchDir;

      alpha        = tr / dot(searchDir, h0eps0searchDir);

      // --- improved psi1 in "convex" part of Hilbert space
      psi1.plus_assign_ax (alpha, searchDir);
      waves.setOrthogonal (&psi1, nValStates, iSpin, ik,
                           SxPW::DONT_NORMALIZE, 0.l);  // TODO: goes to end?

      // --- new gradient in "convex" part of Hilbert space
      r           -= alpha * h0eps0searchDir;   // recursion
      waves.setOrthogonal (&r, nValStates, iSpin, ik,
                           SxPW::DONT_NORMALIZE, 0.l);
      r2 = r.normSqr ();
      if (write)  {
         // --- variables to be written to output files
         lhsVal = checkSternheimer (psi1, h1psi0, false);  // |LHS| of Sthmr. eq.
         double eVal = sternFunctional (psi1, h1psi0);

         // --- append values to convergence files
         SxFileIO::appendToFile
         (   SxString(it)                + "\t"  // iteration
           + SxString(lhsVal, "%30.27f") + "\n"  // value
         , lhsFile);
         SxFileIO::appendToFile
         (   SxString(it)                + "\t"  // iteration
           + SxString(sqrt(r2), "%60.57f")   + "\n"  // value
         , rFile);
         SxFileIO::appendToFile
         (   SxString(it)                + "\t"  // iteration
           + SxString(eVal, "%15.12f")   + "\n"  // value
         , eFile);
      }

      // --- get new search direction
      K            = preconditioner * r;
      tr           = dot(r, K);  // = r M^-1 r
      beta         = -dot(K, h0eps0searchDir) / dot(searchDir, h0eps0searchDir);

      searchDir    = K + beta * searchDir;
      waves.setOrthogonal (&searchDir, nValStates, iSpin, ik,
                           SxPW::DONT_NORMALIZE, 0.l);
   }

   if (checkSternheimerEq)  {
      sxprintf ("    ->  ");
      checkSternheimer (psi1, h1psi0);
      sxprintf ("   (after %d iter.)\n", it-1);
   }

   if (it >= maxSternSteps)
      sxprintf ("WARNING: maximum number of steps reached!\n");

   return psi1;
}

double SxExxSolver::checkSternheimer (const PsiRef &psi1,
                                      const PsiRef &h1psi0,
                                      bool  dump)
{
   SX_CHECK (dynamic_cast<SxPW *>(wavesPtr.getPtr ()));
   SxPW &waves = dynamic_cast<SxPW &>(*wavesPtr);
   SX_CHECK (dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ()));
   SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian &>(*hamPtr);
   int     i      = h1psi0.auxData.iOther[0];
   int     iSpin  = psi1.auxData.iSpin;
   int     ik     = psi1.auxData.ik;
   double eps0 = fermi.eps(i, iSpin, ik);

   PsiG null = H * psi1 - eps0 * psi1 + h1psi0;

   waves.setOrthogonal (&null, nValStates, iSpin, ik,
                        SxPW::DONT_NORMALIZE, 0.l);

   double absNull    = null.norm ();

   if (dump)  sxprintf ("Sternheimer test: %15.12f", absNull);

   return absNull;
}

double SxExxSolver::sternFunctional (const PsiRef &psi1, const PsiRef &h1psi0)
{
   SX_CHECK (dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ()));
   SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian &>(*hamPtr);
   int     i      = h1psi0.auxData.iOther[0];
   int     iSpin  = psi1.auxData.iSpin;
   int     ik     = psi1.auxData.ik;
   double  eps0   = fermi.eps(i,iSpin,ik);

   PsiG    h0eps0psi1 = H * psi1 - eps0 * psi1;
   double  p, q;

   // --- 1/2 < phi^1_i | (H^0 - eps^0_i) | phi^1_i > =: p
   p = 0.5 * dot(psi1, h0eps0psi1).re;

   // --- Re {< phi_i^1 | H1 | phi_i^0 >} =: q
   q = dot(psi1, h1psi0).re;

   return p + q;
}

PsiG SxExxSolver::greensfunction (const PsiRef &h1psi0)
{
   cout << "Muh!" << endl;
   SX_CHECK(wavesPtr);
   SxGkBasis &Gk = wavesPtr->getGkBasis();

   //SVN_HEAD;
   SX_CLOCK      (Timer::GreensFunction);
   SX_CHECK (h1psi0.getBasisPtr() == &(Gk(h1psi0.auxData.ik)),
             (size_t)h1psi0.getBasisPtr(),
             (size_t)&(Gk(h1psi0.auxData.ik)));

   int          i      = h1psi0.auxData.iOther[0];
   int          iSpin  = h1psi0.auxData.iSpin;  // TODO: ugly
   int          ik     = h1psi0.auxData.ik;     // TODO: ugly
   double       epsV   = fermi.eps(i,iSpin,ik);
   double       epsC;
   int          nc     = fermi.getNConductionBands(0,ik);
   PsiG psi1(Gk(ik));
   int          cIdx;
   SxComplex16  elem;
   SX_CHECK (wavesPtr);
   SxPWSet &waves = *wavesPtr;

   // --- initialise |psi1>
   psi1.set (0.);
   psi1.auxData.iOther[0] = i;
   psi1.auxData.iSpin     = char(iSpin);
   psi1.auxData.ik        = ik;

   // --- Greens function
   for (int ic = 0; ic < nc; ic++)  {
//      checkSternheimer (psi1, h1psi0);

      cIdx = fermi.getConductionBandIdx (ic,iSpin,ik);
      PsiRef psiC = waves(cIdx,iSpin,ik);
      epsC = fermi.eps(cIdx,iSpin,ik);

      // < psi0 | H1 | psi1>
      elem = dot(psiC, h1psi0);

      // G H1 | psi1 >
      psi1.plus_assign_ax (elem / (epsV - epsC), psiC);
   }

   checkSternheimer (psi1, h1psi0);

   return  psi1;
}

SxMeshG SxExxSolver::computeRhoInducedVogl (const SxMeshG &vExG,
                                            bool writeConvergence)
{
   SX_CHECK (dynamic_cast<SxPW *>(wavesPtr.getPtr ()));
   SxPW &waves = dynamic_cast<SxPW &>(*wavesPtr);
   SxGkBasis &Gk = waves.getGkBasis();
   SX_CHECK (dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ()));
   SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian &>(*hamPtr);
   SX_CHECK (dynamic_cast<SxExx*>(H.xcPtr.getPtr ()));
   SxExx &exx = dynamic_cast<SxExx&>(*H.xcPtr);
   cout << "SxExxSolver::computeRhoInducedVogl ... starting" << endl;

   //SX_CHECK (H.fastENl == true, H.fastENl);

   int iSpin = 0;  // spin polarisation not yet implemented
   int ik, nk = waves.getNk();
   int mu, muIdx;

   PsiG h1psi0G, h1psi0Gamp;
   PsiG psi1R, psi1Ramp;
   PsiR vExR = ( R | vExG );
   PsiR resR(R);
   resR.set (0.);

   PrecWeights kWeight;

   double fftFactorRtoG = 1. / R.fft3d.scaleFor;  // = sqrt (omega)
   double omega = structure.cell.volume;          // = omega
   double prefactorLocal    = 2. * omega;
   double prefactorNonlocal = 2. * omega * fftFactorRtoG;

   if (!fockStorage)  exx.eX = 0.;

   for (ik = 0; ik < nk; ik++)  {
      kWeight = Gk.weights(ik);

      if (!fockStorage)  {
         if (useNonsymFock)  exx.fockPtr->compute (ik, waves, fermi);
         else                exx.symFockPtr->compute2 (ik, waves, fermi);
      }

      for (mu = 0; mu < nValStates; mu++)  {
         muIdx = fermi.getValenceBandIdx(mu,iSpin,ik);
         PsiRef psi0G = waves(muIdx,iSpin,ik);
         PsiR psi0R = ( R | psi0G );

         if (fockStorage && useNonsymFock)  {
            h1psi0G    = 2. * prefactorLocal    * ( Gk(ik) | (vExR * psi0R ) )
                            - prefactorNonlocal * ( *(exx.fock(ik)) * psi0G );

            h1psi0Gamp =    - prefactorNonlocal * ( *(exx.fock(ik)) * psi0G );
         }  else if (!fockStorage && useNonsymFock)  {
            h1psi0G    = 2. * prefactorLocal    * ( Gk(ik) | (vExR * psi0R ) )
                            - prefactorNonlocal * ( *(exx.fockPtr) * psi0G );

            h1psi0Gamp =    - prefactorNonlocal * ( *(exx.fockPtr) * psi0G );

            // --- exchange contribution to total energy
            exx.eX += kWeight * dot((*(exx.fockPtr) * psi0G), psi0G).re;
            exx.eXc = exx.eExchange = exx.eX;
         }  else if (fockStorage && !useNonsymFock)  {
            h1psi0G    = 2. * prefactorLocal    * ( Gk(ik) | (vExR * psi0R) )
                            - prefactorNonlocal * ( *(exx.symFock(ik)) * psi0G );

            h1psi0Gamp =    - prefactorNonlocal * ( *(exx.symFock(ik)) * psi0G );
         }  else if (!fockStorage && !useNonsymFock)  {
            h1psi0G    = 2. * prefactorLocal    * ( Gk(ik) | (vExR * psi0R) )
                            - prefactorNonlocal * ( *(exx.symFockPtr) * psi0G );

            h1psi0Gamp =    - prefactorNonlocal * ( *(exx.symFockPtr) * psi0G );

            // --- exchange contribution to total energy
            exx.eX += kWeight
                    * dot((*(exx.symFockPtr) * psi0G), psi0G).re;
            exx.eXc = exx.eExchange = exx.eX;
         }  else  {
            sxprintf ("Error: Code should never get here.\n");
            SX_QUIT;
         }

         // TODO: ugly
         h1psi0G.auxData.iOther[0] = muIdx;
         h1psi0G.auxData.iSpin     = char(iSpin);
         h1psi0G.auxData.ik        = ik;
         // TODO: ugly
         h1psi0Gamp.auxData.iOther[0] = muIdx;
         h1psi0Gamp.auxData.iSpin     = char(iSpin);
         h1psi0Gamp.auxData.ik        = ik;

         if (useGreensfunction)  {
            waves1(mu,iSpin,ik)    = greensfunction (h1psi0G);
            waves1amp(mu,iSpin,ik) = greensfunction (h1psi0Gamp);
         }  else  {
cout << "|waves1 - vorher| = " << waves1(mu,iSpin,ik).norm () << endl;
cout << "|waves1amp - vorher| = " << waves1amp(mu,iSpin,ik).norm () << endl;

           if (!writeConvergence)  {
              switch (sternMethod)  {
                 case StDesc   :  waves1(mu,iSpin,ik)
                                  <<= sternheimerStDesc (h1psi0G,
                                                          waves1(mu,iSpin,ik));
                                  waves1amp(mu,iSpin,ik)
                                  <<= sternheimerStDesc (h1psi0Gamp,
                                                       waves1amp(mu,iSpin,ik));
                                  break;
                 case ConjGrad :  waves1(mu,iSpin,ik)
                                  <<= sternheimerConjGrad (h1psi0G,
                                                          waves1(mu,iSpin,ik));
                                  waves1amp(mu,iSpin,ik)
                                  <<= sternheimerConjGrad (h1psi0Gamp,
                                                       waves1amp(mu,iSpin,ik));
                                  break;
                 case PrecCG   :  waves1(mu,iSpin,ik)
                                  <<= sternheimerPrecCG (h1psi0G,
                                                         waves1(mu,iSpin,ik));
                                  waves1amp(mu,iSpin,ik)
                                     <<= sternheimerPrecCG (h1psi0Gamp,
                                                       waves1amp(mu,iSpin,ik));
                                  break;
                 default       :  sxprintf ("Error: Code should never go "
                                            "here.\n");
                                  SX_QUIT;
              }
           }  else  {  // writeConvergence == true
              switch (sternMethod)  {
                 case StDesc   :  waves1(mu,iSpin,ik)
                                  <<= sternheimerStDesc_write (h1psi0G,
                                                      waves1(mu,iSpin,ik), mu);
                                  waves1amp(mu,iSpin,ik)
                                  <<= sternheimerStDesc (h1psi0Gamp,
                                                       waves1amp(mu,iSpin,ik));
                                  break;
                 case ConjGrad :  waves1(mu,iSpin,ik)
                                  <<= sternheimerConjGrad_write (h1psi0G,
                                                       waves1(mu,iSpin,ik), mu);
                                  waves1amp(mu,iSpin,ik)
                                  <<= sternheimerConjGrad (h1psi0Gamp,
                                                       waves1amp(mu,iSpin,ik));
                                  break;
                 case PrecCG   :  waves1(mu,iSpin,ik)
                                  <<= sternheimerPrecCG_write (h1psi0G,
                                                       waves1(mu,iSpin,ik), mu);
                                  waves1amp(mu,iSpin,ik)
                                     <<= sternheimerPrecCG (h1psi0Gamp,
                                                       waves1amp(mu,iSpin,ik));
                                  break;
                 default       :  sxprintf ("Error: Code should never go "
                                            "here.\n");
                                  SX_QUIT;
              }
           }  // :if (!writeConvergence)

cout << "|waves1 - nachher| = " << waves1(mu,iSpin,ik).norm () << endl;
cout << "|waves1amp - nachher| = " << waves1amp(mu,iSpin,ik).norm () << endl;
         }

         psi1R    = ( R | waves1(mu,iSpin,ik) );
         psi1Ramp = ( R | waves1amp(mu,iSpin,ik) );

         resR += kWeight * ( psi0R.conj() * psi1R + psi0R * psi1Ramp.conj() );
      }  // :mu
   }  // :ik

   // --- check the correct norm of rho
   //double f = sqrt((double)(R.getMeshSize())) / fftFactorRtoG;
   //cout << "f = " << f << endl;
   //psi0R /= f;
   //cout << "|psi0R| = " << psi0R.norm () << endl;  // should be 1
   //resR /= (f * f * prefactorLocal);
   //cout << "|resR| = " << resR.norm () << endl;    // should be the
   //SX_EXIT;                                        // correct norm of rho

   if (!useAbdullahSymmetry)  resR = R.symmetrize (resR);

   cout << "SxExxSolver::computeRhoInducedVogl ... done" << endl;

   return ( G | resR );
}

SxMeshG SxExxSolver::computeRhoInduced (const SxMeshG &vExG,
                                        bool writeConvergence)
{
   SX_CHECK (dynamic_cast<SxPW *>(wavesPtr.getPtr ()));
   SxPW &waves = dynamic_cast<SxPW &>(*wavesPtr);
   SxGkBasis &Gk = waves.getGkBasis();
   SX_CHECK (dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ()));
   SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian &>(*hamPtr);
   SX_CHECK (dynamic_cast<SxExx*>(H.xcPtr.getPtr ()));
   SxExx &exx = dynamic_cast<SxExx&>(*H.xcPtr);
   cout << "SxExxSolver::computeRhoInduced ... starting" << endl;

   //SX_CHECK (H.fastENl == true, H.fastENl);

   int iSpin = 0;  // spin polarisation not yet implemented
   int ik, nk = waves.getNk();
   int mu, muIdx;

   PsiG h1psi0G;
   PsiR psi1R;
   PsiR vExR = ( R | vExG );
   PsiR resR(0., R);

   PrecWeights kWeight;

   double fftFactorRtoG = 1. / R.fft3d.scaleFor;  // = sqrt (omega)
   double omega = structure.cell.volume;          // = omega
   double prefactorLocal    = 2. * omega;
   double prefactorNonlocal = 2. * omega * fftFactorRtoG;

   if (!fockStorage)  exx.eX = 0.;

   for (ik = 0; ik < nk; ik++)  {
      kWeight = Gk.weights(ik);

      if (!fockStorage)  {
         if (useNonsymFock)  exx.fockPtr->compute (ik, waves, fermi);
         else                exx.symFockPtr->compute2 (ik, waves, fermi);
      }

      for (mu = 0; mu < nValStates; mu++)  {
         muIdx = fermi.getValenceBandIdx(mu,iSpin,ik);
         PsiRef psi0G = waves(muIdx,iSpin,ik);
         PsiR psi0R = ( R | psi0G );

         if (fockStorage && useNonsymFock)  {
            h1psi0G =   prefactorLocal    * ( Gk(ik) | (vExR * psi0R) )
                      - prefactorNonlocal * ( *(exx.fock(ik)) * psi0G );
         }  else if (!fockStorage && useNonsymFock)  {
            h1psi0G =   prefactorLocal    * ( Gk(ik) | (vExR * psi0R) )
                      - prefactorNonlocal * ( *(exx.fockPtr) * psi0G );

            // --- exchange contribution to total energy
            exx.eX  += kWeight * dot(*(exx.fockPtr) * psi0G, psi0G).re;
            exx.eXc  = exx.eExchange = exx.eX;
         }  else if (fockStorage && !useNonsymFock)  {
            h1psi0G =   prefactorLocal    * ( Gk(ik) | (vExR * psi0R) )
                      - prefactorNonlocal * ( *(exx.symFock(ik)) * psi0G );
         }  else if (!fockStorage && !useNonsymFock)  {
            h1psi0G =   prefactorLocal    * ( Gk(ik) | (vExR * psi0R) )
                      - prefactorNonlocal * ( *(exx.symFockPtr) * psi0G );

            // --- exchange contribution to total energy
            exx.eX  += kWeight * dot((*(exx.symFockPtr) * psi0G), psi0G).re;
            exx.eXc  = exx.eExchange = exx.eX;
         }  else  {
            sxprintf ("Error: Code should never get here.\n");
            SX_QUIT;
         }

         // TODO: ugly
         h1psi0G.auxData.iOther[0] = muIdx;
         h1psi0G.auxData.iSpin     = char(iSpin);
         h1psi0G.auxData.ik        = ik;

         if (useGreensfunction)  {
            waves1(mu,iSpin,ik) = greensfunction (h1psi0G);
         }  else  {
            if (!writeConvergence)  {
               switch (sternMethod)  {
                  case StDesc   :  waves1(mu,iSpin,ik)
                                   <<= sternheimerStDesc (h1psi0G,
                                                          waves1(mu,iSpin,ik));
                                   break;
                  case ConjGrad :  waves1(mu,iSpin,ik)
                                   <<= sternheimerConjGrad (h1psi0G,
                                                          waves1(mu,iSpin,ik));
                                   break;
                  case PrecCG   :  waves1(mu,iSpin,ik)
                                   <<= sternheimerPrecCG (h1psi0G,
                                                          waves1(mu,iSpin,ik));
                                   break;
                  default       :  sxprintf ("Error: Code should never go "
                                             "here.\n");
                                   SX_QUIT;
               }
            }  else  {  // writeConvergence == true
               switch (sternMethod)  {
                  case StDesc   :  waves1(mu,iSpin,ik)
                                   <<= sternheimerStDesc_write (h1psi0G,
                                                       waves1(mu,iSpin,ik), mu);
                                   break;
                  case ConjGrad :  waves1(mu,iSpin,ik)
                                   <<= sternheimerConjGrad_write (h1psi0G,
                                                       waves1(mu,iSpin,ik), mu);
                                   break;
                  case PrecCG   :  waves1(mu,iSpin,ik)
                                   <<= sternheimerPrecCG_write (h1psi0G,
                                                       waves1(mu,iSpin,ik), mu);
                                   break;
                  default       :  sxprintf ("Error: Code should never go "
                                             "here.\n");
                                   SX_QUIT;
               }
            }  // :if (!writeConvergence)
         }  // :if (useGreensfunction)

         psi1R  = ( R | waves1(mu,iSpin,ik) );

         resR  += kWeight * ( psi0R.conj() * psi1R + psi0R * psi1R.conj() );
      }  // :mu
   }  // :ik

   if (!useAbdullahSymmetry)  resR = R.symmetrize (resR);

   cout << "SxExxSolver::computeRhoInduced ... done" << endl;

   return ( G | resR );
}

// --- dump plot along the (1,1,1)-direction in case of fcc - for debugging
void SxExxSolver::writeR111Plot (const SxString &filename,
                                 const SxVector<double> &vecR)
{
   SX_CHECK (vecR.getBasisPtr() == &R,
             (size_t)vecR.getBasisPtr(), (size_t)&R);
   SX_CHECK (vecR.getSize() == R.getMeshSize(),
             vecR.getSize(), R.getMeshSize());

   double volSqrt = sqrt (structure.cell.volume);
   SxVector<SxComplex16> vecG;
   SxVector<double>      vecR111;

   vecG    = (G | vecR) * volSqrt;
   vecR111 = G.GtoR111 (vecG, 75, 25, structure.cell);

   SxTextIO (filename).writeXYPlot (vecR111);
}

// --- dump vector in G basis - for debugging
void SxExxSolver::writeGPlot (const SxString &filename,
                              const SxVector<double> &vecR)
{
   SX_CHECK (vecR.getBasisPtr() == &R,
             (size_t)vecR.getBasisPtr(), (size_t)&R);
   SX_CHECK (vecR.getSize() == R.getMeshSize(),
             vecR.getSize(), R.getMeshSize());

   double volSqrt = sqrt (structure.cell.volume);
   writeGPlot (filename, (G | vecR) * volSqrt);
}

void SxExxSolver::writeGPlot (const SxString &filename,
                              const SxVector<SxComplex16> &vecG)
{
//   SX_CHECK (vecG.getBasisPtr() == &G, (int)vecG.getBasisPtr(), (int)&G);
//   SX_CHECK (vecG.getSize() == G.ng, vecG.getSize(), G.ng);

   int  g, ng = (int)vecG.getSize();
   SxVector<double> vecGRe(ng), vecGIm(ng);

   for (g = 0; g < ng; g++)  {
      vecGRe(g) = vecG(g).re;
      vecGIm(g) = vecG(g).im;
   }

   SxTextIO(filename + "-re.dat").writeXYPlot (vecGRe);
   SxTextIO(filename + "-im.dat").writeXYPlot (vecGIm);
}

void SxExxSolver::relaxEXXPotentialLin (const SxSymbolTable *cmd, int step)
{
   SX_CHECK (dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ()));
   SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian &>(*hamPtr);
   SX_CHECK (dynamic_cast<SxExx*>(H.xcPtr.getPtr ()));
   SxExx &exx = dynamic_cast<SxExx&>(*H.xcPtr);
   cout << "SxExxSolver::relaxEXXPotentialLin ... starting" << endl;

   SX_CHECK (cmd);

   // --- get input parameters
   int             maxRelaxXSteps=-1;  // Don't confuse with member "maxSteps"!
   double          alpha=-1.;
   bool            writeRelaxPlots = false;
   bool            writeInnerSteps = false;
   bool            writeFirstInnerSteps = false;
   SxSymbolTable  *relax = NULL;

   try  {
      relax                =  cmd->getGroup("relaxXPotential");
      maxRelaxXSteps       =  relax->get("maxRelaxXSteps")->toInt();
      alpha                =  relax->get("xMixing")->toReal();
      writeInnerSteps      = (relax->contains("writeInnerSteps"))
                           ?  relax->get("writeInnerSteps")->toAttribute()
                           :  false;
      writeFirstInnerSteps = (relax->contains("writeFirstInnerSteps"))
                           ?  relax->get("writeFirstInnerSteps")->toAttribute()
                           :  false;
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

   // --- dump parameters
   cout << SX_SEPARATOR;
   sxprintf ("| Relax potential w.r.t. given wavefunctions\n");
   cout << SX_SEPARATOR;
   sxprintf ("|   number of steps:   %d\n", maxRelaxXSteps);
   sxprintf ("|   mixing parameter:  %g\n", alpha);
   if (writeInnerSteps)  {
      sxprintf ("|\n");
      sxprintf ("|   Dumps plot of every step.\n");
   }
   if (writeFirstInnerSteps)  {
      sxprintf ("|\n");
      sxprintf ("|   Dumps plot of every step within the 1st outer loop.\n");
   }

   // --- validate parameters
   if (writeInnerSteps && writeFirstInnerSteps)  {
      cout << SX_SEPARATOR;
      sxprintf ("| Error: Either use flag \"writeInnerSteps\" or "
                "\"writeFirstInnerSteps\"\n"
                "|        or none of them. Both at a time is not meaningful. "
                "Sorry.\n");
      SX_QUIT;
   }

   SxMeshG     rhoIndG;
   SxMeshG     vExNewG;
   SX_FFT_REAL fftFactorGtoR = R.fft3d.scaleFor;  // = 1 / sqrt(Omega)

   // --- write results of this loop to output file(s)?
   if (writeInnerSteps)                    writeRelaxPlots = true;
   if (writeFirstInnerSteps && step == 1)  writeRelaxPlots = true;

   // --- the inner loop
   for (int it = 1; it <= maxRelaxXSteps; it++)  {
      cout << SX_SEPARATOR;
      sxprintf ("|   Iteration: outer loop -> %d, inner loop -> %d\n",
                step-1, it);

      // --- get induced charge density
      if (mimicVoglChi)
         rhoIndG = computeRhoInducedVogl (exx.vExLoopG(0));
      else
         SX_EXIT;

      // --- linear mixing
      vExNewG = exx.vExLoopG(0) + alpha * rhoIndG;

      // --- cut high frequencies
      exx.vExLoopG(0).set(0.);
      exx.vExLoopG(0)( SxIdx(1,exx.nGChi-1) ) <<= vExNewG( SxIdx(1,exx.nGChi-1) );

      // --- write potential along (1,1,1) direction
      if (writeRelaxPlots)  {
         exx.vEx(0) = ( R | exx.vExLoopG(0) ).real() * fftFactorGtoR;
         exx.vEx(0).setBasis (&R);
         exx.vEx(0) = R.symmetrize (exx.vEx(0));
         writeR111Plot (SxString("vEx-")+(step-1)+"-"+it+".dat", exx.vEx(0));
      }
   }

   // --- get symmetrised EXX potential
   exx.vEx(0) = ( R | exx.vExLoopG(0) ).real() * fftFactorGtoR;
   exx.vEx(0).setBasis (&R);
   exx.vEx(0) = R.symmetrize (exx.vEx(0));

   exx.vXc(0) <<= exx.vEx(0);

   cout << "SxExxSolver::relaxEXXPotentialLin ... done" << endl;
}

void SxExxSolver::relaxEXXPotential (const SxSymbolTable *cmd, int step)
{
   SX_CLOCK (Timer::RelaxEXX);
   SX_CHECK (dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ()));
   SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian &>(*hamPtr);
   SX_CHECK (dynamic_cast<SxExx*>(H.xcPtr.getPtr ()));
   SxExx &exx = dynamic_cast<SxExx&>(*H.xcPtr);
   cout << "SxExxSolver::relaxEXXPotential ... starting" << endl;

   SX_CHECK (cmd);

   // --- get input parameters
   SxSymbolTable          *relax = NULL;
   int                     maxRelaxXSteps;  // Don't confuse with 'maxSteps'!
   int                     minRelaxXSteps;
   double                  epsX;
   double                  alpha;
   int                     mixingMethod, nPulaySteps = -1;
   bool                    useKerker = false;
   SxString                mixingStr;
   SxRhoMixer::MixerType   mixerType = SxRhoMixer::Pulay;
   double                  kerkerScaling = -1., kerkerDamping = -1.,
                           spinMixing = -1.;
   bool                    writeRelaxPlots = false;
   bool                    writeInnerSteps = false;
   bool                    writeFirstInnerSteps = false;

   try  {
      relax                =  cmd->getGroup("relaxXPotential");
      maxRelaxXSteps       =  relax->get("maxRelaxXSteps")->toInt();
      minRelaxXSteps       = (relax->contains("minRelaxXSteps"))
                           ?  relax->get("minRelaxXSteps")->toInt()
                           :  -1;
      epsX                 = (relax->contains("epsX"))
                           ?  relax->get("epsX")->toReal()
                           :  -1.;
      alpha                =  relax->get("xMixing")->toReal();
      mixingMethod         =  relax->get("mixingMethod")->toInt();
      spinMixing           =  1.;  // spin polarisation not yet implemented
      writeInnerSteps      = (relax->contains("writeInnerSteps"))
                           ?  relax->get("writeInnerSteps")->toAttribute()
                           :  false;
      writeFirstInnerSteps = (relax->contains("writeFirstInnerSteps"))
                           ?  relax->get("writeFirstInnerSteps")->toAttribute()
                           :  false;

      switch (mixingMethod)  {
         case 0 :  mixerType = SxRhoMixer::Linear;
                   useKerker = false;
                   mixingStr = "linear mixer";
                   nPulaySteps = 1;
                   break;
         case 1 :  mixerType = SxRhoMixer::Linear;
                   useKerker = true;
                   mixingStr = "precond. linear mixer";
                   nPulaySteps = 1;
                   kerkerScaling = relax->get("kerkerScaling")->toReal();
                   kerkerDamping = relax->get("kerkerDamping")->toReal();
                   break;
         case 2 :  mixerType = SxRhoMixer::Pulay;
                   useKerker = false;
                   mixingStr = "Pulay mixer";
                   nPulaySteps = relax->get("nPulaySteps")->toInt();
                   break;
         case 3 :  mixerType = SxRhoMixer::Pulay;
                   useKerker = true;
                   mixingStr = "precond. Pulay mixer";
                   nPulaySteps = relax->get("nPulaySteps")->toInt();
                   kerkerScaling = relax->get("kerkerScaling")->toReal();
                   kerkerDamping = relax->get("kerkerDamping")->toReal();
                   break;
         default : sxprintf ("| Unknown mixing method: %d. Sorry.\n",
                             mixingMethod);
                   SX_QUIT;
      }
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

   // --- dump input parameters
   cout << SX_SEPARATOR;
   sxprintf ("| Relax potential w.r.t. given wavefunctions\n");
   cout << SX_SEPARATOR;
   sxprintf ("|   max. number of steps:   %d\n", maxRelaxXSteps);
   if (minRelaxXSteps > 0)  {
      sxprintf ("|   min. number of steps:   %d\n", minRelaxXSteps);
   }  else  {
      sxprintf ("|   max. number of steps:   not given\n");
   }
   if (epsX > 0.)  {
      sxprintf ("|   epsX:                   %g\n", epsX);
   }  else  {
      sxprintf ("|   epsX:                   not given\n");
   }
   if (writeInnerSteps)  {
      sxprintf ("|\n");
      sxprintf ("|   Dumps plot of every step.\n");
   }
   if (writeFirstInnerSteps)  {
      sxprintf ("|\n");
      sxprintf ("|   Dumps plot of every step within the 1st outer loop.\n");
   }
   sxprintf ("|\n");
   sxprintf ("|   Mixer:\n");
   sxprintf ("|     mixing scheme:      %s\n", mixingStr.ascii());
   sxprintf ("|     mixing parameter:   %g\n", alpha);
   if (nSpin == 2)  {  // Can't happen currently, of course.
      sxprintf ("|     spin mixing:        %g\n", spinMixing);
   }
   if (mixerType == SxRhoMixer::Pulay)  {
      sxprintf ("|     Pulay steps:        %d\n", nPulaySteps);
   }
   if (useKerker)  {
      sxprintf ("|     Kerker scaling:     %g\n", kerkerScaling);
      sxprintf ("|     Kerker damping:     %g\n", kerkerDamping);
   }
   cout.flush();

   // --- validate parameters
   if (writeInnerSteps && writeFirstInnerSteps)  {
      cout << SX_SEPARATOR;
      sxprintf ("| Error: Either use flag \"writeInnerSteps\" or "
                "\"writeFirstInnerSteps\"\n"
                "|        or none of them. Both at a time is not meaningful. "
                "Sorry.\n");
      SX_QUIT;
   }

   if (minRelaxXSteps > maxRelaxXSteps)  {
      cout << SX_SEPARATOR;
      sxprintf ("| Error:  minRelaxXSteps = %d > %d = maxRelaxXSteps\n",
                minRelaxXSteps, maxRelaxXSteps);
      SX_QUIT;
   }

   SxMeshG     rhoIndG;
   SxMeshG     vExNewG;
   SX_FFT_REAL fftFactorGtoR = R.fft3d.scaleFor;  // = 1 / sqrt(Omega)
   double      rhoIndAbs;
   VxcR        vExLoopR(nSpin);
   RhoR        rhoIndR(nSpin);
   RhoG        rhoIndGSpin(nSpin);
   VxcR        vExNewR(nSpin);
   SxString    rhoIndFile = "rhoInd.dat";
   SxString    dVExLoopFile = "dVExLoop.dat";
   SxMeshG     vExLoopGOld, dVExLoopG;
   double      dVExLoopGNorm;
   SX_FFT_REAL rhoNormFac  = 0.5 * fftFactorGtoR
                                 / sqrt ((double)R.getMeshSize());
   SX_FFT_REAL dVExNormFac = 1.0 / (fftFactorGtoR
                                 * sqrt ((double)R.getMeshSize()));

   // --- write results of this loop to output file(s)?
   if (writeInnerSteps)                    writeRelaxPlots = true;
   if (writeFirstInnerSteps && step == 1)  writeRelaxPlots = true;

   // --- mark outer loop
   SxFileIO::appendToFile (SxString("# --- outer loop -> ")+step+"\n", rhoIndFile);
   SxFileIO::appendToFile (SxString("# --- outer loop -> ")+step+"\n", dVExLoopFile);

   // --- set up mixer
   
   // --- mixing variant 1a
   /*
   SxRhoMixer mixer (mixerType, alpha, nPulaySteps, spinMixing);
   if (useKerker)
      mixer.preconditioner.useKerker (kerkerScaling, kerkerDamping);
   else {
      mixer.rhoMixing = 1.;
      mixer.preconditioner.scaling = alpha;
   }
   mixer.setNormModus (SxRhoMixer::RenormOff);
   */
   

   // --- mixing variant 2
   SxRhoMixerG mixer ((SxRhoMixerG::MixerType)(int)mixerType,
                      alpha, nPulaySteps, spinMixing);
   mixer.setNormModus (SxRhoMixerG::RenormOff);
   if (useKerker)  {
      SxRhoMixerG::SxKerker Kerker(G, kerkerScaling, kerkerDamping);
      mixer.setPreconditioner (Kerker);
   }

   //SxPulay<SxVector<PrecCoeffG>> pulayG;
   SxVector<PrecCoeffG> rhoIndGOld;

   // --- the inner loop
   for (int it = 1; it <= maxRelaxXSteps; it++)  {
      cout << SX_SEPARATOR;
      sxprintf ("|   Iteration: outer loop -> %d, inner loop -> %d\n",
                step-1, it);

      // --- get induced charge density
      if (mimicVoglChi)
         rhoIndG = computeRhoInducedVogl (exx.vExLoopG(0), (step==1 && it==1));
      else
         rhoIndG = computeRhoInduced (exx.vExLoopG(0), (step==1 && it==1));

      // --- mixing

      // --- mixing variant 1: R-mixer
      /*
      vExLoopR(0) = ( R | exx.vExLoopG(0) );
      rhoIndR(0)  = ( R | rhoIndG );

      // --- mixing variant 1a: general mixer in R
      mixer.addRhoIn (vExLoopR);
      mixer.addResidue (rhoIndR);
      vExNewR = mixer.getMixedRho ();

      // --- mixing variant 1b: linear mixer in R
cout << "Mixing: linear mit Fourier ..." << endl;
      vExNewR(0) = vExLoopR(0) + alpha * rhoIndR(0);

      // mixing variant 1
      vExNewG = ( G | vExNewR(0) );
      */

      // --- mixing variant 2: G mixer
cout << "Mixing: SxRhoMixerG ohne Fourier ..." << endl;
      rhoIndGSpin(0) = rhoIndG * 1.;
      mixer.addRhoInG (exx.vExLoopG);
      mixer.addResidueG (rhoIndGSpin);
      vExNewG = mixer.getMixedRhoG()(0);
      // mixing variant 2b: linear in G space
      //vExNewG = exx.vExLoopG(0) + alpha * rhoIndG;
      /*
      if (it > 1)  {
         pulayG.addStepRes (exx.vExLoopG(0)-vExLoopGOld, rhoIndG - rhoIndGOld);
         pulayG.compute (rhoIndG);
         vExNewG = pulayG.makeStep (exx.vExLoopG(0))
                 + alpha * pulayG.getMixedResidue (rhoIndG);
         if (pulayG.steps.getSize () >= nPulaySteps)  {
            pulayG.steps.removeFirst ();
            pulayG.residues.removeFirst ();
         }
      } else {
         // linear mixing
         vExNewG = exx.vExLoopG(0) + alpha * rhoIndG;
      }
      rhoIndGOld = rhoIndG;
      */

      // --- store old exchange potential
      vExLoopGOld = exx.vExLoopG(0) * 1.;

      // --- cut high frequencies
      exx.vExLoopG(0).set(0.);
      exx.vExLoopG(0)( SxIdx(1,exx.nGChi-1) ) <<= vExNewG( SxIdx(1,exx.nGChi-1) );

      // --- write potential along (1,1,1) direction
      if (writeRelaxPlots)  {
         exx.vEx(0) = ( R | exx.vExLoopG(0) ).real() * fftFactorGtoR;
         exx.vEx(0).setBasis (&R);
         exx.vEx(0) = R.symmetrize (exx.vEx(0));
         writeR111Plot (SxString("vEx-")+(step-1)+"-"+it+".dat", exx.vEx(0));
      }

      // --- write to rhoInd convergence file
      rhoIndAbs = rhoNormFac * rhoIndG.norm ();
      SxFileIO::appendToFile (SxString(rhoIndAbs) + "\n", rhoIndFile);

      // --- compute norm of the exchange potential
      dVExLoopG     = exx.vExLoopG(0) - vExLoopGOld;
      dVExLoopGNorm = dVExNormFac * dVExLoopG.norm ();
      SxFileIO::appendToFile (SxString(dVExLoopGNorm) + "\n", dVExLoopFile);
   }

   // --- get symmetrised EXX potential
   exx.vEx(0) = ( R | exx.vExLoopG(0) ).real() * fftFactorGtoR;
   exx.vEx(0).setBasis (&R);
   exx.vEx(0) = R.symmetrize (exx.vEx(0));

   exx.vXc(0) <<= exx.vEx(0);

   cout << "SxExxSolver::relaxEXXPotential ... done" << endl;
}

void SxExxSolver::updateEXXCPotential (const SxSymbolTable *cmd, int step)
{
   SX_CHECK (dynamic_cast<SxPW *>(wavesPtr.getPtr ()));
   SxPW &waves = dynamic_cast<SxPW &>(*wavesPtr);
   SX_CHECK (dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ()));
   SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian &>(*hamPtr);
   SX_CHECK (dynamic_cast<SxExx*>(H.xcPtr.getPtr ()));
   SxExx &exx = dynamic_cast<SxExx&>(*H.xcPtr);
   cout << "SxExxSolver::updateEXXCPotential ... starting" << endl;

   SxXC::XCFunctional xcFunctional = exx.xcFunctional;

   bool calcXo = exx.getExchangeStatus ();
   // --- calculate x and c within the specified approximation
   switch (xcFunctional)  {
      case SxXC::EXX        : exx.vCor(0).resize (R.getMeshSize());
                              exx.vCor(0).set (0.);
                              exx.eCorrelation = exx.eC = 0.;
                              break;
      case SxXC::EXX_LDA    : exx.rhoXcR(0) = H.rho(0);
                              exx.disableExchange ();
                              exx.computeLDA ();
                              break;
      case SxXC::EXX_PBE    : SX_EXIT; // vEx/vCor not set up
                              //exx.updateGGA_PBE_WB (H.rho(0));
                              break;
      case SxXC::EXX_PBE_WB : SX_EXIT; // vEx/vCor not set up
                              //exx.updateGGA_PBE_WB (H.rho(0));
                              break;
      default : sxprintf ("| You must choose one of the following "
                          "exchange-correlation functionals:\n"
                          "|\n"
                          "|        EXX, EXX_LDA, EXX_PBE or EXX_PBE_WB\n");
                SX_QUIT;
   }
   if (calcXo) exx.enableExchange ();

   // --- overwrite vEx and vXc by EXX calculation
//   relaxEXXPotentialLin (cmd, step);
   relaxEXXPotential (cmd, step);
   if (fockStorage)  {
      if (useNonsymFock)  exx.computeFockEnergy (waves, fermi);
      else                exx.computeSymFockEnergy (waves, fermi);
   }

   // --- add correlation part of vXc
   exx.vXc(0) += exx.vCor(0);

   // --- add correlation part of eXc
   switch (xcFunctional)  {
      case SxXC::EXX        : break;  // nothing to be done
      case SxXC::EXX_LDA    : exx.eXc += exx.eC;
                              break;
      case SxXC::EXX_PBE    : exx.eXc += exx.eCorrelation;
                              break;
      case SxXC::EXX_PBE_WB : exx.eXc += exx.eCorrelation;
                              break;
      default : sxprintf ("| Error: The program never should get here.\n");
                SX_EXIT;
   }

   // --- set (G=0)-component to zero
   exx.shiftXC ();

   // --- compute the exchange/correlation potential energy
   exx.computeXCPotEnergy (H.rho.rhoR);

   cout << "SxExxSolver::updateEXXCPotential ... done" << endl;
}

void SxExxSolver::gsEXX (const SxSymbolTable *cmd, bool calc)
{
   //SVN_HEAD;
   SX_CHECK (dynamic_cast<SxPW *>(wavesPtr.getPtr ()));
   SxPW &waves = dynamic_cast<SxPW &>(*wavesPtr);
   SxGkBasis &Gk = waves.getGkBasis();
   SX_CHECK (dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ()));
   SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian &>(*hamPtr);
   SX_CHECK (dynamic_cast<SxExx*>(H.xcPtr.getPtr ()));
   SxExx &exx = dynamic_cast<SxExx&>(*H.xcPtr);
   SX_CLOCK (Timer::ExxSolver);
   SX_CHECK (cmd);

   cout << SX_SEPARATOR;
   cout << "| Exact Exchange Formalism -- Groundstate Based Method (gsEXX)\n";
   cout << SX_SEPARATOR;

   // --- get input parameters
   SX_START_TIMER (Timer::ExxStartup);
   int             maxStepsEXX;  // These parameters belong to the outer loop.
   PrecEnergy      dEnergyEXX;   // "maxSteps" and "dEnergy" are controlled by
                                 // the inner loop (group "relaxRho").
   SxSymbolTable  *relaxRho = NULL, *stern = NULL, *direct = NULL;
   bool            withInterimResults;

   try  {
      maxStepsEXX           = (cmd->contains("maxSteps"))
                            ?  cmd->get("maxSteps")->toInt()
                            :  1000000;  // approximately unlimited
      dEnergyEXX            =  cmd->get("dEnergy")->toReal();
      relaxRho              =  cmd->getGroup("relaxRho");
      useGreensfunction     = (cmd->contains("useGreensfunction"))
                            ?  cmd->get("useGreensfunction")->toAttribute()
                            :  false;
      mimicVoglChi          = (cmd->contains("mimicVoglChi"))
                            ?  cmd->get("mimicVoglChi")->toAttribute()
                            :  false;
      useAbdullahSymmetry   = (cmd->contains("useAbdullahSymmetry"))
                            ?  cmd->get("useAbdullahSymmetry")->toAttribute()
                            :  false;
      useNonsymFock         = (cmd->contains("useNonsymFock"))
                            ?  cmd->get("useNonsymFock")->toAttribute()
                            :  true;  // to be changed into 'false'!!
      fockStorage           = (cmd->contains("fockStorage"))
                            ?  cmd->get("fockStorage")->toAttribute()
                            :  true;
      withInterimResults    = (cmd->contains("withInterimResults"))
                            ?  cmd->get("withInterimResults")->toAttribute()
                            :  false;
      stern                 = (cmd->containsGroup("Sternheimer"))
                            ?  cmd->getGroup("Sternheimer")
                            :  NULL;
      direct                = (cmd->containsGroup("fullDiag"))
                            ?  cmd->getGroup("fullDiag")
                            :  NULL;
      if (stern)  {
         minSternSteps      = (stern->contains("minSternSteps"))
                            ?  stern->get("minSternSteps")->toInt()
                            :  -1;
         maxSternSteps      = (stern->contains("maxSternSteps"))
                            ?  stern->get("maxSternSteps")->toInt()
                            :  1000000;  // approximately unlimited
         epsResidue         = (stern->contains("epsResidue"))
                            ?  stern->get("epsResidue")->toReal()
                            :  -1.;
         sternMethod        = (stern->contains("method"))
                            ? (SternMethod)(stern->get("method")->toInt())
                            :  NoMethod;
         precondType        = (stern->contains("preconditioner"))
                            ?  stern->get("preconditioner")->toInt()
                            :  -1;
         checkSternheimerEq = (stern->contains("checkEquation"))
                            ?  stern->get("checkEquation")->toAttribute()
                            :  false;
      }  else  {
         minSternSteps = maxSternSteps = precondType = -1;
         epsResidue = -1.;
         checkSternheimerEq = false;
      }
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

   // --- validate input parameters
   if (dEnergyEXX <= 0.)  {
      sxprintf ("| Input parameter \"dEnergy\" should be greater than 0. "
                "Please change the\n"
                "| input file and try again!\n");
      SX_QUIT;
   }

   if (nStatesChi < nStates && useGreensfunction)  {
      sxprintf ("| Error: \"nEmptyStatesChi\" must be no smaller than "
                "\"nEmptyStates\".\n"
                "         Please increase it in your input file and try "
                "again!\n");
      SX_QUIT;
   }

   if (nSpin > 1)  {
      sxprintf ("| Spin polarization not yet implemented. Sorry.\n");
      SX_QUIT;
   }

   if (useGreensfunction && !direct)  {
      sxprintf ("| If you want to use the Greensfunction, you must provide\n"
                "| direct diagonal diagonalisation. Add a group \"fullDiag\"\n"
                "| into the \"gsEXX\" group and try again!.\n");
      SX_QUIT;
   }

   if (!useGreensfunction && direct)  {
      sxprintf ("| Group \"fullDiag\" found within the \"gsEXX\" group. "
                "Maybe, you want\n"
                "| to use the Greensfunctions instead of the Sternheimer "
                "equation? If so,\n"
                "| please add the flag \"useGreensfunction\", or else remove "
                "the group \"fullDiag\".\n");
      SX_QUIT;
   }

   if (useGreensfunction && stern)  {
      sxprintf ("| Either you use the Sternheimer equation or the "
                "Greensfunction,\n"
                "| both at the same time is not meaningful. Please remove "
                "either the\n"
                "| flag \"useGreensfunction\" or the group \"fullDiag\" - but "
                "not both of\n"
                "| them :-) - and try again, please!\n");
      SX_QUIT;
   }

   if (!useGreensfunction && !stern)  {
      sxprintf ("| You have to use either the Sternheimer equation or the "
                "Greensfunction\n"
                "| in this formalism. Add the group \"Sternheimer\" or the "
                "flag\n"
                "| \"useGreensfunction\" into the \"gsEXX\" group and try "
                "again, please!\n");
      SX_QUIT;
   }

   if (stern)  {
      if (maxSternSteps >= 1000000 && epsResidue < 0.)  {
         sxprintf ("| Error: Group \"Sternheimer\" must contain "
                   "\"maxSternSteps\" or\n"
                   "|        \"epsResidue\" or both. Sorry.\n");
         SX_QUIT;
      }

      if (maxSternSteps < minSternSteps)  {
         sxprintf ("| Error: \"maxSternSteps\" should be greater than "
                   "\"minSternSteps\",\n"
                   "|        or what's your opinion, huh?\n");
         SX_QUIT;
      }

      if ((int)sternMethod < 1 || (int)sternMethod > 3)  {
         sxprintf ("| Error: Not a valid minimisation method chosen within "
                   "\"Sternheimer\"\n"
                   "|        group (method = %d). Sorry.\n",
                   (int)sternMethod);
         SX_QUIT;
      }

      if (sternMethod == PrecCG)  {
         if (precondType < 0 || precondType > 1)  {
            sxprintf ("| Error: Not a valid preconditioner chosen within "
                      "\"Sternheimer\"\n"
                      "|        group. Possible choices are \"PAYNE\" and "
                      "\"ARIAS\".\n");
            SX_QUIT;
         }
      }

      if (sternMethod != PrecCG && precondType != -1)  {
         sxprintf ("| Error: In \"Sternheimer\" group: Preconditioning only "
                   "works for\n"
                   "|        Preconditioned Conjugate Gradient method. "
                   "Sorry.\n");
         SX_QUIT;
      }
   }

   // --- dump parameters
   {
      const SxSymbolTable *top = cmd->topLevel ();
      int      nEmptyStates    = SxHamiltonian::getNEmptyStates (top);
      int      nEmptyStatesChi = SxHamiltonian::getNEmptyStatesChi (top);
      nValStates               = fermi.getNValenceBands(0,0);
      double   eCut            = G.getECut (top);
      double   eCutChi         = top->getGroup("basis")
                                    ->get("eCutChi")->toReal();
      sxprintf ("|   SCF electronic minimisation\n");
      sxprintf ("|\n");
      sxprintf ("|      max. number of steps:      %d\n", maxStepsEXX);
      sxprintf ("|      dEnergy (conv. criterion): %g\n", dEnergyEXX);
      sxprintf ("|\n");
      sxprintf ("|      eCut:                      %g Ry\n", eCut);
      sxprintf ("|      eCutChi:                   %g Ry\n", eCutChi);
      sxprintf ("|\n");
      sxprintf ("|      nStates:                   %d\n", nStates);
      sxprintf ("|      nEmptyStates:              %d\n", nEmptyStates);
      sxprintf ("|      nStatesChi:                %d\n", nStatesChi);
      sxprintf ("|      nEmptyStatesChi:           %d\n", nEmptyStatesChi);
      sxprintf ("|      number of valence bands:   %d\n", nValStates);
      sxprintf ("|\n");
      sxprintf ("|      useAbdullahSymmetry:       ");
      if (useAbdullahSymmetry)  sxprintf ("yes\n");
      else                      sxprintf ("no\n");
      sxprintf ("|      Fock operator symmetrised: ");
      if (useNonsymFock)        sxprintf ("no\n");
      else                      sxprintf ("yes\n");
      sxprintf ("|      mimicVoglChi:              ");
      if (mimicVoglChi)         sxprintf ("yes\n");
      else                      sxprintf ("no\n");
      sxprintf ("|      useGreensfunction:         ");
      if (useGreensfunction)    sxprintf ("yes\n");
      else                      sxprintf ("no\n");
      sxprintf ("|      fockStorage:               ");
      if (fockStorage)          sxprintf ("yes\n");
      else                      sxprintf ("no\n");
      sxprintf ("|\n");
      if (stern)  {
         sxprintf ("|   Solving Sternheimer equation:\n");
         sxprintf ("|\n");
         sxprintf ("|      minSternSteps:             ");
         if (minSternSteps > -1)  sxprintf ("%d\n", minSternSteps);
         else                     sxprintf ("not given\n");
         sxprintf ("|      maxSternSteps:             ");
         if (maxSternSteps > -1)  sxprintf ("%d\n", maxSternSteps);
         else                     sxprintf ("not given\n");
         sxprintf ("|      epsResidue:                ");
         if (epsResidue > -1.)    sxprintf ("%g\n", epsResidue);
         else                     sxprintf ("not given\n");
         sxprintf ("|      method:                    ");
         switch (sternMethod)  {
            case StDesc   :  sxprintf ("steepest descent\n");
                             break;
            case ConjGrad :  sxprintf ("conjugate gradient\n");
                             break;
            case PrecCG   :  sxprintf ("precond. conjugate gradient\n");
                             break;
            default       :  sxprintf ("\n\n"
                                       "Error: "
                                       "Program should never get here.\n");
                             SX_QUIT;
         }
         if (sternMethod == PrecCG)  {
            sxprintf ("|      preconditioner:            ");
            switch (precondType)  {
               case 0  :  sxprintf ("Payne\n");
                          break;
               case 1  :  sxprintf ("Arias\n");
                          break;
               default :  sxprintf ("\n\n"
                                    "Error: Program should never get here.\n");
                          SX_QUIT;
            }
         }
         if (checkSternheimerEq)  {
            sxprintf ("|\n");
            sxprintf ("|      The Sternheimer equation will be checked each "
                      "time and a warning\n"
                      "|      will be dumped where required.\n");
         }
      }  else  {
         sxprintf ("|   Not solving Sternheimer equation.\n");
      }
      cout << SX_SEPARATOR;  cout.flush ();
   }

   if (!calc)  return;

   int             it;
   SxMeshG         vExG;
   SxSymbolTable  *relaxRhoCmd = NULL, *directCmd = NULL;
   PrecEnergy      eTot, eTotOld = SX_HUGE;
   typedef         SxPWHamiltonian::Contrib HContrib;
   SX_FFT_REAL     fftFactorRtoG = 1. / R.fft3d.scaleFor;  // = sqrt(Omega)
   int             nStatesOrig = nStates;
   SxString        logStr;
   int             nk = wavesPtr->getNk();

   // --- clean up possible previous files
   SX_OUTPUT_FILE (SxString("rhoInd.dat"));
   SX_OUTPUT_FILE (SxString("dVExLoop.dat"));

   // --- initialise correction functions
   SX_CHECK  (nValStates > 0, nValStates);
   SX_CHECK (nSpin == wavesPtr->getNSpin(), nSpin, wavesPtr->getNSpin());
   waves1 = SxPW (nValStates, nSpin, wavesPtr->getGkBasisPtr());
   waves1.setZero();
   if (mimicVoglChi)  {
      waves1amp = SxPW (nValStates, nSpin, wavesPtr->getGkBasisPtr());
      waves1amp.setZero();
   }

   // --- just for describing procedure in log file
   if (useGreensfunction)  logStr = "still";  else  logStr = "entering";

   // --- resize vExt, used as buffer for vXc
   H.vExtR.resize (nSpin);
   H.vExtR(0).resize (R.getMeshSize());

   // --- initialise exchange (Fock) operator
   int nGFock = exx.nGChi;

   if (fockStorage)  {
      if (useNonsymFock)  {
         exx.fock.resize (nk);
         for (int ik = 0; ik < nk; ik++)
            exx.fock(ik) = SxPtr<SxFockGk>::create (Gk, nGFock, R.cell);
      }  else  {
         exx.symFock.resize (nk);
         const SxSymbolTable *top = cmd->topLevel ();
         for (int ik = 0; ik < nk; ik++)
            exx.symFock(ik)
               = SxPtr<SxSymFockGk>::create (Gk, nGFock, top);
      }
   }  else  {  // no Fock storage
      if (useNonsymFock)  {
         exx.fockPtr = SxPtr<SxFockGk>::create (Gk, nGFock, R.cell);
      }  else  {
         const SxSymbolTable *top = cmd->topLevel ();
         exx.symFockPtr
            = SxPtr<SxSymFockGk>::create (Gk, nGFock, top);
      }
   }

   // --- get starting point for xc-potential
   exx.vXc(0).setBasis (&R);
   exx.shiftXC ();
   exx.vExLoopG(0) = ( G | exx.vXc(0) ) * fftFactorRtoG;

   // --- dump (shifted) input vX along (1,1,1) direction
   try  {
      SxBinIO io (SxString("vXC-0.sxb"), SxBinIO::BINARY_WRITE_ONLY);
      io.writeMesh (exx.vXc(0) * HA2EV, R.cell, R.getMesh());
      io.setMode (SxBinIO::WRITE_DATA);
      io.writeMesh (exx.vXc(0) * HA2EV, R.cell, R.getMesh());
      io.close ();
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }
   writeR111Plot (SxString("vEx-0.dat"), exx.vXc(0));

   // --- use vExtR instead of vXc
   H.contrib = (HContrib)((H.contrib 
                           & ~SxPWHamiltonian::CALC_XC)  // &~ => unset
                          | SxPWHamiltonian::CALC_EXT);  // | => (set)
   H.vExtActsOnNuclei = false;

   SX_STOP_TIMER (Timer::ExxStartup);
   // --- self-consistent (SC) gsEXX loop
   for (it = 1; it <= maxStepsEXX; it++)  {

      // --- (0) provide "complete" basis set, if required
      if (useGreensfunction)  {
         cout << SX_SEPARATOR;
         sxprintf ("| EXX -- entering iteration %d "
                   "-- getting \"complete\" basis\n", it);

         setNStates (nStatesChi);

         for (directCmd  = direct->begin();
              directCmd != NULL;
              directCmd  = directCmd->nextSibling())
         {
            execute (directCmd, calc);
         }
      }
rhoOld = H.rho.rhoR(0) * 1.;

      // --- (1) get xc-potential
      cout << SX_SEPARATOR;
      sxprintf ("| EXX -- %s iteration %d -- updating potentials\n",
                logStr.ascii(), it);

      //H.fastENl = true;         // some minimisers change that to 'false'
      waves.orthonormalize ();  // just to be on the safe side

      // --- calculate exchange (Fock) operator
      if (fockStorage)  {
         if (useNonsymFock)
            for (int ik = 0; ik < nk; ik++)
               exx.fock(ik)->compute (ik, waves, fermi);
         else
            for (int ik = 0; ik < nk; ik++)
               exx.symFock(ik)->compute2 (ik, waves, fermi);
      }

      // --- the formalism's heart is hidden here:
      updateEXXCPotential (cmd, it);

      if (useGreensfunction)  setNStates (nStatesOrig);


      // --- (2) get waves and rho consistent with xc-potential
      cout << SX_SEPARATOR;
      sxprintf ("| EXX -- still iteration %d -- relax waves and rho\n", it);

      H.vExtR(0) <<= exx.vXc(0);

      for (relaxRhoCmd  = relaxRho->begin();
           relaxRhoCmd != NULL;
           relaxRhoCmd  = relaxRhoCmd->nextSibling())
      {
         execute (relaxRhoCmd, calc);
      }

      // --- dump vEx along (1,1,1) direction
      writeR111Plot (SxString("vEx-")+it+".dat", exx.vEx(0));

      // --- total energy
      eTot = H.eTotal;
      sxprintf ("| EXX -- finishing iteration %d -- total energy: %15.12g H\n",
                it, eTot);

      // --- write out interim results
      if (withInterimResults)  {
         cout << SX_SEPARATOR;
         cout << "| Saving interim rho and vXC\n";
         cout << SX_SEPARATOR;

         // --- write rho
         cout << "|   Saving charge density                  ...  ";
         cout.flush ();
         H.rho.writeRho (SxString("rho-")+it+".sxb");  cout << "done\n";

         // --- write vXc
         cout << "|   Saving exchange-correlation potential  ...  ";
         cout.flush ();
         try  {
            SxBinIO io (SxString("vXC-")+it+".sxb", SxBinIO::BINARY_WRITE_ONLY);
            io.writeMesh (exx.vXc(0) * HA2EV, R.cell, R.getMesh());
            io.setMode (SxBinIO::WRITE_DATA);
            io.writeMesh (exx.vXc(0) * HA2EV, R.cell, R.getMesh());
            io.close ();
         }  catch (SxException e)  {
            e.print ();
            SX_QUIT;
         }
         cout << "done\n";
      }

      // --- energy converged?
      if (fabs(eTot - eTotOld) < dEnergyEXX)  {
         cout << SX_SEPARATOR;
         H.printEnergies ();
         fermi.printOccupation ();
         sxprintf ("| EXX self-consistent run: convergence reached.\n");
         break;
      }

      // --- maximum step number exceeded?
      if (it == maxStepsEXX)  {
         cout << SX_SEPARATOR;
         H.printEnergies ();
         fermi.printOccupation ();
         sxprintf ("| WARNING:  EXX self-consistent run: Maximum number of "
                   "steps exceeded.\n"
                   "|           Total energy convergence not yet reached.\n");
         break;
      }

      eTotOld = eTot;
   }  // :it

   // --- write vXc to output file
   cout << SX_SEPARATOR;
   cout << "| Saving exchange-correlation potential  ...  "; cout.flush ();
   try  {
      SxBinIO io ("vXC.sxb", SxBinIO::BINARY_WRITE_ONLY);
      io.writeMesh (exx.vXc(0), R.cell, R.getMesh());
      io.setMode (SxBinIO::WRITE_DATA);
      io.writeMesh (exx.vXc(0), R.cell, R.getMesh());
      io.close ();
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }
   cout << "done" << endl;
}

void SxExxSolver::stEXX (const SxSymbolTable *cmd, bool calc)
{
   SX_CHECK (wavesPtr);
   SxGkBasis &Gk = wavesPtr->getGkBasis();
   SX_CHECK (dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ()));
   SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian &>(*hamPtr);
   SX_CHECK (dynamic_cast<SxExx*>(H.xcPtr.getPtr ()));
   SxExx &exx = dynamic_cast<SxExx&>(*H.xcPtr);
   SX_CLOCK (Timer::ExxSolver);
   SX_CHECK (cmd);

   cout << SX_SEPARATOR;
   cout << "| Exact Exchange Formalism -- Standard Method" << endl;
   cout << SX_SEPARATOR;

   // --- get input parameters
   int             maxStepsEXX;   // These parameters belong to the outer loop.
   PrecEnergy      dEnergyEXX;    // "maxSteps" and "dEnergy" are controlled by
                                  // the inner loop, i.e., "relaxRho".
   bool            relaxation = false;
   SxSymbolTable  *diag = NULL, *relax = NULL;

   try  {
      maxStepsEXX     = (cmd->contains("maxSteps"))
                      ?  cmd->get("maxSteps")->toInt()
                      :  1000000;  // unlimited
      dEnergyEXX      =  cmd->get("dEnergy")->toReal();
      diag            =  cmd->getGroup("fullDiag");
      exx.useCorrectChi = (cmd->contains("useCorrectChi"))
                      ?  cmd->get("useCorrectChi")->toAttribute()
                      :  false;
      if (cmd->containsGroup("relaxRho"))  {
         relax      = cmd->getGroup("relaxRho");
         relaxation = true;
      }
   }  catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   // --- validate parameters
   if (dEnergyEXX <= 0.)  {
      sxprintf ("| Input parameter \"dEnergy\" should be greater than 0. "
                "Change input file\n"
                "| and try again!\n");
      SX_QUIT;
   }

   if (nStatesChi < nStates)  {
      sxprintf ("| Error: \"nEmptyStatesChi\" must be no smaller than "
                "\"nEmptyStates\".\n"
                "       Please increase it in your input file and try "
                "again!\n");
      SX_QUIT;
   }

   if (!relaxation)  {
      sxprintf ("| WARNING:  No relaxation of charge density. This may "
                "result in a slower\n"
                "|           convergence. Think about involving the "
                "optional group \"relaxRho\".\n"
                "|           The final result, however, should not be "
                "influenced.\n");
      cout << SX_SEPARATOR;
   }

   if (nSpin > 1)  {
      sxprintf ("| Spin polarization not yet implemented. Sorry.\n");
      SX_QUIT;
   }

   // --- dump parameters
   sxprintf ("|   max. number of steps:      %d\n", maxStepsEXX);
   sxprintf ("|   dEnergy (conv. criterion): %g\n", dEnergyEXX);
   sxprintf ("|   SCF electronic minimisation\n");
   cout << SX_SEPARATOR;  cout.flush ();

   if (!calc)  return;

   // --- resize vExt, used as buffer for vXc
   H.vExtR.resize(nSpin);
   H.vExtR(0).resize (exx.vXc(0).getSize());

   // --- dump input vXc along (1,1,1) direction
   exx.vXc(0).setBasis (&R);
   exx.shiftXC ();
   writeR111Plot (SxString("vEx-in-shifted.dat"), exx.vXc(0));

   // --- self-consistent (SC) exact exchange loop
   int             it;
   SxSymbolTable  *diagCmd = NULL, *relaxCmd = NULL;
   PrecEnergy      eTot, eTotOld = SX_HUGE;
   typedef         SxPWHamiltonian::Contrib HContrib;
   int             nStatesOrig = nStates;

   // --- initialise exchange (Fock) operator
   int nGChi = exx.nGChi;
//   exx.fockPtr = SxPtr<SxFockG>::create (Gk, nGChi, SxFockG::FCC);
   exx.fockPtr = SxPtr<SxFockGk>::create (Gk, nGChi, R.cell);
//cout << "exx_loop: nGChi = " << nGChi << endl;

   for (it = 0; it < maxStepsEXX; it++)  {

      // --- (1) bandstructure-like run to get full basis set
      setNStates (nStatesChi);
cout << "Gk(0).ng = " << Gk(0).ng << endl;
cout << "nStatesChi = " << nStatesChi << endl;

      cout << SX_SEPARATOR;
      sxprintf ("| EXX -- entering iteration %d "
              "-- getting \"complete\" basis\n", it+1);

      H.contrib = (HContrib)(( H.contrib
                             & ~SxPWHamiltonian::CALC_XC)   // &~ => unset
                             | SxPWHamiltonian::CALC_EXT);  // | => set
      H.vExtActsOnNuclei = false;
      H.vExtR(0) <<= exx.vXc(0);

      nStates = nStatesChi;

      for (diagCmd  = diag->begin();
           diagCmd != NULL;
           diagCmd  = diagCmd->nextSibling())
      {
         execute (diagCmd, calc);
      }

/*
einschleifentest (*(exx.fockPtr));
EXIT;
H.computeVExUnsym (waves, fermi);
H.fastENl = true;
maxSteps = 500;
computeLHSminusRHS (exx.vExUnsym(0), *(exx.fockPtr));
//EXIT;
computeRhoIndSternheimer (*(exx.fockPtr));
EXIT;
SxVector<PrecCoeffG> stoer (G.ng, 1.);
stoer.setBasis (&G);
SxVector<PrecCoeffG> resGreensVor = rhoIndCompare (exx.vExUnsym(0));
SxVector<PrecCoeffG> resSternVor  = rhoIndCompareSternheimer (exx.vExUnsym(0));
cout << SX_SEPARATOR;
cout << "probier (Stern) = " << resSternVor << endl;
cout << "|probier (Stern)| = " << resSternVor.norm () << endl;
SxVector<PrecCoeffG> diffResVor = resGreensVor - resSternVor;
cout << SX_SEPARATOR;
cout << "diffResVor = " << diffResVor << endl;
cout << "|diffResVor| = " << diffResVor.norm () << endl;
*/

//      exx.vXc(0) <<= exx.vExtR(0);
      H.contrib = (HContrib)((H.contrib
                             |  SxPWHamiltonian::CALC_XC)   // | => set
                             & ~SxPWHamiltonian::CALC_EXT); // &~ => unset

      // --- (2) get EXX potential
      cout << SX_SEPARATOR;
      sxprintf ("| EXX -- still iteration %d -- updating potentials\n", it+1);

      // --- The actual "EXX formalism" is solely called, i.e. "hidden",
      //     within the following lines:
      H.update (fermi);

/*
SxVector<PrecCoeffG> resGreens = rhoIndCompare (exx.vExUnsym(0));
//H.contrib = (HContrib)(  H.contrib
//                       ^ SxPWHamiltonian::CALC_XC     // ^ = XOR (-)
//                       | SxPWHamiltonian::CALC_EXT);  // | = OR  (+)
//SxVector<PrecCoeffG> resStern  = rhoIndCompareSternheimer (exx.vExUnsym(0));
SxVector<PrecCoeffG> diffRes = resGreens - resGreensVor;
cout << "diffRes = " << diffRes << endl;
cout << "|diffRes| = " << diffRes.norm () << endl;
EXIT;
*/

      nStates = nStatesOrig;

      setNStates (nStates);

      // --- (3) relax charge density
      if (relaxation)  {
         cout << SX_SEPARATOR;
         sxprintf ("| EXX -- still iteration %d -- relaxing charge density\n",
                 it+1);

         H.contrib = (HContrib)((  H.contrib
                                & ~SxPWHamiltonian::CALC_XC)   // &~ => unset
                                | SxPWHamiltonian::CALC_EXT);  // | => set
         H.vExtActsOnNuclei = false;
         H.vExtR(0) <<= exx.vXc(0);         

         for (relaxCmd  = relax->begin();
              relaxCmd != NULL;
              relaxCmd  = relaxCmd->nextSibling())
         {
            execute (relaxCmd, calc);
         }

         H.contrib = (HContrib)((H.contrib
                                |  SxPWHamiltonian::CALC_XC)    // | => set
                                & ~SxPWHamiltonian::CALC_EXT);  // &~ => unset
      }

      // --- write vXc to output file
      try  {
         SxBinIO io ("vXC.sxb", SxBinIO::BINARY_WRITE_ONLY);
         io.writeMesh (exx.vXc(0), R.cell, R.getMesh());
         io.setMode (SxBinIO::WRITE_DATA);
         io.writeMesh (exx.vXc(0), R.cell, R.getMesh());
         io.close ();
      }  catch (SxException e)  {
         e.print ();
         SX_QUIT;
      }

      // --- dump vEx along (1,1,1) direction
      writeR111Plot (SxString("vEx-")+(it+1)+".dat", exx.vEx(0));

      // --- total energy
      eTot = H.eTotal;
      sxprintf ("| EXX -- finishing iteration %d -- total energy: %15.12g H\n",
              it+1, eTot);

      // --- energy convergence?
      if (fabs(eTot - eTotOld) < dEnergyEXX)  {
         cout << SX_SEPARATOR;
         H.printEnergies ();
         fermi.printOccupation ();
         sxprintf ("| EXX self-consistent run: convergence reached.\n");
         break;
      }

      if (it >= maxStepsEXX-1)  {
         cout << SX_SEPARATOR;
         H.printEnergies ();
         fermi.printOccupation ();
         sxprintf ("| WARNING:  EXX self-consistent run: Maximum number of "
                   "steps exceeded.\n"
                   "|           Convergence not yet reached.\n");
         break;
      }

      eTotOld = eTot;
   }  // :it
}

void SxExxSolver::odp (const SxSymbolTable *cmd, bool calc)
{
   SX_CLOCK (Timer::OrbitalDependentPotential);
   SX_CHECK (wavesPtr);
   SxGkBasis &Gk = wavesPtr->getGkBasis();
   SX_CHECK (cmd);
   SX_CHECK (dynamic_cast<SxPWHamiltonian *>(hamPtr.getPtr ()));
   SxPWHamiltonian &H = dynamic_cast<SxPWHamiltonian &>(*hamPtr);
   SX_CHECK (dynamic_cast<SxExx*>(H.xcPtr.getPtr ()));
   SxExx &exx = dynamic_cast<SxExx&>(*H.xcPtr);

   SxString method;

   cout << SX_SEPARATOR;
   if (H.xcPtr->xcFunctional == SxXC::KLI)  {
      cout << "| KLI -- approximation method to OEP / EXX" << endl;
      method = "KLI";
   }
   else if (H.xcPtr->xcFunctional == SxXC::Slater)  {
      cout << "| Electronic minimisation using  S l a t e r  potential" << endl;
      method = "Slater";
   }
   cout << SX_SEPARATOR;

   // --- get input parameters
   int             maxStepsODP=0; // These parameters belong to the outer loop.
   PrecEnergy      dEnergyODP=0;  // "maxSteps" and "dEnergy" are controlled by
   SxSymbolTable  *relax = NULL;  // the inner loop, i.e., "relaxRho".
   bool            suele = false;

   try  {
      maxStepsODP = (cmd->contains("maxSteps"))
                  ?  cmd->get("maxSteps")->toInt()
                  :  1000000;  // unlimited
      dEnergyODP  =  cmd->get("dEnergy")->toReal();
      relax       =  cmd->getGroup("relaxRho");
      suele       = (cmd->contains("sueleGuess"))
                  ?  cmd->get("sueleGuess")->toAttribute()
                  :  false;
   }  catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

   // --- validate parameters
   if (dEnergyODP <= 0.)  {
      sxprintf ("| Input parameter \"dEnergy\" should be greater than 0. "
                "Change input file\n"
                "| and try again!\n");
      SX_QUIT;
   }

   if (nSpin > 1)  {
      sxprintf ("| Spin polarisation not yet implemented. Sorry.\n");
      SX_QUIT;
   }

   // --- dump parameters
   sxprintf ("|   max. number of steps:      %d\n", maxStepsODP);
   sxprintf ("|   dEnergy (conv. criterion): %g Hartree\n", dEnergyODP);
   if (suele)
      sxprintf ("|   Formular of Suele et al. used for initial guess at vXC."
                "\n");
   sxprintf ("|   SCF electronic minimisation\n");
   cout << SX_SEPARATOR;  cout.flush ();

   if (!calc)  return;

   // --- resize vExt, used here as a buffer for vXc
   H.vExtR.resize (nSpin);
   H.vExtR(0).resize (exx.vXc(0).getSize());

   // --- dump input vXc along (1,1,1) direction
   exx.vXc(0).setBasis (&R);
   writeR111Plot (SxString("vEx-0-in.dat"), exx.vXc(0));

   // --- shift initial potential by setting its average to 0
   exx.shiftXC ();

   // --- initialise non-local exchange (Fock) operator
   int nGFock = exx.nGChi;  // TODO: change the name "eCutChi"
                          //       or provide an alternative
//   exx.fockPtr = SxPtr<SxFockG>::create (Gk, nGFock, SxFockG::FCC);
   exx.fockPtr = SxPtr<SxFockGk>::create (Gk, nGFock, R.cell);
cout << "odp_loop: nGFock = " << nGFock << endl;

   // --- initial guess by mixing with Slater potential,
   //     see J. Chem. Phys. 112, p. 7355 (2000), Eq. (14) by Suele et al.
   if (suele)  {
      VxcR    vSlater (exx.getVSlater (H.rho(0), getWaves(), fermi));
      SxMeshG vG = ( G | vSlater(0) );
      vG(0)      = 0.;
      vSlater(0) = ( R | vG );
      exx.vXc(0) = 0.5 * ( exx.vXc(0) + 2./3. * vSlater(0) );
   }

   // --- dump initial vEx along (1,1,1) direction
   writeR111Plot (SxString("vEx-0-init.dat"), exx.vXc(0));

   // --- self-consistent orbital dependend potential (ODP) loop
   int             it;
   SxSymbolTable  *relaxCmd = NULL;
   PrecEnergy      eTot, eTotOld = SX_HUGE;
   typedef         SxPWHamiltonian::Contrib HContrib;

   for (it = 0; it < maxStepsODP; it++)  {

      // --- (1) update ODP potential
      cout << SX_SEPARATOR;
      sxprintf ("| %s -- entering iteration %d -- updating potential\n",
                method.ascii(), it+1);

      // --- The actual ODP formalism is called by the following line:
      H.update (fermi);

      // --- (2) relax charge density
      cout << SX_SEPARATOR;
      sxprintf ("| %s -- still iteration %d -- relaxing charge density\n",
                method.ascii(), it+1);

      H.contrib = (HContrib)((H.contrib
                             & ~SxPWHamiltonian::CALC_XC)    // => unset
                             |  SxPWHamiltonian::CALC_EXT);  // => set
      H.vExtActsOnNuclei = false;
      H.vExtR(0) <<= exx.vXc(0);

      for (relaxCmd  = relax->begin();
           relaxCmd != NULL;
           relaxCmd  = relaxCmd->nextSibling())
      {
         execute (relaxCmd, calc);
      }

      H.contrib = (HContrib)((H.contrib
                             |  SxPWHamiltonian::CALC_XC)    // => set
                             & ~SxPWHamiltonian::CALC_EXT);  // => unset

      // --- write vXc to output file
      try  {
         SxBinIO io ("vXC.sxb", SxBinIO::BINARY_WRITE_ONLY);
//         SxBinIO io (SxString("vXC-")+(it+1)+".sxb", SxBinIO::BINARY_WRITE_ONLY);
         io.writeMesh (exx.vXc(0), R.cell, R.getMesh());
         io.setMode (SxBinIO::WRITE_DATA);
         io.writeMesh (exx.vXc(0), R.cell, R.getMesh());
         io.close ();
      }  catch (SxException e)  {
         e.print ();
         SX_QUIT;
      }

      // --- dump vEx along (1,1,1) direction
      writeR111Plot (SxString("vEx-")+(it+1)+".dat", exx.vXc(0));

      // --- total energy
      eTot = H.eTotal;
      sxprintf ("| %s -- finishing iteration %d -- total energy: %17.14g H\n",
                method.ascii(), it+1, eTot);

      // --- energy convergence?
      if (fabs(eTot - eTotOld) < dEnergyODP)  {
         cout << SX_SEPARATOR;
         H.printEnergies ();
         fermi.printOccupation ();
         sxprintf ("| %s self-consistent run: convergence reached.\n",
                   method.ascii());
         break;
      }

      if (it >= maxStepsODP-1)  {
         cout << SX_SEPARATOR;
         H.printEnergies ();
         fermi.printOccupation ();
         sxprintf ("| WARNING:  %s self-consistent run: Maximum number of "
                   "steps exceeded.\n"
                   "|           Convergence not yet reached.\n",
                   method.ascii());
         break;
      }

      eTotOld = eTot;
   }  // :it
}

