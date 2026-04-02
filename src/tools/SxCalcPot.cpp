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

#include <SxCLI.h>
#include <SxBinIO.h>
#include <SxHamSolver.h>
#include <SxPWHamiltonian.h>
#include <SxPAWHamiltonian.h>
#include <SxDipoleCorrZ.h>


int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.authors = "C. Freysoldt";
   cli.preUsageMessage =
      "This add-on calculates the electrostatic, exchange-correlation, and "
      "effective potentials.";
   SxString rhoFile =
      cli.option ("-r|--rho","file","density file")
      .toString ("rho.sxb");

   SxFFT::quickFFTPlanner ();
   SxFFT::plannerCLI (cli);
   SxString inFile =
      cli.option ("-i|--input","file","S/PHI/nX input file")
      .toString ("input.sx");

   SxString structFile =
      cli.option ("-s|--structure","file","S/PHI/nX input file for structure")
      .toString ("");
   cli.last ().defaultValue =
      "default: from input file = no structure optimization";

   bool elStat =
      cli.option ("--elstat", "write electrostatic potential to vElStat.sxb")
      .toBool ();

   bool xc =
      cli.option ("--xc", "write exchange-correlation potential to vXC.sxb")
      .toBool ();

   int xcDebugGroup = cli.newGroup ("Exchange or correlation only");
   bool onlyX =
      cli.option ("--only-exchange", "write exchange potential to vXC.sxb")
      .toBool ();

   bool onlyC =
      cli.option ("--only-correlation",
                  "write correlation potential to vXC.sxb").toBool ();
   cli.newGroup("Effective potential");
   cli.excludeGroup (xcDebugGroup);

   bool eff = cli.option ("--eff", "write effective potential to vEff.sxb")
              .toBool ();
   cli.setGroup (cli.generalGroup);

   bool enforceDipole =
      cli.option ("--dipole","use dipole correction").toBool ();

   cli.version ("2.0");
   cli.finalize ();

   if (onlyX && onlyC)  {
      cout << "Cannot write both only exchange and only correlation." << endl;
      SX_QUIT;
   }


   // --- read input file
   SxParser parser;
   SxParser::Table table = parser.read (inFile);

   SxHamSolver pot;
   // read structure
   if (structFile.getSize () > 0)
      pot.structure = SxAtomicStructure(&*
                      SxParser ().read (structFile,"std/structure.std") );

   pot.init (rhoFile, table);
   SxKPoints kp(pot.structure.cell, &*table);
   double eCut = SxGBasis::getECut (&*table);
   pot.setupHam (rhoFile, table,
                 SxPtr<SxGkBasis>::create(kp, pot.G, eCut, true));

   SxHamiltonian *hamPtr = pot.hamPtr.getPtr ();
   SxPWHamiltonian  *ham = dynamic_cast<SxPWHamiltonian*>(hamPtr);
   SxPAWHamiltonian *paw = dynamic_cast<SxPAWHamiltonian*>(hamPtr);
   if (!ham && !paw) SX_EXIT;

   if (ham) ham->computeESelf ();

   bool recompute = false;
   if (onlyX || onlyC) {
      xc = true;
      recompute = true;
      if (paw)  {
         cout << "PAW only exchange + correlation" << endl;
         SX_QUIT;
      }
      SX_CHECK(ham);
      if (onlyC) ham->xcPtr->disableExchange ();
      if (onlyX) ham->xcPtr->disableCorrelation ();
   }
   // switch on dipole correction if wanted
   if (enforceDipole) {
      if (ham)
         ham->dipoleCorrection = true;
      else if (paw)
         paw->dipoleCorr = SxPtr<SxDipoleCorrZ>::create ();
      recompute = true;
   }

   // recompute minimalistic, if necessary
   if (paw && recompute)  {
      paw->hContrib = 0;
      if (eff || xc)
         paw->hContrib |= SxPAWHamiltonian::CalcXc;
      if (eff || elStat)
         paw->hContrib |= SxPAWHamiltonian::CalcHartree;
      paw->computeRhoTerms ();
   } else if (ham && recompute)  {
      long int contrib = 0;
      if (eff || xc)  {
         if (onlyC)
            contrib = ham->CALC_C;
         else if (onlyX)
            contrib = ham->CALC_X;
         else
            contrib = ham->CALC_XC;
      }
      if (eff || elStat)
         contrib += ham->CALC_HARTREE + ham->CALC_LOC;
      ham->contrib = (SxPWHamiltonian::Contrib)contrib;
      ham->compute (SxFermi (), true);
   }

   // Define generic vXc and vEff
   RhoR &vXc = ham ? ham->xcPtr->vXc
                   : paw->xcPtr->vXc;
   RhoR &vEff = ham ? ham->vEffR
                    : paw->vPS;

   if (eff) SxRho(vEff).writeRho ("vEff.sxb");
   if (xc)  SxRho(vXc) .writeRho ("vXC.sxb");
   if (elStat)  {
      if (recompute && !(xc || eff))
         // vEff = electrostatic only
         SxRho(vEff(0)).writeRho ("vElStat.sxb");
      else
         // vEff = electrostatic + xc
         SxRho(vEff(0) - vXc(0)).writeRho ("vElStat.sxb");
   }
   if (ham)  {
      ham->printHartree = true;
      ham->printEnergies ();
   }

   return 0;
}

