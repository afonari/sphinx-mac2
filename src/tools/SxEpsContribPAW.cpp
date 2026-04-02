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
#include <SxStickyFilter.h>
#include <SxHamSolver.h>
#include <SxPAWHamiltonian.h>
#include <SxPWHamiltonian.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.authors = "B. Lange";
   cli.preUsageMessage =
      "This tool calculates contributions of the energy eigenvalues of the Hamiltonian";
   SxString rhoFile = 
      cli.option ("-r|--rho","file","density file")
      .toString ("rho.sxb");

   SxString inFile =
      cli.option ("-i|--input","file","S/PHI/nX input file")
      .toString ("input.sx");

   SxString wavesFile =
      cli.option ("-w|--waves","file","waves file")
      .toString ("waves.sxb");

   SxFFT::plannerCLI (cli);
   cli.version ("0.1");
   cli.finalize ();

   // --- read input file
   SxParser parser;
   SxParser::Table table = parser.read (inFile);
   
   // The HamSolver
   SxHamSolver hamSolver (wavesFile, rhoFile, table);

   double result = hamSolver.hamPtr->getEnergy(*hamSolver.wavesPtr,hamSolver.fermi);

   cout << SX_SEPARATOR;
   cout << "Etot is " << result << endl;
   cout << SX_SEPARATOR;

   // ---Now specify different Energy contributions
   SxPtr<SxPAWHamiltonian> HPtr = SxPtr<SxPAWHamiltonian> (hamSolver.hamPtr);
   SxPAWHamiltonian &H = *HPtr;
   const SxVecRef<SxComplex16> &Psi = (*hamSolver.wavesPtr)(0,0);
   long int hContrib;
   SxVector<double> energies;
  
   // total energy
   hContrib = SxPAWHamiltonian::CalcAll;
   H.hContrib = hContrib;
   energies = HA2EV * Psi.overlap (H * Psi).diag(); //inefficient
   cout << "Total energies (eV): " << endl;
   energies.print(true);

   // Hartree energy
   hContrib = SxPAWHamiltonian::CalcHartree;
   H.hContrib = hContrib;
   energies = HA2EV * Psi.overlap (H * Psi).diag(); //inefficient
   cout << "Hartree energies (eV): " << endl;
   energies.print(true);

   // Kinetic energy
   hContrib = SxPAWHamiltonian::CalcKin;
   H.hContrib = hContrib;
   energies = HA2EV * Psi.overlap (H * Psi).diag(); //inefficient
   cout << "Kinetic energies (eV): " << endl;
   energies.print(true);
   
   // Nuc energy
   hContrib = SxPAWHamiltonian::CalcNuc;
   H.hContrib = hContrib;
   energies = HA2EV * Psi.overlap (H * Psi).diag(); //inefficient
   cout << "Nuc energies (eV): " << endl;
   energies.print(true);
   
   // Bar energy
   hContrib = SxPAWHamiltonian::CalcBar;
   H.hContrib = hContrib;
   energies = HA2EV * Psi.overlap (H * Psi).diag(); //inefficient
   cout << "Bar energies (eV): " << endl;
   energies.print(true);

   // XC energy
   hContrib = SxPAWHamiltonian::CalcXc;
   H.hContrib = hContrib;
   energies = HA2EV * Psi.overlap (H * Psi).diag(); //inefficient
   cout << "XC energies (eV): " << endl;
   energies.print(true);
}
