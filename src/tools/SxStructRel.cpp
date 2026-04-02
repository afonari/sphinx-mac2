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
#include <SxCliStandard.h>
#include <SxConfig.h>
#include <SxAtomicStructure.h>
#include <SxTextIO.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage =
      "This add-on prints a structure in relative coordinates.";
   cli.authors = "C. Freysoldt";

   // standard structure reading options
   bool sxbFile    = sxbOption   (cli);
   SxString inFile = inputOption (cli, sxbFile);
   
   SxString outFile
      = cli.option ("-o","filename", "output file name (screen otherwise)")
        .toString ("");

   cli.finalize ();

   // --- read input
   SxAtomicStructure structure = structureFromInput (inFile, sxbFile);
   
   // --- Relative coordinates
   if (! structure.isPeriodic ())  {
      cout << "Structure from file '" << inFile << "' is not periodic." << endl;
      cout << "Nonperiodic structures do not have relative coordinates.\n";
      SX_QUIT;
   }
   
   // --- output
   if (outFile.getSize () > 0)
      structure.fprintRel (SxTextIO (outFile).getFp ());
   else
      structure.fprintRel (stdout);

}

