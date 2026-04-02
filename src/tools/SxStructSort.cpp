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
#include <SxAtomicStructure.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage =
      "This add-on sorts the atoms according to ascending z, y, x.";
   cli.authors = "C. Freysoldt";

   // int binGroup = cli.newGroup ("binary files");
   // standard --sxb option
   bool sxbFile = sxbOption (cli);
   cli.newGroup ("sx files");
   
   // standard --input option
   SxString inFile = inputOption (cli, sxbFile);
   
   // standard output options (-o, --outsxb)
   SxStructOut out (cli, sxbFile);

   cli.finalize ();

   // --- read input
   SxAtomicStructure structure = structureFromInput (inFile, sxbFile);
   // --- sort
   structure.sort ();
   
   // --- output
   out.write (structure);
}

