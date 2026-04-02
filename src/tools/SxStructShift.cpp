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
#include <SxConfig.h>
#include <SxBinIO.h>
#include <SxCliStandard.h>
#include <SxAtomicStructure.h>
#include <SxStickyFilter.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage =
      "This add-on allows to shift all atoms within the structure.";
   cli.authors = "C. Freysoldt";

   // ---  definition
   SxVector3<double> trans(cli.option ("--by|--vector","vector",
                                       "translation vector").toList3 ());
   
   bool relative = cli.option ("-r|--relative", "assume relative coordinates")
                   .toBool ();

   bool wrap = cli.option ("--wrap", "keep atoms within cell").toBool ();

   // standard sxb option
   int binGroup = cli.newGroup ("binary files");
   bool sxbFile = sxbOption (cli);
   cli.newGroup ("sx files");
   cli.excludeGroup (binGroup);
   bool keepMovable
      = cli.option("-m|--keep-movable","retain mobility of atoms as "
                   "it is defined in the input structure")
        .toBool ();
   cli.setGroup(SxCLI::generalGroup);
   
   // standard input file option
   SxString inFile = inputOption (cli, sxbFile);
   
   // standard output options (-o, --outsxb)
   SxStructOut out (cli, sxbFile);

   cli.finalize ();

   // --- read input
   SxAtomicStructure structure = structureFromInput (inFile, sxbFile, keepMovable);
   
   // --- Relative coordinates
   if (! structure.isPeriodic () && relative)  {
      cout << "Nonperiodic structures do not have relative coordinates.\n";
      cout << "--relative flag can't be used for the structure from." << endl;
      cout << "file '" << inFile << "'." << endl;
      SX_QUIT;
   }
   
   if (relative) structure.cell.changeToCar (&trans);
   
   // shift structure
   structure += trans;

   if (wrap) structure %= structure.cell;

   // --- output
   out.write (structure);

}

