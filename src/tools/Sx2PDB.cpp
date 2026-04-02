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

#include <SxUtil.h>
#include <SxConstants.h>
#include <SxList.h>
#include <SxPrecision.h>
#include <SxVector3.h>
#include <SxMatrix3.h>
#include <SxString.h>
#include <SxPDBFast.h>
#include <SxCLI.h>
#include <SxCliStandard.h>
#include <SxBinIO.h>
#include <SxAtomicStructure.h>
#include <SxSpeciesData.h>

int main (int argc, char **argv)
{
   // --- parse command line options
   SxCLI cli (argc, argv);
   
   // --- input and output files
   // standard structure reading options
   bool sxbFile    = sxbOption   (cli);
   SxString inFile = inputOption (cli, sxbFile);

   SxString outFile 
      = cli.option("-o","file","PDB output file to be written")
        .toString("input.pdb");
   
   // --- repetition
   cli.option ("-r|--repeat","AxBxC",
               "with numbers A, B, C, and x being an 'x'. "
               "Repeat the cell A times in first (x) direction, B times in "
               "second (y) direction and C times in 3rd (z) direction");
   SxVector3<int> xyzRepeat (1,1,1);
   if (cli.last ().exists ())
      xyzRepeat = SxVector3<int> (cli.last ().toIntList3 ("x,"));
   cli.last ().defaultValue = "default: 1x1x1";
   cli.finalize ();
   
   SxAtomicStructure structure = structureFromInput (inFile, sxbFile);
   
   initSPHInXMath ();

   if (xyzRepeat.product () > 1)
      structure = structure.repeat (xyzRepeat);
   const SxArray<SxString> &chemName = structure.getElements ();
   SxPDBFast(outFile).write (structure, chemName);
   return 0;
}

