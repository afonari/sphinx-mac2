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
#include <SxAims.h>
#include <SxCLI.h>
#include <SxCliStandard.h>
#include <SxAtomicStructure.h>
#include <SxSpeciesData.h>

int main (int argc, char **argv)
{
   // --- parse command line options
   SxCLI cli (argc, argv);
   
   // --- input and output files
   bool sxbFile = sxbOption (cli); // standard --sxb option
   SxString inFile = inputOption (cli, sxbFile); // standard -i|--input option

   SxString outFile = cli.option("-o","file","FHIaims structure file to "
                                 "be written")
                      .toString("geometry.in");

   cli.finalize ();
   
   // --- read from input file
   initSPHInXMath ();

   SxAtomicStructure structure = structureFromInput (inFile, sxbFile);
   const SxArray<SxString> &chemName = structure.getElements ();
   
   SxAims aims (outFile);
   aims.write(structure, chemName);
   
   return 0;
}

