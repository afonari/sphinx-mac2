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
#include <SxBinIO.h>
#include <SxAtomicStructure.h>

template <class T>
SxList<T> &operator<< (SxList<T> &list, const SxArray<T> &array)
{
   for (int i = 0; i < array.getSize (); ++i)
      list.append(array(i));
   return list;
}

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage =
      "This add-on interpolates the atomic positions between two structures.";
   cli.authors = "C. Freysoldt";

   // ---  definition
   bool wrap = cli.option ("--wrap", "keep atoms within cell").toBool ();

   SxString inFile1
      = cli.option ("-i|--initial", "input file", "S/PHI/nX structure file "
                    "for first structure")
        .toString ();

   SxString inFile2
      = cli.option ("-f|--final", "input file",
                    "S/PHI/nX structure file for second structure")
        .toString ();

   double x = cli.option ("-x|--mixing", "value",
                          "linear structure mixing parameter: 0=initial 1=final")
              .toDouble (0.5);

   bool keepMovable
      = cli.option("-m|--keep-movable","retain mobility of atoms as "
                   "it is defined in the input structure")
        .toBool ();

   // standard output options (-o, --outsxb)
   SxStructOut out (cli);

   cli.finalize ();

   // --- read input
   SxAtomicStructure structure = structureFromInput (inFile1, false, keepMovable),
                     end = structureFromInput (inFile2);

   if (structure.getNAtoms () != end.getNAtoms ()
       || (structure.atomInfo->nAtoms - end.atomInfo->nAtoms).sqr ().sum () > 0)
   {
      cout << "Structure mismatch - must be the same number of atoms in each species"
           << endl;
      SX_QUIT;
   }
   end.atomInfo = structure.atomInfo;

   SxAtomicStructure delta = end - structure;
   cout << "transDir = " << delta.coordRef () << ";" << endl;

   structure += x * delta;

   if (wrap) structure %= structure.cell;

   // --- output
   out.write (structure);
}
