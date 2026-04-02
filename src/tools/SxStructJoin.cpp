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
      "This add-on allows to join the atoms of two structures. The final "
      "structure will be the first structure plus all the atoms of the second "
      "structure, with their cartesian coordinates. The origin of the second "
      "structure can be optionally shifted with respect to the first one.";
   cli.authors = "C. Freysoldt";

   // ---  definition
   SxVector3<double> trans(0.,0.,0.);
   
   SxCLI::CliArg *opt = &cli.option ("--by|--vector","vector", 
                                     "translation vector");
   if (opt->exists ())  {
      trans = SxVector3<double>(opt->toList3 ());
   }
   opt->required (false);
   
   bool wrap = cli.option ("--wrap", "keep atoms within cell").toBool ();

   SxString substrateFile 
      = cli.option ("-i|--substrate", "input file", "S/PHI/nX structure file "
                    "for first structure (determines cell)")
        .toString ("input.sx");
   
   SxString adsorbateFile 
      = cli.option ("-j|--join|--adsorbate", "input file", 
                    "S/PHI/nX structure file for second structure")
        .toString ();
   
   // standard output options (-o, --outsxb)
   SxStructOut out (cli);

   cli.finalize ();

   // --- read input
   SxAtomicStructure substrate = structureFromInput (substrateFile),
                     adsorbate = structureFromInput (adsorbateFile);
   // shift adsorbate
   adsorbate += trans;

   // --- join structures
   const SxArray<SxString> &chemName = substrate.getElements (),
                           &adsorbChemName = adsorbate.getElements ();
   int is, ia, joinedSpecies;
   int nSpecies = substrate.getNSpecies ();
   SxAtomicStructure structure(substrate, SxAtomicStructure::Copy);
   SxList<SxString> chemNames;
   chemNames << chemName;
   structure.startCreation ();
   for (is = 0; is < adsorbate.getNSpecies (); ++is)  {
      cout << "Joining " << adsorbChemName(is) << "..." << endl;
      joinedSpecies = int(chemNames.findPos (adsorbChemName(is)));
      if (joinedSpecies < 0)  {
         joinedSpecies = nSpecies++;
         structure.newSpecies ();
         chemNames << adsorbChemName(is);
      }
      for (ia = 0; ia < adsorbate.getNAtoms(is); ++ia)
         structure.addAtom (joinedSpecies,adsorbate(is,ia));
   }
   structure.endCreation ();
   structure.atomInfo->meta.attach (SxAtomicStructure::Elements,
                                    SxArray<SxString> (chemNames));

   if (wrap) structure %= structure.cell;

   // --- output
   out.write (structure);
}
