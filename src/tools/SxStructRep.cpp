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
#include <SxCliStandard.h>
#include <SxAtomicStructure.h>
#include <SxGrid.h>
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
   SxMatrix3<int> repMatrix;

   int simpleRep = cli.newGroup ("simple repetition");

   SxVector3<int> diag (cli.option ("-r|--repeat","3 numbers",
                                     "repetition along current cell vectors")
                            .toIntList3 ("x,"));

   int trueRep = cli.newGroup ("3dimensional repetition");
   cli.excludeGroup (simpleRep);
   SxVector3<int> col1 (cli.option ("--a1","vector",
                                     "new a1 in relative coordinates")
                         .toIntList3 ()),
                  col2 (cli.option ("--a2","vector",
                                     "new a2 in relative coordinates")
                         .toIntList3 ()),
                  col3 (cli.option ("--a3","vector",
                                     "new a3 in relative coordinates")
                         .toIntList3 ());
   if (!cli.error)  {
      if (cli.groupAvailable (simpleRep))  {
         if (diag.product () <= 0)  {
            cout << "Illegal repetition factors " << diag << endl;
            cli.setError ();
         }
      } else if (cli.groupAvailable (trueRep))  {
         repMatrix = SxMatrix3<int> (col1, col2, col3).transpose ();
         cout << repMatrix << endl;
         if (repMatrix.determinant () == 0)  {
            cout << "Illegal repetition matrix with zero determinant." << endl;
            cli.setError ();
         }
      } else {
         cout << "no repetition given!" << endl;
         cli.setError ();
      }
   }
   cli.setGroup (cli.generalGroup);
   
   bool wrap = cli.option ("--wrap", "keep atoms within cell").toBool ();

   // standard structure reading options
   bool sxbFile    = sxbOption   (cli);
   SxString inFile = inputOption (cli, sxbFile);

   // standard output options (-o, --outsxb)
   SxStructOut out (cli, sxbFile);

   bool labels = cli.option ("-l|--labels", "transfer labels").toBool ();
   bool movable = cli.option ("-m|--movable", "transfer movable").toBool ();
   cli.last ().defaultValue = "default: same format as input file";

   cli.finalize ();

   // --- read input
   SxAtomicStructure structure = structureFromInput (inFile, sxbFile, movable);

   SxAtomicStructure newStructure 
      = cli.groupAvailable (simpleRep) ? structure.repeat (diag)
                                       : structure.repeat (repMatrix);

   if (wrap) newStructure %= newStructure.cell;

   cout << "nSym = " << newStructure.cell.symGroupPtr->getSize () << endl;

   if (labels || movable)  {
      SxGrid grid(structure,10);
      SxConstPtr<SxAtomInfo> info = structure.match (grid, newStructure);
      newStructure.replaceInfo (info);
   }
   if (labels && structure.hasLabels ())  {
      const SxArray<SxString> &oldLabels = structure.getLabels ();
      SxArray<SxString> newLabels(newStructure.getNAtoms ());
      SX_CHECK (newStructure.atomInfo->parentMap.getSize () > 0);
      SX_LOOP(iNew)
         newLabels(iNew) = oldLabels(newStructure.atomInfo->parentMap(iNew));
      newStructure.atomInfo->meta.update (SxAtomicStructure::Labels, newLabels);
   }
   if (movable)  {
      SxConstPtr<SxAtomInfo> &info = newStructure.atomInfo;
      SX_CHECK (info->parentMap.getSize () > 0);
      SxArray<SxVector3<int> > &oldMovable
         = structure.atomInfo->meta.get (SxAtomicStructure::StickyFilter);
      SxArray<SxVector3<int> > newMovable(newStructure.getNAtoms ());
      SX_LOOP(iNew)
         newMovable(iNew) = oldMovable(info->parentMap(iNew));
      newStructure.atomInfo->meta.update (SxAtomicStructure::StickyFilter,
                                          newMovable);
   }

   // --- output
   out.write (newStructure);
}

