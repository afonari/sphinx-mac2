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
#include <SxNeighbors.h>
#include <SxGrid.h>
#include <SxFileIO.h>

SxVector3<double> toCoords (const SxString &s)
{
   SxList<SxString> coords
      = s.substitute ("{","")
         .substitute ("}","")
         .substitute (" ","")
         .tokenize (',');
   if (coords.getSize () != 3)  {
      throw SxException (("Cannot extract coordinates from " + s).ascii ());
   }
   // => might throw an exception, too
   Coord pos(coords(0).toDouble (),
             coords(1).toDouble (),
             coords(2).toDouble ());
   return pos;
}

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage =
      "This add-on modifies the structure according to the structure patch "
      "file as produced by sxstructdiff.";
   cli.authors = "C. Freysoldt";

   SxString inFile
      = cli.option ("-i|--input", "input file",
                    "take original input file")
        .toString ("input.sx");

   SxString outFile
      = cli.option ("-o","filename", "output file name (screen otherwise)")
        .toString ("");

   SxString patchFile
      = cli.option ("-p","filename", "patch file name")
        .toString ();

   double distMax = cli.option ("-d|--dmax","distance",
                                "tolerance for finding atoms").toDouble (0.);

   bool labels = cli.option ("-l|--labels", "transfer labels").toBool ();

   bool keepMovable
      = cli.option("-m|--keep-movable","retain mobility of atoms as "
                   "it is defined in the input structure")
        .toBool ();
   cli.finalize ();

   // --- read input
   SxAtomicStructure structure = structureFromInput (inFile, false, keepMovable);
   SxArray<SxString> chemNames = structure.getElements ();
   if (!structure.hasLabels ()) labels=false;

   SxArray<bool> keep(structure.getNAtoms ());
   keep.set (true);
   SxArray<SxString> patch = SxString(SxFileIO::readBinary (patchFile,-1))
                             .tokenize ('\n');
   SxGrid grid(structure, 10);

   SxArray<SxList<Coord> > addedAtoms(structure.getNSpecies ());
   for (int ip = 0; ip < patch.getSize (); ++ip)  {
      if (patch(ip).getSize () == 0) continue;
      if (patch(ip)(0) != '>') continue;
      if (patch(ip).contains ("> new"))  {
         // --- additional atom
         SxList<SxString> split = patch(ip).tokenize ('@');
         if (split.getSize () != 2)  {
            if (split.getSize () > 0)  {
               cout << "Cannot parse this line (ignored):" << endl;
               cout << patch(ip) << endl;
            }
            continue;
         }
         SxString element = split(0).right ("new").substitute (" ","");
         int is = (int)chemNames.findPos (element);
         if (is == -1)  {
            is = (int)chemNames.getSize ();
            chemNames.resize (is + 1, true);
            chemNames(is) = element;
            addedAtoms.resize (is + 1, true);
         }

         Coord pos;
         try {
            pos = toCoords (split(1));
         } catch (SxException e) {
            e.print ();
            continue;
         }

         addedAtoms(is) << pos;
         cout << "Adding " << element << " atom @ " << pos << endl;
      }

      // --- atoms were deleted, shifted, or changed species
      SxList<SxString> split = patch(ip).tokenize (':');
      if (split.getSize () != 2) {
         if (split.getSize () > 2)  {
            cout << "Cannot parse this line (ignored):" << endl;
            cout << patch(ip) << endl;
         }
         continue;
      }

      // --- find the relevant atom
      Coord pos;
      try {
         pos = toCoords (split(0).right ("@"));
      } catch (SxException e)  {
         e.print ();
         continue;
      }

      int iTlAtom = structure.find (pos, grid);
      if (iTlAtom == -1 && distMax > 0.)  {
         // --- try atoms nearby
         SxNeighbors nn;
         nn.compute (grid, structure, pos, distMax,
                     SxNeighbors::StoreIdx | SxNeighbors::StoreAbs
                     | SxNeighbors::IncludeZeroDistance);
         if (nn.getSize () == 1) iTlAtom = nn.idx(0);
         else if (nn.getSize () > 0) cout << nn.absPositions << endl;
      }
      if (iTlAtom == -1)  {
         cout << "Could not find atom @ " << pos << endl;
         continue;
      }

      // --- now apply the patch
      //
      // --- delete atom
      if (split(1).contains ("delete"))  {
         keep(iTlAtom) = false;
         cout << "Deleting atom " << (iTlAtom+1) << " @ "
              << structure.getAtom (iTlAtom) << endl;
         continue;
      }

      // --- shift atom
      if (split(1).contains ("shift"))  {
         Coord by;
         try {
            by = toCoords (split(1).right ("shift"));
         } catch (SxException e)  {
            e.print ();
            continue;
         }
         cout << "Shifting atom " << (iTlAtom+1) << " @ "
              << structure.getAtom (iTlAtom)
              << " by " << by
              << endl;
         structure.ref (iTlAtom) += by;
         continue;
      }

      // --- change species (delete and add)
      if (split(1).contains ("new species"))  {
         keep(iTlAtom) = false;
         if (keepMovable)  {
            SxVector3<int> movable
               = structure.atomInfo->meta.get<SxArray<SxVector3<int> > >
                (SxAtomicStructure::StickyFilter)(iTlAtom);
            if (movable != SxVector3<int>(1, 1, 1))
               cout << "Warning: atom " << (iTlAtom+1) << "@ "
                    << structure.getAtom (iTlAtom) << " becomes movable upon changing species"
                    << endl;
         }

         SxString element = split(1).right ("species").substitute (" ","");
         int is = (int)chemNames.findPos (element);
         if (is == -1)  {
            is = (int)chemNames.getSize ();
            chemNames.resize (is + 1, true);
            chemNames(is) = element;
            addedAtoms.resize (is + 1, true);
         }
         addedAtoms(is) << structure.getAtom (iTlAtom);
         cout << "Changing species for atom " << (iTlAtom+1)
              << " @ " << structure.getAtom (iTlAtom) << " from "
              << chemNames(structure.getISpecies(iTlAtom)) << " to "
              << element << endl;
      }
   }

   // --- copy atoms that have not been deleted
   SxAtomicStructure result(structure.cell);
   SxStack<SxString> newLabels;
   const SxArray<SxString> *oldLabels = labels ? &structure.getLabels () : NULL;
   for (int is = 0, iTl = 0; is < chemNames.getSize (); ++is)  {
      result.newSpecies ();
      if (is  < structure.getNSpecies ())  {
         for (int ia = 0; ia < structure.getNAtoms (is); ++ia, ++iTl)  {
            if (keep(iTl))  {
               result.addAtom (structure.getAtom (is, ia));
               if (oldLabels) newLabels << (*oldLabels)(iTl);
            }
         }
      }
      for (int ia = 0; ia < addedAtoms(is).getSize (); ++ia)  {
         result.addAtom (addedAtoms(is)(ia));
         newLabels << "";
      }
   }
   result.endCreation ();
   result.atomInfo->meta.attach (SxAtomicStructure::Elements, chemNames);
   if (labels)  {
      result.atomInfo->meta.update (SxAtomicStructure::Labels,
                                    SxArray<SxString> (newLabels));
   }
   if (keepMovable)  {
      const SxArray<SxVector3<int> > &orig
         = structure.atomInfo->meta.get (SxAtomicStructure::StickyFilter);
      SxArray<SxVector3<int> > stickyArray(result.getNAtoms ());
      stickyArray.set (SxVector3<int> (1, 1, 1));
      for (int is = 0; is < structure.getNSpecies (); ++is)  {
         for (int ia = 0, ja = 0; ia < structure.getNAtoms (is); ++ia)  {
            int iTl = structure.getIAtom (is, ia);
            if (keep(iTl))
               stickyArray (result.getIAtom (is, ja++)) = orig(iTl);
         }
      }
      result.atomInfo->meta.attach (SxAtomicStructure::StickyFilter,
                                    stickyArray);
   }

   // --- output
   writeStructure (result, outFile);

}
