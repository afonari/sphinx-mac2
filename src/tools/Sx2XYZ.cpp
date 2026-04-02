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
#include <SxXYZ.h>
#include <SxCLI.h>
#include <SxCliStandard.h>
#include <SxAtomicStructure.h>
#include <SxSpeciesData.h>
#include <SxSimpleParser.h>
#include <SxIO.h>

int main (int argc, char **argv)
{
   // --- parse command line options
   SxCLI cli (argc, argv);

   // --- input and output files

   // standard structure reading options
   bool sxbFile    = sxbOption   (cli);
   SxString inFile = inputOption (cli, sxbFile);

   SxString outFile = cli.option("-o","file","XYZ output file to be written")
                      .toString("input.xyz");
   bool withAMat = cli.option ("-a", "prepend lattice vectors (a1, a2, a3)")
                   .toBool();
   bool trajectory = cli.option ("--trajectory", "file is a trajectory")
                     .toBool ();
   int interpolate = cli.option ("--interpolate", "number of points",
                        "trajectory interpolation with additional points")
                     .toInt (0,0);

   int simpleRep = cli.newGroup ("simple repetition");
   SxVector3<int> repeat(0,0,0);
   cli.option ("-r|--repeat","3 numbers",
               "repetition along current cell vectors");
   if (cli.last ().exists ())
      repeat = SxVector3<int> (cli.last ().toIntList3 ("x,"));

   int trueRep = cli.newGroup ("3dimensional repetition");
   SxMatrix3<int> repMatrix;
   repMatrix.set (0);
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
         if (repeat.product () <= 0)  {
            cout << "Illegal repetition factors " << repeat << endl;
            cli.setError ();
         }
      } else if (cli.groupAvailable (trueRep))  {
         repMatrix = SxMatrix3<int> (col1, col2, col3).transpose ();
         cout << repMatrix << endl;
         if (repMatrix.determinant () == 0)  {
            cout << "Illegal repetition matrix with zero determinant." << endl;
            cli.setError ();
         }
      }
   }
   //cli.setGroup (cli.generalGroup);
   bool doWrap = cli.option ("--wrap", "keep atoms within cell").toBool ();

   cli.finalize ();

   initSPHInXMath ();

   // --- non-trajectories
   if (!trajectory)  {
      // --- read input
      SxAtomicStructure structure = structureFromInput (inFile, sxbFile);

      // --- repeat
      if (repeat.product () > 0)
         structure = structure.repeat (repeat);
      else if (repMatrix.determinant () != 0)
         structure = structure.repeat (repMatrix);
      if (doWrap) structure %= structure.cell;

      // --- write
      const SxArray<SxString> &chemName = structure.getElements ();
      SxXYZ (outFile).write (structure, chemName, withAMat);
      return 0;
   }

   // --- trajectories
   SxAtomicStructure structure;
   SxArray<SxString> chemName;
   FILE *fp =sxfopen (outFile, "w");
   // --- read input
   SxParser parser;
   SxParser::Table table = parser.read (inFile, "std/forces.std");

   // --- read elements
   chemName = SxSpeciesData::getElements (&*table);
   SX_LOOP(is) chemName(is) = chemName(is).toUpper ();

   // --- now loop over single frames/steps
   int it = 0;
   SxAtomicStructure wrap, oldStructure;
   SYMBOLPARSE(&*table)  {
      FOREACH_SYMBOLGROUP("structure")  {
         // --- read new structure
         structure = SxAtomicStructure(SYMBOLGROUP_TABLE);

         // --- repeat
         if (repeat.product () > 0)
            structure = structure.repeat (repeat);
         else if (repMatrix.determinant () != 0)
            structure = structure.repeat (repMatrix);

         // --- wrap
         if (doWrap)  {
            if (it == 0)
               wrap = (structure % structure.cell) - structure;
            structure += wrap;
         }
         SX_CHECK(structure.getNSpecies () == chemName.getSize (),
                  structure.getNSpecies (), chemName.getSize ());

         // --- interpolate and write interpolated steps
         if (interpolate > 0 && oldStructure.getNAtoms () > 0)  {
            for (int ip = 1; ip <= interpolate; ++ip)  {
               double x = double(ip)/(interpolate+1);
               fprintf (fp, "%d\nStep %g\n", structure.getNAtoms (),
                        it + x);
               SX_LOOP2(is,ia)  {
                  const Coord &newPos = structure.getAtom (is,ia);
                  Coord pos = oldStructure.getAtom (is,ia);
                  pos += x * (newPos - pos);
                  fprintf (fp, "%s %.8f %.8f %.8f\n", chemName(is).ascii (),
                           pos(0) / A2B, pos(1) / A2B, pos(2) / A2B);
               }
            }
         }
         oldStructure = structure;

         // --- write
         fprintf (fp, "%d\nStep %d\n", structure.getNAtoms (), ++it);
         SX_LOOP2(is,ia)  {
            const Coord &pos = structure.getAtom (is,ia);
            fprintf (fp, "%s %.8f %.8f %.8f\n", chemName(is).ascii (),
                     pos(0) / A2B, pos(1) / A2B, pos(2) / A2B);
         }
      }
   }
   fclose (fp);
   return 0;
}

