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
#include <SxTypes.h>
#include <SxCell.h>
#include <SxFFT2d.h>
#include <SxFileParser.h>

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.preUsageMessage =
      "This tool finds isosurface positions above the a0 x a1 mesh.";

   double target = cli.option ("--value|--target|-t", "value",
                               "target value of the isosurface")
                 .toDouble ();

   SxString outFile
      = cli.option ("-o","filename", "output file name (screen otherwise)")
        .toString ("");

   SxFFT::plannerCLI (cli);

   SxString meshFile =
      cli.option ("-i","mesh file","netcdf mesh file (one mesh only)")
      .toString ();


   int iGroup = cli.newGroup ("points from file");
   SxString pointsFile = cli.option ("-p|--points", "file", "files with x y")
                         .toString ("");
   bool relative = cli.option ("--relative", "assume relative coordinates")
                   .toBool ();
   int zGroup = cli.newGroup ("zRange");
   double zMin = cli.option ("--zmin","val", "min. z value").toDouble ();
   double zMax = cli.option ("--zmax","val", "max. z value").toDouble ();
   cli.setGroup (cli.generalGroup);


   cli.finalize ();

   // --- read input

   SxMesh3D mesh;
   SxCell cell;
   SxMeshR data;
   try {
      SxBinIO io (meshFile, SxBinIO::BINARY_READ_ONLY);
      data = io.readMesh (&cell, &mesh)(0);
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   cell.setup ();

   // --- points to interpolate
   SxList<Coord> rVec;
   if (cli.groupAvailable (iGroup) && pointsFile.getSize () > 0)  {
      SxFileParser fp(pointsFile);
      Coord xyz;
      xyz(2) = 0.;
      while (!feof(fp.fp))  {
         fp >> xyz(0) >> xyz(1);
         fp.skipWhite ();
         if (!relative) {
            xyz(2) = 0.;
            xyz = cell.carToRel (xyz);
         }
         rVec << xyz;
      }
   }

   // output
   FILE *file = (outFile.getSize () > 0) ? sxfopen(outFile,"w") : stdout;

   SxMesh3D mesh2D(mesh(0), mesh(1), 1);
   SxVector<double> isoPos(mesh2D.getSize ());
   SxVector3<int> i; // fft index
   for (i(0) = 0; i(0) < mesh(0); i(0)++)  {
      for (i(1) = 0; i(1) < mesh(1); i(1)++)  {
         i(2) = mesh(2) - 1;
         double prev = data(mesh.getMeshIdx(i,SxMesh3D::Positive));
         for (i(2) = 0; i(2) < mesh(2); i(2)++)  {
            ssize_t idx = mesh.getMeshIdx (i, SxMesh3D::Positive);
            double val = data(idx);
            if (  (val <= target && prev > target)
                ||(val >= target && prev < target))
            {
               double c = 0.5 * (prev + data(mesh.getMeshIdx(i(0),i(1),i(2)+1,SxMesh3D::Unknown)) - 2. * val);
               double x = (val - target)/(val-prev);
               double gamma = c / (prev - val);
               if (fabs(c) > 1e-8)  {
                  double g1 = 0.5 * (gamma - 1.);
                  double D = sqrt(g1 * g1 + gamma * x);
                  //if (gamma > 0. || gamma < -1.)
                     x = (g1+D)/gamma;
                  //else
                  //   x = (g1 - D)/gamma;
               }
               x = i(2) - x;
               SxVector3<double> pos(i(0), i(1), x);
               pos /= mesh;
               pos = cell.relToCar (pos);
               if (!cli.groupAvailable (zGroup)
                   || (pos(2) >= zMin && pos(2) <= zMax))
               {
                  isoPos(mesh2D.getMeshIdx(i(0),i(1),0,SxMesh3D::Positive)) = pos(2);
                  if (!cli.groupAvailable (iGroup))
                     sxfprintf (file, "%.6f %.6f %.6f\n", pos(0), pos(1), pos(2));
               }
            }
            prev = val;
         }
      }
   }

   if (cli.groupAvailable (iGroup) && rVec.getSize () > 0)  {
      SX_VALIDATE_VECTOR(isoPos);
      SxFFT2d fft2d(SxFFT::Reverse, mesh(0), mesh(1));
      SxVector<SxComplex16> isoPosR = isoPos, isoPosG (isoPos.getSize ());
      SX_VALIDATE_VECTOR(isoPosR);
      fft2d.fftReverse ((int)isoPos.getSize (), isoPosR.elements, isoPosG.elements);
      SX_VALIDATE_VECTOR(isoPosG);
      for (SxList<Coord>::Iterator it = rVec.begin ();
           it != rVec.end (); ++it)
      {
         double posZ = 0.;
         SX_LOOP(ig)  {
            Coord G = mesh2D.getMeshVec(ig, SxMesh3D::Origin);
            SxComplex16 phase = SxComplex16::phase (TWO_PI * (G ^ *it));
            posZ += (isoPosG(ig) * phase).re;
         }
         //posZ /= double(isoPos.getSize ());
         Coord xy = cell.relToCar (*it);
         sxfprintf (file, "%.6f %.6f %.6f\n", xy(0), xy(1), posZ);
      }
   }

   if (file != stdout) fclose (file);

}

