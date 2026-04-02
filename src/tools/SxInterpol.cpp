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
#include <SxUtil.h>
#include <SxVector.h>
#include <SxTextIO.h>
#include <SxParser.h>
#include <SxAtomicStructure.h>
#include <SxGkBasis.h>
#include <SxRadBasis.h>
#include <SxPAWPot.h>
#include <SxCubicSpline.h>

double noise ();

typedef SxArray<SxVector<double> > SxSplineType;

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);
   cli.authors = "B.Lange";
   cli.preUsageMessage = "Least square fit of natural cubic spline";
   int dim = cli.option ("--nPoints", "number", "number of points")
             .toInt (10);

   cli.newGroup ("Testing");
   int fitGroup = cli.newGroup ("External fit");
   SxString inFile = cli.option ("--input", "file", "data file").toString ();
   cli.setGroup (cli.generalGroup);
   cli.version ("0.1");
   cli.finalize ();

   bool externalFit = cli.groupAvailable (fitGroup);


   SxVector<double> xVector, xRes;
   SxVector<double> yVector, yRes;
   int dimRes;
   double deltaRes;

   if (externalFit)  {
      ifstream data(inFile.ascii ());
      if (!data)  {
         cout << "Cannot open " << inFile << endl;
         SX_QUIT;
      }
      SxStack<double> xIn, yIn;
      while (data)  {
         double x,y;
         data >> x >> y;
         xIn << x; 
         yIn << y;
      }
      xVector = SxVector<double> (xIn);
      yVector = SxVector<double> (yIn);
      int nData = int(xVector.getSize ());

      // --- set up x values for support and sampling
      double xMin = xVector.minval ();
      double xMax = xVector.maxval ();
      dimRes = 10 * nData;
      xRes.resize (dimRes);
      double dx = (xMax - xMin) / (dimRes - 1);
      for (int i = 0; i < dimRes; ++i)  {
         xRes(i) = i * dx + xMin;
      }

   } else {
      // setup yVectorian
      double delta = 10 / double(dim);
      xVector.resize (dim);
      yVector.resize (dim);
      for (int i = 0; i < dim; i++)   {
         xVector(i) = i * delta;
         yVector(i) = (i*delta - 5.0) * (i*delta - 5.0);
      }
      SxTextIO ("FuncIN.dat").writeXYPlot(xVector,yVector);

      dimRes = 10 * dim;
      deltaRes = (xVector(dim-1) - xVector(0)) / dimRes;
      xRes.resize (dimRes);
      for (int i = 0; i < dimRes; i++)
         xRes(i) = i * deltaRes + xVector(0);
   }

   SxCubicSpline spline (xVector, yVector, SxCubicSpline::Natural);
   yRes = spline.getY(xRes);

   for (int i = 0; i < yRes.getSize(); i++)  {
      yRes(i) += noise ();
   }

   SxTextIO ("FuncSplineNoise.dat").writeXYPlot(xRes,yRes);

   SxCubicSpline splineFit (xRes, yRes, xVector, SxCubicSpline::Natural,
                            SxCubicSpline::Normal);
   SxVector<double> yFit = splineFit.getYFit ();

   SxTextIO ("FuncFit.dat").writeXYPlot(xVector,yFit);

   return 0;
}

double noise ()
{
   double x = SxRandom::get() * 2 * PI;
   double result = 0.05 * sin(x);
   return result;
}

