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
#include <SxPW.h>
#include <SxBinIO.h>
#include <SxFermi.h>
#include <SxPerturbK.h>


int main (int argc, char **argv)
{

   initSPHInXMath ();

   SxCLI cli (argc,argv);
   cli.authors = "C. Freysoldt";
   cli.preUsageMessage = "This add-on calculates d(eps)/dk.\n";
   /*
      We do this to test against numerical derivatives:
      sxepsslope ... --path | sed -e'1,/Start/d' | xmgrace -nxy -
      Data->transformations->Differences->Central differences for first set 

      Should match if k-direction is along x, y, or z.
      */

   SxString wavesFile = cli.option ("-w|--waves","file","waves file")
                        .toString("waves.sxb");

   int n = cli.option ("-n","index","state").toInt (false,1) - 1;
   bool shortTangent = cli.option ("--shortTangent", "request short tangents")
                       .toBool ();
   bool longTangent  =  cli.option ("--longTangent", "request long tangents")
                       .toBool ();
   bool printTangent = shortTangent || longTangent;

   SxPerturbK kp (cli);

   bool path
      = cli.option ("--path", "first col in output is path length").toBool ();

   int precision
      = cli.option ("--prec", "digits", "set number of digits for output")
        .toInt (6);

   cli.version ("1.2");
   cli.finalize ();

   cout.precision (precision);


   SxPW waves;
   SxFermi fermi;

   SxAtomicStructure structure;

   try {
      SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
      waves.read (io);
      fermi.read (io);
      structure.read(io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   SxGkBasis &gkBasis = waves.getGkBasis();
   gkBasis.changeTau(structure);
   fermi.kpPtr = &gkBasis;
   
   int ik, idir, jk;
   // int nStates = waves.getNStates ();
   int nk = waves.getNk ();
   int iSpin = 0;

   SxArray<SxVector3<double> > npn(nk);

   // kp is SxPerturbK set up from cli (see above) and gkBasis
   kp.setup (gkBasis);
   
   for (ik = 0; ik < nk; ik++)  {
      (cout << "ik = " << ik << endl).flush ();
      npn(ik) = kp.getMatrixElements (waves(n,iSpin,ik),
                                      waves(n,iSpin,ik)).toVector3 ();
   } // ik

   kp.printTimer ();
   
   // --- calculate lines E_nj(k) = <nj|p|nj> * (k - k_j) + eps(n,jk)
   cout << "Start" << endl;
   double pathLength = 0.;
   for (ik = 0; ik < nk; ik++)  {
      if (path)  {
         cout << pathLength;
         if (ik + 1 < nk)
            pathLength += (gkBasis.kVec(ik+1) - gkBasis.kVec(ik)).norm ();
      } else
         cout << ik;
      cout << "\t" << fermi.eps(n,iSpin,ik) * HA2EV;
      if (printTangent)  {
         for (jk = 0; jk < nk; jk++)  {
            if (shortTangent && abs(ik-jk) > 1)  {
               cout << "\t0.0";
            } else {
               SxVector3<double> dk = gkBasis.kVec(ik) - gkBasis.kVec(jk);
               double tangentVal = (npn(jk) ^ dk) + fermi.eps(n,iSpin,jk);
               cout << "\t" << (tangentVal * HA2EV);
            }
         }
      } else {
         for (idir = 0; idir < 3; idir++)
            cout << "\t" << npn(ik)(idir) * HA2EV;
      }
      cout << endl;
   }
   
}

