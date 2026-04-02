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

#include <SxParser.h>
#include <SxRadBasis.h>
#include <SxAtomicOrbitals.h>
#include <SxCLI.h>
#include <SxRadialAtom.h>
#include <SxEigensystem.h>
#include <SxTextIO.h>

int main (int argc, char** argv)
{
   initSPHInXMath ();

   cout.precision(10);

   // Command line parsing
   SxCLI cli (argc,argv);

   // Define Author
   cli.authors = "B. Lange";

   // What does the program do ?
   cli.preUsageMessage = "AO Basis joining tool.";

   SxString basisFile = cli.option ("-b|--basis", "file", "AOBasis input file")
                                    .toString ("basis.sx");
   SxString outputFile = cli.option ("-o|--out", "file", "AOBasis sxb file")
                                    .toString ("basis.sxb");
   enum OrthoMode { DiagLaplace = 0, DiagR2 = 1 };
   int orthoMode = cli.option ("--laplace|--r2",
                               "diagonalize l-channel Laplacian/r^2 matrix")
                   .toChoice ();
   cli.finalize();

   SxParser parser;
   SxConstPtr<SxSymbolTable> table = parser.read (basisFile);
   SxAtomicOrbitals orbitals;
   
   try   {
         SxSymbolTable *basisGroup = table->getGroup("AOBasis");
         orbitals.setup(basisGroup);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   SxConstPtr<SxRadBasis> radBasisPtr = orbitals.getBasis ();
   
   if (orthoMode > -1)  {
      for (int is = 0; is < orbitals.getNSpecies (); ++is)  {
         SxVector<double> r2 = radBasisPtr->radFunc(is).sqr ();
         for (int l = 0; l < orbitals.getLMax (is); ++l)  {
            int n = orbitals.getFuncPerL(is,l);
            if (n <= 1) continue;
            SxVector<double> S(n,n), R(n,n);
            SxVector<double> phi(r2.getSize (), n);
            for (int i = 0; i < n; ++i)  {
               phi.colRef (i) <<= orbitals.getFuncL(is,l,i);
               for (int j = 0; j < n; ++j)  {
                  S(j,i) = tr(orbitals.getFuncL(is,l,i) * orbitals.getFuncL(is,l,j));
                  if (orthoMode == DiagR2)
                     R(j,i) = tr(orbitals.getFuncL(is,l,i) 
                                 * orbitals.getFuncL(is,l,j) * r2);
                  else {
                     SX_CHECK (orthoMode == DiagLaplace);
                     SxVector<double> L;
                     L = SxRadialAtom::laplace (orbitals.getFuncL(is,l,j));
                     R(j,i) = -tr(orbitals.getFuncL(is,l,i) * L);
                  }
               }
            }
            cout << "is=" << is << "; l = " << l << "; n=" << n << endl;
            cout << "S=" << S << endl;
            SxVector<double> L = S.choleskyDecomposition ().inverse ();
            cout << "R=" << R << endl;
            R = L.transpose () ^ R ^ L;
            cout << "R(ortho)=" << R << endl;
            SxSymEigensystem<double> eig (R);
            eig.vecs = L ^ eig.vecs;
            cout << "eigenvals R: " << eig.vals << endl;
            phi = phi ^ eig.vecs;
            for (int i = 0; i < n; ++i)  {
               orbitals.getFuncL(is,l,i) <<= phi.colRef (i);
            }
         }
         for (int ipt = 0; ipt < orbitals.getNOrbTypes(is); ++ipt)  {
            SxTextIO ("orbital-" + SxString(is) + "-" + SxString(ipt) + ".dat")
               .writeXYPlot (radBasisPtr->radFunc(is), orbitals(is, ipt), 10, 10);
         }
      }
   }

   orbitals.write(outputFile);
}
