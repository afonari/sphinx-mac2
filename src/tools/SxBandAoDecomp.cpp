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
#include <SxPW.h>
#include <SxBinIO.h>
#include <SxIO.h>
#include <SxAOBasis.h>
#include <SxFermi.h>


int main (int argc, char **argv)
{

   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   SxCLI cli (argc,argv);
   cli.authors = "C. Freysoldt";
   cli.preUsageMessage = "This is the Mulliken-like band analysis tool.";

   SxString wavesFile = cli.option ("-w|--waves","file","waves file")
                        .toString("waves.sxb");

   SxString inFile = cli.option ("-i|--input", "file", "S/PHI/nX input file")
                     .toString ("input.sx");

   SxList<int> states, kPoints;
   // --- selection by state number
   states = cli.option ("-n|--states","list",
                        "states (starting from 1) to be selected. "
                        "They are given as a ','-separated list of numbers or "
                        "ranges a-b, e.g. '-n 1,5-8,10' would select the states"
                        " 1, 5, 6, 7, 8, and 10.")
            .toIdxList ();
   cli.last ().defaultValue = "default: all states";
   kPoints = cli.option ("-k|--kpoints","list",
                         "k-points (starting from 1) to be "
                         "selected. The list is given like for the -n option")
             .toIdxList ();
   cli.last ().defaultValue = "default: all k-points";

   bool readStructFromInput = cli.option("--counterpoise",
         "Read structure from input file (allows for "
         "counterpoise-like corrections of the BSSE)").toBool ();

   double filter = cli.option("--filter","min occ",
                              "set minimum occupation to print")
                   .toDouble (0.01);

   SxString outFile = cli.option("-o|--output","file","(packed) output file")
                      .toString ("");

   bool recompPhaseFactors
      = cli.option ("--recompute-phases",
                    "recompute phase factors: slower, but less memory")
        .toBool ();
   bool averageK
      = cli.option ("--averageK",
                    "average over k-points")
        .toBool ();

   SxString basisFile = cli.option ("--basis", "file", "AOBasis file")
                        .toString ("");

   cli.finalize ();

   SxPtr<SxGkBasis> gkBasisPtr;
   SxFermi fermi;

   FILE *output = (outFile.getSize () > 0) ? sxfopen(outFile, "w") : NULL;
   SxAtomicStructure structure;
   try {
      SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
      gkBasisPtr = SxPtr<SxGkBasis>::create (io, false);
      if (!readStructFromInput) structure.read(io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   SxPW waves(wavesFile, SxPW::ReadOneByOne);
   waves.setGkBasisPtr (gkBasisPtr);
   SxGkBasis &gkBasis = *gkBasisPtr;

   // read pseudopotential data
   cout << "Read PP data" << endl;
   SxParser parser;
   SxParser::Table table = parser.read (inFile);
   if(readStructFromInput) {
      structure = SxAtomicStructure(&*table);
      structure.readElements (&*table);
   }
   const SxArray<SxString> &chemName = structure.getElements ();

   // setup phase factors
   cout << "Set up phase factors" << endl;
   if (recompPhaseFactors)
      for (int ik = 0; ik < gkBasis.getNk (); ++ik)
         gkBasis(ik).memMode = SxGBasis::SaveMemory;
   gkBasis.changeTau(structure);

   cout << "Setting up ao basis" << endl;
   SxPtr<SxAOBasis> aoBasisPtr = aoBasisFromInput (table, basisFile, gkBasis);
   SxAOBasis &aoBasis = *aoBasisPtr;

   // setup radial basis and aoBasis
   aoBasis.setOverlapCaching (SxAOBasis::CacheCurrentK);
   aoBasis.setInvOverlapCaching (SxAOBasis::CacheCurrentK);

   // --- setup map io -> isal
   int nOrb = aoBasis.getNOrb ();
   SxArray<int> map(nOrb); // map orbital to m-summed orbitals
   int nSal = 0;
   SxArray<SxArray<char> > lMu(structure.getNSpecies ());
   SxString lNames="spdfghijklm";
   for (int io = 0; io < nOrb; ++io)  {
     SxAOBasis::OrbitalIndex &sao = aoBasis.orbitalMap(io);
     map(io) = nSal;
     //cout << io << "-> " << isal << endl;
     int l = aoBasis.refOrbMap(sao.is)(sao.io).l;
     int m = aoBasis.refOrbMap(sao.is)(sao.io).m;
     if (l == m) {
        nSal++;
        if (sao.ia == 0)  {
           ssize_t il = lMu(sao.is).getSize ();
           lMu(sao.is).resize (il + 1, true);
           lMu(sao.is)(il) = lNames(l);
        }
     }
   }

   int iSpin = 0, nk = gkBasis.getNk ();
   int nStates = waves.getNStates ();
   PsiG psi, psiAO, allPsi, allPsiAO;
   SxVector<double> normSal(nSal);

   int idxK, ik, idxState, iState;

   if (kPoints.getSize () == 0)
      for (ik = 0; ik < nk; ++ik) kPoints << ik;
   if (states.getSize () == 0)
      for (iState = 0; iState < nStates; ++iState) states << iState;

   if (!averageK)  {
      for (idxK = 0; idxK < kPoints.getSize (); ++idxK)  {
         ik = kPoints(idxK);
         cout << "ik = " << (ik+1) << endl;
         allPsi.reformat(gkBasis(ik).ng, states.getSize ());
         allPsi.setBasis(&gkBasis(ik));
         for (idxState = 0; idxState < states.getSize (); ++idxState)
            allPsi.colRef(idxState) <<= waves(states(idxState),iSpin,ik);
         allPsiAO = (aoBasis | allPsi);
         allPsiAO.reshape(nOrb, states.getSize ());

         for (idxState = 0; idxState < states.getSize (); ++idxState)  {
            iState = states(idxState);
            cout << "i = " << (iState + 1) << endl;
            if (output) fprintf(output,"%u",iState+1);
            // psi = waves(iState,iSpin,ik);
            // psiAO = (aoBasis | psi);
            psiAO = allPsiAO.colRef(idxState);
            normSal.set (0.);
            {
               PsiG normOrb
                  = psiAO.conj () * (aoBasis.getInverseOverlap (ik) ^ psiAO);
               cout << normOrb.real () << endl;
               for (int io = 0; io < nOrb; ++io)
                  normSal(map(io)) += normOrb(io).re; // imaginary parts sum to 0
            }
            // --- output
            for (int is = 0, iSal = 0; is < structure.getNSpecies (); ++is)  {
               for (int ia = 0; ia < structure.getNAtoms(is); ++ia)  {
                  for (int il = 0; il < lMu(is).getSize (); ++il, ++iSal)  {
                     double shellNorm = normSal(iSal);
                     if (output) fprintf(output,"\t%8.6f",shellNorm);
                     if (fabs(shellNorm) >= filter)  {
                        cout << chemName(is) << (ia + 1);
                        cout << ' ' << lMu(is)(il) << "-shell: "
                             << shellNorm << endl;
                     }
                  }
               }
            }
            cout << "Total AO norm: " << normSal.sum () << endl;
            if (output) fprintf(output,"\n");
         }
      }
   } else  {
      for (idxState = 0; idxState < states.getSize (); ++idxState)  {
         iState = states(idxState);
         cout << "i = " << (iState + 1) << endl;
         if (output) fprintf(output,"%u",iState+1);
         normSal.set (0.);
         for (idxK = 0; idxK < kPoints.getSize (); ++idxK)  {
            ik = kPoints(idxK);
            psi = waves(states(idxState),iSpin,ik);
            psiAO = (aoBasis | psi);
            SxVector<PrecCoeffG> normOrb
               = psiAO.conj () * (aoBasis.getInverseOverlap (ik) ^ psiAO);
            for (int io = 0; io < nOrb; ++io)
               normSal(map(io)) += gkBasis.weights(ik) * normOrb(io).re;
         }
         // --- output
         for (int is = 0, iSal = 0; is < structure.getNSpecies (); ++is)  {
            for (int ia = 0; ia < structure.getNAtoms(is); ++ia)  {
               for (int il = 0; il < lMu(is).getSize (); ++il, ++iSal)  {
                  double shellNorm = normSal(iSal);
                  if (output) fprintf(output,"\t%8.6f",shellNorm);
                  if (fabs(shellNorm) >= filter)  {
                     cout << chemName(is) << (ia + 1);
                     cout << ' ' << lMu(is)(il) << "-shell: "
                          << shellNorm << endl;
                  }
               }
            }
         }
         cout << "Total AO norm: " << normSal.sum () << endl;
         if (output) fprintf(output,"\n");
      }
   }

   if (output) fclose(output);

   return 0;

}

