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
#include <SxAOBasis.h>
#include <SxFermi.h>

SxVector<SxComplex16>
getDensityMatrix (const SxAOBasis          &ao,
                  const PsiG               &waves,
                  const Focc               &focc)
{
   int iSpin = waves.auxData.iSpin;
   int ik = waves.auxData.ik;

   int nStates = int(waves.getNCols ());
   if (nStates == 0) nStates = 1;

   // --- get number of occupied states
   int nOccStates = 0;
   PrecFocc lastFocc = focc(0,iSpin,ik);
   bool sorted = true;
   for (int i = 0; i < nStates; ++i)  {
      if (focc(i,iSpin,ik) > lastFocc + 1e-6) sorted = false;
      if (fabs(lastFocc = focc(i,iSpin,ik)) > 1e-5) ++nOccStates;
   }
   SX_CHECK (nOccStates > 0);
   
   // --- get waves/focc for occupied states
   SxVecRef<PrecCoeffG> occWaves;
   SxVecRef<PrecFocc> myFocc;
   int ng = int(waves.getNRows ());
   if (sorted)  {
      // standard case: packed storage
      occWaves = waves(SxIdx(0, ng * nOccStates - 1));
      occWaves.reshape (ng, nOccStates);
      occWaves.auxData = waves.auxData;
      myFocc = focc(iSpin,ik)(SxIdx(0,nOccStates-1));
   } else {
      if (2 * nOccStates < nStates)  {
         // pack occupied states
         occWaves = PsiG(ng, nOccStates);
         occWaves.auxData = waves.auxData;
         myFocc = SxVector<PrecFocc>(nOccStates);
         int iState = 0;
         for (int i = 0; i < nStates; ++i)  {
            if (focc(i,iSpin,ik) > 1e-5)  {
               occWaves.colRef(iState) <<= waves.colRef(i);
               myFocc(iState++) = focc(i,iSpin,ik);
            }
         }
      } else {
         // do not care about the unnecessary states
         occWaves = waves;
         nOccStates = nStates;
         myFocc = focc(iSpin,ik);
      }
   }
   
   // get (S^-1) <mu|psi>
   SxAOBasis::TPsi aoPsi;
   //               S^{-1}           < mu | psi >
   aoPsi = ao.getInverseOverlap(ik) ^ (ao | occWaves);

   // release occWaves memory
   occWaves.unref ();

   // multiply aoPsi with sqrt(focc(iState))
   for (int i = 0; i < nOccStates; ++i)
      aoPsi.colRef (i) *= (myFocc(i) > 0. ? SxComplex16(1.) : I )
                        * sqrt(fabs(myFocc(i)));
   
   // return <mu|psi(i)>focc(i)<psi(i)|mu>
   return (aoPsi ^ aoPsi.adjoint ());
}

int main (int argc, char **argv)
{

   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   SxCLI cli (argc,argv);
   cli.authors = "C. Freysoldt";
   cli.preUsageMessage = "This is the Mulliken analysis tool.";

   SxString wavesFile = cli.option ("-w|--waves","file","waves file")
                        .toString("waves.sxb");

   SxString inFile = cli.option ("-i|--input", "file", "S/PHI/nX input file")
                     .toString ("input.sx");

   enum SpinMode { Alpha, Beta, Rho, Spin } mode = Rho;
   cli.option ("--alpha|--beta|--rho|--spin", "mode for spin");
   cli.last ().defaultValue = "rho";
   if (cli.last ().exists ())  {
      mode = SpinMode (cli.last ().toChoice ());
   };
   bool recompPhaseFactors
      = cli.option ("--recompute-phases",
                    "recompute phase factors: slower, but less memory")
        .toBool ();

   SxString basisFile = cli.option ("--basis", "file", "AOBasis file")
                        .toString ("");

   cli.version ("1.0");
   cli.finalize ();

   SxPtr<SxGkBasis> gkBasisPtr;
   SxFermi fermi;

   SxAtomicStructure structure;

   try {
      SxBinIO io (wavesFile, SxBinIO::BINARY_READ_ONLY);
      gkBasisPtr = SxPtr<SxGkBasis>::create (io, false);
      fermi.read (io);
      fermi.kpPtr = &*gkBasisPtr;
      structure.read(io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   SxPW waves(wavesFile, SxPW::ReadOnDemand);
   waves.setGkBasisPtr (gkBasisPtr);
   SxGkBasis &gkBasis = *gkBasisPtr;

   // read pseudopotential data
   cout << "Read input file" << endl;
   SxParser parser;
   SxParser::Table table = parser.read (inFile);
   SxSpeciesData sData(&*table);

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

   SxAOBasis::TPsi psiAO, psi;
   int iSpin = 0, nk = gkBasis.getNk ();
   int nOrb = aoBasis.getNOrb ();

   SxVector<SxComplex16> P, S;
   SxVector<double> M(nOrb, nOrb);
   M.set (0.);
   cout << "Computing full Mulliken matrix..." << endl;
   for (int ik = 0; ik < nk; ++ik)  {
      cout << "Computing S..." << endl;
      S = aoBasis.getOverlap(ik);
      cout << "Computing P..." << endl;
      for (iSpin = 0; iSpin < waves.getNSpin (); ++iSpin)  {
         if (mode == Alpha && iSpin == 1) continue;
         if (mode == Beta  && iSpin == 0) continue;
         P = getDensityMatrix (aoBasis, waves(iSpin, ik), fermi.focc);
         if (mode == Spin && iSpin == 1)  {
            M.plus_assign_ax(-gkBasis.weights(ik),(S.conj () * P).real ());
         } else {
            M.plus_assign_ax(gkBasis.weights(ik),(S.conj () * P).real ());
         }
      }
   }

   cout << "Condensing full Mulliken matrix to l-channels..." << endl;
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

   SxVector<double> MsumM(nSal, nSal);
   MsumM.set (0.);
   
   {
      int iSal,jSal,iOrb,jOrb;
      SxAOBasis::OrbitalIndex orb;
      for (iOrb = 0; iOrb < nOrb; ++iOrb)  {
         iSal = map(iOrb);
         for (jOrb = 0; jOrb < nOrb; ++jOrb)  {
            jSal = map(jOrb);
            MsumM(iSal,jSal) += M(iOrb,jOrb);
         }
      }
   }

   // --- output
   {
      double totalN = 0.;
      for (int is = 0, iSal = 0; is < structure.getNSpecies (); ++is)  {
         for (int ia = 0; ia < structure.getNAtoms(is); ++ia)  {
            double charge = sData.valenceCharge(is), nMull;
            SX_LOOP(il)  {
               charge -= (nMull = MsumM.colRef(iSal).sum ());
               cout << sData.chemName(is) << (ia + 1);
               cout << ' ' << lMu(is)(il) << "-shell: "
                    << SxString(nMull, "%.3lf") << endl;;
               totalN += nMull;
               iSal++;
            }
            cout << sData.chemName(is) << (ia + 1);
            if (mode == Rho)  {
               cout << SxString(charge, " Mulliken charge: %.3lf") << endl;
            } else {
               cout << " total occupation: " 
                    << SxString(sData.valenceCharge(is) - charge, "%.3lf")
                    << endl;
            }
         }
      }
      cout << "Total nEl " << SxString(totalN, "%.3lf") << endl;
   }

   return 0;

}
