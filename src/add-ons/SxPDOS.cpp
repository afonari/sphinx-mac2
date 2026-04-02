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

#include <SxPDOS.h>
#ifndef SX_STANDALONE
#include <SxSpectrum.h>
#include <SxHamSolver.h>
#include <SxSimpleParser.h>

SxPDOS::SxPDOS (const SxConstPtr<SxPartialWaveBasis> &pBasisIn,
                const SxConstPtr<SxPAWPot> &pawPotIn)
   : pBasisPtr(pBasisIn), pawPotPtr (pawPotIn),
     broadening (0.1 / HA2EV),
     eMin (-1.), eMax(-1.),
     nPerPeak (15),
     verbose (false)
{
   outFile = "pdos";
}

void SxPDOS::compute (const SxPAWSet &waves, const Eps &eps)
{
   const SxPartialWaveBasis &pBasis = *pBasisPtr;
   const SxPAWPot &pawPot = *pawPotPtr;   
   const SxAtomicStructure &structure = waves.getGkBasis ().getTau ();
   const SxVector<double> &weights = waves.getGkBasis ().weights;

   // --- calculate spectrum
   SxSpectrum pdos;
   pdos.nPerPeak = nPerPeak;
   SxComplex16 kWeight, pWeight;

   int nChannels = 1;
   SxArray<int> offset (structure.getNSpecies ());
   if (atoms.getSize () == 0)  {
      SX_LOOP(is) {
         offset(is) = nChannels;
         nChannels += pawPot.lMax(is) + 1;
      }
   } else {
      offset.set (0);
      for (auto iAtomTl : atoms)  {
         int is = structure.getISpecies (iAtomTl);
         if (offset(is) == 0)  {
            cout << "pdos for species " << (is + 1)
                 << " (" << structure.getElements ()(is) 
                 << ") L=0.." << pawPot.lMax(is)
                 << " in columns " << (nChannels + 1);
            offset(is) = nChannels;
            nChannels += pawPot.lMax(is) + 1;
            cout << "..." << nChannels 
                 << endl;

         }
      }
   }

   int nSpin = eps.getNSpin ();
   int nk = eps.getNk ();
   for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
      pdos.set (eMin, eMax, broadening, nChannels);
      int ipBasis = 0;

      // --- check atom number
      int nAtom = structure.getNAtoms ();
      for (int iAtom = 0; iAtom < nAtom; ++iAtom)  {
         int iSpecies = structure.getISpecies (iAtom);
         int npt = pawPot.getNProjType (iSpecies);
         bool lastAtom = (iAtom == nAtom - 1);

         int aWrite = 0;
         double aWeight = 0.;
         if (atoms.getSize () == 0) aWeight = 1.;
         else if (atoms.contains (iAtom))  {
            aWrite = 1;
            aWeight = 1.;
         }

         SxSpectrum pdosA;
         pdosA.nPerPeak = nPerPeak;
         if (aWrite) 
            pdosA.set (eMin, eMax, broadening, pawPot.lMax(iSpecies)+2);

         // ---integration
         if (verbose)
            cout << "Species " << (iSpecies+1) << " has rPAW="
                 << pawPot.rPAW(iSpecies) << endl;

         for (int ipt = 0; ipt < npt; ++ipt)  {
            int li = pawPot.lPhi(iSpecies)(ipt);
            int jpBasis = ipBasis;
            for (int jpt = ipt; jpt < npt; ++jpt)  {
               int lj = pawPot.lPhi(iSpecies)(jpt);
               if (lj == li)  {
                  double integral = pawPot.omegaPAW(iSpecies)(ipt, jpt);
                  if (verbose)
                     cout << "ipt " << ipt << " jpt " << jpt 
                          << " integral=" << integral << endl;

                  // --- add peaks
                  for (int ik = 0; ik < nk; ++ik)  {
                     SX_MPI_LEVEL("waves-k");
                     if (!SxLoopMPI::myWork(ik)) continue;
                     kWeight = 2./double(nSpin) * weights(ik) * aWeight;
                     for (int in = 0; in < eps.getNStates(ik); ++in)  {

                        // --- calculate projected weight
                        SxVector<SxComplex16> projpsi = pBasis | waves (in, iSpin, ik);
                        pWeight = 0.;
                        // --- sum <psi|p_j><p_i|psi> over m (im = m + l)
                        SxComplex16 sum = 0.;
                        for (int im = 0; im < 2 * li + 1; ++im)  {
                           SxComplex16 element = projpsi(jpBasis + im).conj () 
                                               * projpsi(ipBasis + im);
                           sum += element;
                        }
                        if (jpt == ipt) pWeight += integral * sum;
                        else pWeight += integral * (sum + sum.conj());
                        double w = kWeight * pWeight;
                        pdos.addPeak (eps(in,iSpin,ik), w, 0);
                        pdos.addPeak (eps(in,iSpin,ik), w, offset(iSpecies) + li);
                        if (aWrite)  {
                           pdosA.addPeak (eps(in,iSpin,ik), w, 0);
                           pdosA.addPeak (eps(in,iSpin,ik), w, li + 1);
                        }
                     } //in
                  } //ik
                  jpBasis += 2 * lj + 1;
               } //check jpt
            } //jpt
            ipBasis += 2 * li + 1;
         } //ipt
      
         SX_MPI_SOURCE("waves-k", TaskGroupMaster);
         SX_MPI_TARGET(TopLevel, TaskGroupAll);
         // compute spectrum
         if (aWrite)  {
            pdosA.compute ();
            pdosA.spectra  /= (iSpin == 0 ? 1. : -1.) * HA2EV;
            SxLoopMPI::sum (pdosA.spectra);
            pdosA.energies *= HA2EV;
         }

         if (lastAtom)  {
            pdos.compute ();
            pdos.spectra  /= (iSpin == 0 ? 1. : -1.) * HA2EV;
            SxLoopMPI::sum (pdos.spectra);
            pdos.energies *= HA2EV;
         }

         // --- output
         SX_MPI_MASTER_ONLY  {
            SxString ext = ".dat";
            if (nSpin == 2)  {
               if (iSpin == 0) ext = ".up";
               else ext = ".down";
            }

            if (aWrite)
               pdosA.write (outFile + "Atom" + SxString (iAtom) + ext);

            if (lastAtom)
               pdos.write (outFile + "sum" + ext);
         }

      } //atom
   } //spin
}

void SxPDOS::getRange (const Eps &eps, double shift)
{
   eMin = eps(0,0,0);
   eMax = eps(0,0,0);
   double eMinK,eMaxK;
   SX_LOOP2(iSpin, ik) {
      eMinK = eps(iSpin,ik).minval ();
      if (eMinK < eMin) eMin = eMinK;
      eMaxK = eps(iSpin,ik).maxval ();
      if (eMaxK > eMax) eMax = eMaxK;
   }
   eMin -= shift;
   eMax += shift;
}

void computePDOS (SxHamSolver *pot, const SxSymbolTable *table, bool calc)
{
   SxPtr<SxPAWSet> waves = pot->wavesPtr;
   const SxAtomicStructure &structure = waves->getGkBasis ().getTau ();
   SxPDOS pDos(waves->getPBasis (), pot->potPtr);
   SYMBOLPARSE(table) {
      pDos.outFile << SYMBOLGET("file");
      if (pDos.outFile.tail (4) == ".dat")
         pDos.outFile = pDos.outFile.head (pDos.outFile.getSize () - 4);

      // --- atoms
      if (HAVE_SYMBOL("atoms")) {
         pDos.atoms = SYMBOLGET("atoms")->toIntList ();
         // change from 1..N to 0..N-1
         for (auto ia : pDos.atoms) ia--;
      } else if (HAVE_SYMBOL("elements")) {
         SxList<SxString> elements = SYMBOLGET("elements")->toStringList ();
         const SxArray<SxString> &chemNames = structure.getElements ();
         SX_LOOP(is)  {
            if (elements.contains (chemNames(is)))  {
               for (int ia = 0; ia < structure.getNAtoms ((int)is); ++ia)
                  pDos.atoms.append (structure.getIAtom (is, ia));
            }
         }
      } else {
         FOREACH_SYMBOLGROUP("atomRange")  {
            int from = SYMBOLGET("from"),
                to = SYMBOLGET("to");
            for (int ia = from; ia <= to; ++ia)
               pDos.atoms.append (ia-1); // start from 0 rather than from 1
         }
      }
      if (pDos.atoms.getSize () == 0 && ! SYMBOLGET("sumOnly").toBool ())  {
         for (int ia = 0; ia < structure.getNAtoms (); ++ia)
            pDos.atoms.append (ia);
      }

      pDos.nPerPeak << SYMBOLGET("fine");
      pDos.broadening = (SYMBOLGET("broadening") || 0.1) / HA2EV;

      // --- energy range
      if (HAVE_SYMBOL("range"))  {
         SxList<double> range = SYMBOLGET("range")->toList ();
         SX_CHECK(range.getSize () == 2);
         pDos.eMin = range(0) / HA2EV;
         pDos.eMax = range(1) / HA2EV;
      } else {
         pDos.getRange (pot->fermi.eps, 2. * pDos.broadening);
      }
      pDos.verbose = SYMBOLGET("verbose").toBool ();
   }
   if (pDos.verbose)  {
      cout << "Writing PDOS into '" << pDos.outFile << "*.dat'" << endl;
      cout << "broadening=" << pDos.broadening * HA2EV << " eV" << endl;
      cout << "range=" << pDos.eMin * HA2EV << " ... " << pDos.eMax * HA2EV << " eV" << endl;
      if (pDos.atoms.getSize () > 0)  {
         cout << "atoms=";
         for (auto ia : pDos.atoms)
            cout << (ia + 1) << ",";
         cout << endl;
      }
   }
   if (!calc) return;
   pDos.compute (*waves, pot->fermi.eps);
}

#else

#include <SxCLI.h>

int main (int argc, char **argv)
{
   initSPHInXMath ();
   SxLoopMPI::init (argc, argv);
   
   SxCLI cli(argc, argv);

   cli.preUsageMessage 
      = SxString(" Add-on calculates projected density of states for PAW\n"
                 "  Gaussian broadening applied to each state"
                 " The first column of output is the total PDOS, and"
                 " proceeding data columns are PDOS of respective l.").wrap ();
   cli.authors = "Siyuan Zhang and Christoph Freysoldt";
   
   SxString wavesFile
            = cli.option ("-w|--waves","file","input waves file")
             .toString   ("waves.sxb");
   
   SxString inputFile 
      = cli.option("--input","file","read k-vectors from this input file")
        .toString ("input.sx");
   
   SxString outFile
            = cli.option ("-o","file","prefix of the output file;"
                          ".dat is always added to the end of the output file")
             .toString   ("pdos");
   
   SxList<int> atoms = cli.option ("-a|--atoms", "atom list", 
                                   "atoms of interest for PDOS")
                      .toIdxList ();
   cout << "Entered atoms contain number: " << atoms << endl;
   cli.last ().required (false);
   cli.last ().defaultValue 
            = "default returns a total spectrum but no individuals";

   double broadening = cli.option ("-b|--broad", "energy [eV]", 
                                   "Gaussian weight factor, eV")
                       .toDouble(0.1) / HA2EV;

   double shift = cli.option ("-s|--shift", "energy [eV]", 
          "shift borders of energy interval from (energyMin, energyMax) "
          "to (energyMin-shift, energyMax+shift), eV")
                  .toDouble(2. * broadening * HA2EV) / HA2EV;
   cli.last ().defaultValue = "default: 2 * broadening";

   SxList<double> range = cli.option("--range","range","Emin:Emax (eV)")
                          .toDoubleList ();
   if (range.getSize () != 2 && range.getSize () != 0 && !cli.error)  {
      cout << "Illegal range." << endl;
      cli.setError ();
   }
   
   int nPerPeak = cli.option ("--fine", "integer", 
                              "energy resolution parameter, approximately "
                              "number of points used to display 1/2 Gauss peak")
                  .toInt(15);
   
   cli.finalize ();

   // --- read input file
   SxParser parser;
   SxParser::Table table = parser.read(inputFile);

   SxPtr<SxPAWPot> pawPotPtr =SxPtr<SxPAWPot>::create (&*table);
  
   // --- read waves file (not waves yet)
   SxFermi fermi;
   SxAtomicStructure structure;
   SxPtr<SxGkBasis> gkBasisPtr;
   SxPtr<SxPartialWaveBasis> pBasisPtr;
   SxBinIO io;
   try  {
      io.open (wavesFile, SxBinIO::BINARY_READ_ONLY);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   fermi.read (io);

   structure.read (io);
   gkBasisPtr = SxPtr<SxGkBasis>::create (io, false);
   gkBasisPtr->changeTau (structure);

   // create partial wave basis
   pBasisPtr = SxPtr<SxPartialWaveBasis>::create (pawPotPtr, structure);
   pBasisPtr->createProjBasis (*gkBasisPtr);
   // construct PAW waves container
   SxPAWSet waves (gkBasisPtr, pBasisPtr, fermi.getNStates (), fermi.getNSpin ());
   
   waves.read (io, SxPWSet::KeepGkBasis);
   io.close ();

   SxPDOS pDos (pBasisPtr, pawPotPtr);

   pDos.outFile = outFile;
   pDos.atoms = atoms;
   pDos.nPerPeak = nPerPeak;
   //pDos.verbose = true;

   // --- get minimum and maximum
   if (range.getSize () == 0)  {
      pDos.getRange (fermi.eps, shift);
   } else {
      pDos.eMin = range(0) / HA2EV;
      pDos.eMax = range(1) / HA2EV;
   }

   pDos.compute (waves, fermi.eps);

}
#endif
