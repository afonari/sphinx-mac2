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

#include <SxHubbardU.h>
#include <SxHamiltonian.h>
#include <SxPAWPot.h>
#include <SxSimpleParser.h>
#include <SxHubbardMO.h>
#include <SxEigensystem.h>
#include <SxIO.h>

// Reference 1: SPHInX Hubbard U implementation notes
// Reference 2: M. Cococcioni, S. de Gironcoli, Phys. Rev. B 71, 035105 (2005)

SxHubbardU::SxHubbardU ()
   : energy(0.), eDoubleCounting (0.), verbose(false)
{
   // empty
}

SxVector<SxComplex16>
SxHubbardU::computeIncr (int iSite, const SxVecRef<SxComplex16> &nij)
{
   SX_CHECK (nij.getNCols () == nij.getNRows (),
             nij.getNCols (), nij.getNRows ());
   int nOrbital = (int)nij.getNRows ();
   SX_CHECK(nOrbital > 0, nOrbital);
   SxVector<SxComplex16> res;
   res.reformat (nOrbital, nOrbital);

   SX_VALIDATE_VECTOR(nij);
   //SxVector<SxComplex16>::Eigensystem eig
   //   = SxVector<SxComplex16>(nij).eigensystem ();
   //cout << "nij ev=" << eig.vals << endl;

   double alpha = (alphaBias.getSize () > 0) ? alphaBias(iSite) : 0.;

   //cout << "alpha=" << alpha << endl;
   SX_LOOP(i)  {
      // ref 1, Eq. (1), ref. 2, Eq. (9)  (diagonal-only part)
      energy += 0.5 * Ueff(iSite) * nij(i,i).re;
      SX_LOOP(j)  {
         // ref 1, Eq. (1), ref. 2, Eq. (9)  (full matrix part)
         energy -= 0.5 * Ueff(iSite) * (nij(i,j) * nij(j,i)).re;

         // ref 1, Eq. (9)  (full matrix part)
         res(i,j) = -0.5 * Ueff(iSite) * (nij(i,j) + nij(j,i)).re;

         // ref 1, Eq. (13)
         eDoubleCounting += 0.5 * Ueff(iSite) * (nij(i,j) * nij(i,j)).re;
      }
      // ref 1, Eq. (9)  (diagonal-only part)
      res(i,i) += 0.5 * Ueff(iSite);
      // alpha contribution to energy ref.2 , Eq. (17)
      energy += alpha * nij(i,i).re;
      // ref 2, Eq. (18)
      res(i,i) += alpha;
   }
   siteOccupation(nij.auxData.iSpin, iSite) = nij.real ();
   if (nij.imag ().norm () > 1e-6 * nij.real ().norm ())
      cout << "nij is complex for site " << (iSite + 1) << endl;
   //cout << (eig.vecs.adjoint () ^ res ^  eig.vecs) * HA2EV << endl;
   //sxprintf ("Hubbard site occupation(%d)=%.8f\n", iSite+1, siteOccupation);
   return res;
}

void SxHubbardU::resize (ssize_t nSite)
{
   Ueff.resize (nSite, true);
   if (alphaBias.getSize () > 0)
      alphaBias.resize (nSite, true);
   if (nSite == 0)  {
      atomicSiteMap.resize (0);
   }
}

void SxHubbardU::setZero ()
{
   energy = eDoubleCounting = 0.;
   SX_LOOP2(iSpin,iSite)
      siteOccupation(iSpin,iSite).resize (0);
}

void SxHubbardU::syncMPI ()
{
#ifdef USE_LOOPMPI
   int nSpin = siteOccupation.getNSpin ();
   SxVector<int> haveOcc(nSpin + 1, siteOccupation.getNSite ());
   SX_LOOP(iSite)  {
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)
         haveOcc(iSpin,iSite) = siteOccupation(iSpin,iSite).getSize () > 0 ? 1 : 0;
      haveOcc(nSpin,iSite) = (int)siteOccupation(0,iSite).getNRows ();
   }
   SxLoopMPI::sum (haveOcc);
   // index nSpin contains number of orbitals
   SX_LOOP(iSite) haveOcc(nSpin,iSite) /= haveOcc(0,iSite);
   SX_LOOP2(iSpin,iSite)  {
      int nOrb = haveOcc(nSpin,iSite);
      if (siteOccupation(iSpin,iSite).getNRows () != nOrb)  {
         /// we don't have this occupation, resize and set to zero
         SX_CHECK (siteOccupation(iSpin,iSite).getNRows () == 0,
                   siteOccupation(iSpin,iSite));
         siteOccupation(iSpin,iSite).reformat (nOrb, nOrb);
         siteOccupation(iSpin,iSite).set (0.);
      }
   }
   siteOccupation.sumMPI ();
   SX_LOOP(iSite)  {
      for (int iSpin = 0; iSpin < nSpin; ++iSpin)  {
         siteOccupation(iSpin,iSite) /= double(haveOcc(iSpin,iSite));
      }
   }
#endif
}

double SxHubbardU::getSiteTotal (ssize_t iSite) const
{
   double res=0.;
   // ref. 2, Eq.(7)
   for (int iSpin = 0; iSpin < siteOccupation.getNSpin (); ++iSpin)
      res += siteOccupation(iSpin,iSite).diag ().sum ();
   return res;
}

SxVector<double> SxHubbardU::computeAtom(int iTl, const SxVecRef<double> &Dij,
                                           double fFull)
{
   SX_CHECK (Dij.getNRows () == Dij.getNCols (),
             Dij.getNRows (), Dij.getNCols ());
   ssize_t npl = Dij.getNRows ();
   SxVector<double> Vij (npl, npl);
   Vij.set (0.);
   // --- to take into account non-spinpolarized cases, we do the following:
   // - Dij contains the summation over spins -> scale by 1/2 to get nij
   // - Vij is fine (acts on each spin the same way)
   // - energy must include spin summation. This is done by changing for
   //   each atom first to the single-spin (scale energy by 1/2), do the
   //   incremental compute, and then do the spin summation (multiply by 2).
   //   For the spin-polarized case, fFull is 1, and the scaling does not
   //   do anything
   double fInv = 1. / fFull;

   // rescale energy to single spin
   energy *= fInv;
   eDoubleCounting *= fInv;

   // --- run over orbital types with Hubbard U
   for (int iType = 0; iType < atomicSiteMap(iTl).getSize (); ++iType)  {
      const AtomicSite& site = atomicSiteMap(iTl)(iType);
      // get occupancy matrix
      SxVector<SxComplex16> nik = fInv * site.map->mapPAWtoHubbard (Dij);
      nik.auxData.iSpin = Dij.auxData.iSpin;
      if (verbose) {
         SxSymEigensystem<double> eig(nik);
         cout << "occ eigenspace (site " << (site.iSite+1) << "):" << endl;
         SX_LOOP(mm)
            cout << eig.vals(mm) << ": " << eig.vecs.colRef(mm) << endl;
      }
      // Hubbard U incremental compute and add contrib to PAW Hamiltonian
      site.map->addHubbardToPAW (computeIncr(site.iSite, nik), Vij);
   }
   // rescale energy to full occupation
   energy *= fFull;
   eDoubleCounting *= fFull;
   return Vij;
}

static void inputProblem (const SxSymbolTable *table)
{
   cout << "Hubbard U in " 
        << table->parserFilename
        << " line " << table->parserLineNumber
        << ":";
}

namespace {
   /// Class that implements DFT+U for single valence projector
   class SingleProjU : public SxHubbardU::HubbardPAWMapper
   {
      public:
         /** \brief Where m=-l is located within the intra-atomic
                    all-projector index (ipl)
             @note same as SxPAWPot.offset
         */
         int offset;
         /// Inverse Norm of the projector ref 1, Eq. (7)
         double Pnl;

         /// Constructor
         SingleProjU (int lIn, int offsetIn, double PnlIn = 1.)
            : offset(offsetIn), Pnl(PnlIn)
         { l = lIn; }

         /// Map the PAW density matrix Dij to the occupancy matrix
         virtual
         SxVector<SxComplex16> mapPAWtoHubbard (const SxVecRef<double> &Dij)
         {
            int nOrb = 2 * l + 1;
            SX_CHECK(nOrb > 0, nOrb, l);
            SxVector<SxComplex16> nik (nOrb, nOrb);
            // extract the occupation matrix for this orbital
            // Ref. 1, Eq. (6)
            SX_LOOP2(i,j)
               nik(i,j) = Dij(offset + i, offset + j) * Pnl;
            return nik;
         }

         /// Add the Hubbard orbital Hamiltonian to PAW projectors
         virtual void addHubbardToPAW (const SxVecRef<SxComplex16> &Hij,
                                       SxVecRef<double> &Vij)
         {
            SX_CHECK(Hij.imag ().normSqr () < 1e-20);
            // Ref. 1, Eq. (12)
            SX_LOOP2(i,j)
               Vij(offset + i, offset + j) += Hij(i,j).re * Pnl;
         }
   };
   /// Class that implements DFT+U for PAW sphere projection
   class ProjPAWSphere : public SxHubbardU::HubbardPAWMapper
   {
      public:
         /// species
         int is;

         /// PAW projector
         SxConstPtr<SxPAWPot> pawPot;

         /// Constructor
         ProjPAWSphere (int iSpecies, int lIn, const SxConstPtr<SxPAWPot> &pot)
            : is(iSpecies), pawPot(pot) { l = lIn; }

         /// Map the PAW density matrix Dij to the occupancy matrix
         virtual
         SxVector<SxComplex16> mapPAWtoHubbard (const SxVecRef<double> &Dij)
         {
            int nOrb = 2 * l + 1;
            SX_CHECK(nOrb > 0, nOrb, l);
            SxVector<SxComplex16> nik (nOrb, nOrb);
            nik.set (0.);
            // get occupancy matrix for this l-channel within PAW sphere
            SX_LOOP(ipt) if (pawPot->lPhi(is)(ipt) == l) {
               SX_LOOP(jpt) if (pawPot->lPhi(is)(jpt) == l) {
                  // <phi_i | PAW sphere | phi_j>
                  double omegaIJ = pawPot->omegaPAW(is)(ipt,jpt);
                  ssize_t offset = Dij.getIdx (pawPot->offset(is)(ipt),
                                               pawPot->offset(is)(jpt));
                  const SxVecRef<double,SubMatrix> &subDij
                     = Dij.getRef<SubMatrix> (offset, nOrb, nOrb);
                  nik.plus_assign_ax (omegaIJ, subDij);
               }
            }
            return nik;
         }

         /// Add the Hubbard orbital Hamiltonian to PAW projectors
         virtual void addHubbardToPAW (const SxVecRef<SxComplex16> &Hij,
                                       SxVecRef<double> &Vij)
         {
            SX_CHECK(Hij.imag ().normSqr () < 1e-20);
            int nOrb = 2 * l + 1;
            SX_LOOP(ipt) if (pawPot->lPhi(is)(ipt) == l) {
               SX_LOOP(jpt) if (pawPot->lPhi(is)(jpt) == l) {
                  double omegaIJ = pawPot->omegaPAW(is)(ipt,jpt);
                  ssize_t offset = Vij.getIdx (pawPot->offset(is)(ipt),
                                               pawPot->offset(is)(jpt));
                  Vij.getRef<SubMatrix> (offset, nOrb, nOrb)
                     .plus_assign_ax (omegaIJ, Hij.real ());
               }
            }
         }
   };
}

void SxHubbardU::read (const SxSymbolTable *table,
                       const SxPtr<SxPAWPot> &pawPotPtr,
                       const SxAtomicStructure &structure)
{
   SX_CHECK (pawPotPtr);
   const SxPAWPot &pawPot = *pawPotPtr;
   SxArray<SxStack<AtomicSite> > sites(structure.getNAtoms ());
   int nSites = (int)getSize ();
#ifndef NDEBUG
   if (nSites > 0)  {
      // --- we override the atomic site map
      SX_CHECK(atomicSiteMap.getSize () == 0, atomicSiteMap.getSize ());
      cout << "In " << __FILE__ << ":" << __LINE__ << endl;
      cout << "WARNING: calling SxHubbardU::read with some previous data"
           << endl;
      // this may be ok if you initialize the atomic sites after other
      // sites. If you need to reinitialize, use resize (0) first
   }
#endif
   SxStack<double> Ustack;
   SxStack<double> bias;
   bool useBias = alphaBias.getSize () > 0;
   if (useBias)  {
      SX_LOOP(i) bias << alphaBias(i);
   }
   // for overriding (updating) settings
   SxList<int> updateId;
   SxList<double> updateU, updateBias;

   SYMBOLPARSE(table)  {
      SYMBOLGROUP("HubbardU")  {
         verbose = SYMBOLGET("verbose").toBool ();
         SxPtr<HubbardPAWMapper> map;
         FOREACH_SYMBOLGROUP("site")  {
            double Usite = SYMBOLGET("U");
            Usite /= HA2EV; // units is eV
            bool update  = SYMBOLGET("update").toBool ();

            int is = -1;
            SxString label;

            // --- get species
            if (HAVE_SYMBOL("species"))  {
               is = SYMBOLGET("species");
               is--;
            } else if (HAVE_SYMBOL("element"))  {
               is = pawPot.find (SYMBOLGET("element"));
               if (is < 0)  {
                  cout << "Hubbard U: element '" 
                       << SYMBOLGET("element")->toString ()
                       << "' not found." << endl;
               }
            } else if (HAVE_SYMBOL("label"))  {
               label = SYMBOLGET("label");
               if (!structure.hasLabels ())  {
                  inputProblem (SYMBOLGROUP_TABLE);
                  cout << "No labels defined in structure" << endl;
                  SX_QUIT;
               }
               // --- find species id for this label (and check that there's
               //     a unique one)
               const SxArray<SxString> &labels = structure.getLabels ();
               SX_LOOP(iTl)  {
                  if (labels(iTl) == label)  {
                     if (is == -1)
                        is = structure.getISpecies ((int)iTl);
                     else if (structure.getISpecies ((int)iTl) != is)  {
                        inputProblem (SYMBOLGROUP_TABLE);
                        cout << "Multiple species for same label:"  << endl;
                        SX_LOOP(jTl) if (labels(jTl) == label)
                           cout << "Atom " << (jTl + 1) << ": species "
                                << (structure.getISpecies ((int)jTl) + 1) 
                                << endl;
                        SX_QUIT;
                     }
                  }
               }
            }
            // 
            double alpha = SYMBOLGET("shift") || 0.;
            alpha /= HA2EV;
            if (fabs(alpha) > 1e-10)  {
               if (!useBias)
                  for (int i = 0; i < nSites; ++i)
                     bias << 0.;
               useBias = true;
            }
            if (is < 0)  {
               inputProblem (SYMBOLGROUP_TABLE);
               cout << " incomplete site specification" << endl;
               cout << "Use species, element, or label" << endl;
            }

            SxString siteSpec;

            int l = SYMBOLGET("l") || -1;

            // --- check ipt 
            if (l >= 0)  {
               if (l > pawPot.lMax(is))  {
                  cout << "Illegal projector type l = " << l
                       << " for species " << (is + 1)
                       << ". Valid range is 0.." << pawPot.lMax(is) << endl;
                  SX_QUIT;
               }
               map = SxPtr<ProjPAWSphere>::create(is, l, pawPotPtr);
               if (verbose)  {
                  siteSpec = SxString::sprintf("PAW sphere (l=%d)", l);
               }
            } else {
               int ipt = SYMBOLGET("projectorType");
               ipt--;
               if (ipt < 0 || ipt >= pawPot.getNProjType (is))  {
                  inputProblem (SYMBOLGROUP_TABLE);
                  cout << "Illegal projector type " << (ipt + 1)
                       << " for species " << (is + 1)
                       << ". Valid range is 1.."
                       << pawPot.getNProjType(is) << endl;
                  SX_QUIT;
               }
               l = pawPot.lPhi(is)(ipt);
               int offset = pawPot.offset(is)(ipt);
               // Ref. 1, Eq. (7)
               //double Pnl = 1./tr(pawPot.pPS(is).colRef(ipt).sqr ());
               double Pnl = 1.;

               if (verbose)  {
                  cout << "Pnl=" << Pnl << endl;
                  siteSpec = SxString::sprintf("projector %d (l=%d)", ipt+1, l);
               }
               if (!update)
                  map = SxPtr<SingleProjU>::create(l, offset, Pnl);
            }

            const SxArray<SxString> *labels = NULL;
            if (label.getSize () > 0) labels = &structure.getLabels ();
            SX_LOOP(iAtom)  {
               ssize_t iTl = structure.getIAtom (is, iAtom);
               // in case labels are used: check if label matches 
               if (labels && (*labels)(iTl) != label) continue;
               if (update)  {
                  if (sites(iTl).isEmpty () || sites(iTl).top ().map != map)
                  {
                     inputProblem (SYMBOLGROUP_TABLE);
                     cout << endl;
                     cout << "'update' flag is only allowed directly after a "
                             " generic setup for the same projector type."
                          << endl;
                     SX_QUIT
                  }
                  updateId << sites(iTl).top ().iSite;
                  updateU << Usite;
                  updateBias << alpha;
               } else {
                  // add new atomic site
                  sites(iTl) << AtomicSite({nSites++, map});
                  Ustack     << Usite;
                  if (useBias) bias << alpha;
               }
               if (verbose)  {
                  cout << "Site " << (sites(iTl).top ().iSite + 1) << ": "
                       << siteSpec << " at atom " << (iAtom+1)
                       << " of species " << (is+1)
                       << " (total atom id: " << (iTl+1) << ") with U="
                       << (Usite * HA2EV) << " eV" << endl;
               }
            }
         }
         FOREACH_SYMBOLGROUP("MO")  {
            // create new MO site
            ssize_t iHubMO = moSite.getSize ();
            moSite.resize (iHubMO + 1, true);
            moSite(iHubMO) = SxPtr<SxHubbardMO>::create ((int)nSites);
            // read
            moSite(iHubMO)->read (SYMBOLGROUP_TABLE, structure, pawPotPtr);

            // read U and potential bias
            double Usite = SYMBOLGET("U");
            Usite /= HA2EV; // units is eV
            double alpha = SYMBOLGET("shift") || 0.;
            alpha /= HA2EV;
            if (fabs(alpha) > 1e-10)  {
               if (!useBias)
                  for (int i = 0; i < nSites; ++i)
                     bias << 0.;
               useBias = true;
            }

            // add U's to complete list of U's
            int nSitesMO = moSite(iHubMO)->getNSite ();
            for (int iSite = 0; iSite < nSitesMO; ++iSite)  {
               Ustack << Usite;
               if (useBias) bias << alpha;
            }

            nSites += nSitesMO;
         }
         FOREACH_SYMBOLGROUP("AO")  {
            // create new MO site
            ssize_t iHubMO = moSite.getSize ();
            moSite.resize (iHubMO + 1, true);
            moSite(iHubMO) = SxPtr<SxHubbardMO>::create ((int)nSites);
            // read
            moSite(iHubMO)->readAO (SYMBOLGROUP_TABLE, structure, pawPotPtr);

            // read U and potential bias
            double Usite = SYMBOLGET("U");
            Usite /= HA2EV; // units is eV
            double alpha = SYMBOLGET("shift") || 0.;
            alpha /= HA2EV;
            if (fabs(alpha) > 1e-10)  {
               if (!useBias)
                  for (int i = 0; i < nSites; ++i)
                     bias << 0.;
               useBias = true;
            }

            // add U's to complete list of U's
            int nSitesMO = moSite(iHubMO)->getNSite ();
            for (int iSite = 0; iSite < nSitesMO; ++iSite)  {
               Ustack << Usite;
               if (useBias) bias << alpha;
            }

            nSites += nSitesMO;
         }
      }
   }
   resize (nSites);
   atomicSiteMap.resize (structure.getNAtoms ());
   SX_LOOP(iTl) atomicSiteMap(iTl) = sites(iTl);
   // fill in effective U from stack from the end (to keep the
   // previous ones, if any)
   while (!Ustack.isEmpty ()) Ueff(--nSites) = Ustack.pop ();

   // set defined alphas
   if (useBias) alphaBias = bias;

   // implement updates (requires random access)
   while (updateId.getSize () > 0)  {
      int iSite = updateId(0);
      Ueff(iSite) = updateU(0);
      if (useBias) alphaBias(iSite) = updateBias(0);
      updateId.removeFirst ();
      updateU.removeFirst ();
      updateBias.removeFirst ();
   }
   siteOccupation.resize (SxHamiltonian::getNSpin (table->topLevel ()), (int)getSize ());
}

void SxHubbardU::computeMO (const SxArray<SxPtr<SxBlockDensityMatrix> > &Pij,
                            const SxAtomicStructure &structure)
{
   for (int i = 0; i <moSite.getSize (); i++)  {
      moSite(i)->compute (this, *Pij(i), structure);
   }
}

void SxHubbardU::writeOccupations(const SxString &fileName) const
{
   {
      // --- check if there is anything available
      ssize_t nData = 0;
      SX_LOOP2(iSpin,iSite) nData += siteOccupation(iSpin,iSite).getSize ();
      if (nData == 0) {
         cout << "No Hubbard occupations to write" << endl;
         return;
      }
   }
   SX_CHECK (getSize () > 0, getSize ());
   FILE *fp = sxfopen(fileName,"w");
   SX_LOOP(iSite)  {
      SX_LOOP(iSpin)  {
         SxSymEigensystem<double> eig(siteOccupation(iSpin,iSite));
         SX_LOOP (iVal)  {
            fprintf (fp, "%8.6f   ", eig.vals(iVal));
            SX_LOOP(i) fprintf (fp, "\t% 8.6f", eig.vecs(i,iVal));
            fprintf (fp, "\n");
         }
      }
      fprintf (fp,"\n");
   }
   fclose (fp);
}
