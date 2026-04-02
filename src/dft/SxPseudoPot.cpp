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
#include <SxPseudoPot.h>
#include <SxConstants.h>
#include <SxRadBasis.h>
#include <SxYlm.h>
#include <SxFileParser.h>

SxPseudoPot::SxPseudoPot ()
{
   nlcc = false;
}


SxPseudoPot::SxPseudoPot (const SxSymbolTable *group)
   : SxSpeciesData ()
{
   int l, iSpecies;
   SxSymbolTable *species;
   SxList<int>    intList;
   SxList<double> realList;
   SxString potType;
   try  {
      if (group->getName () != "pseudoPot")
         group = group->getGroup("pseudoPot");
      SxSpeciesData::readSpecies (group);

      // --- pseudoPot
      // --- pseudoPot.species
      nSpecies = group->getGroup("species")->getNItems ("species");
      lMax.resize           (nSpecies);
      lLoc.resize           (nSpecies);
      realSpace.resize      (nSpecies);
      rGauss.resize         (nSpecies);
      pseudoFocc.resize     (nSpecies);
      foccAtomicRho.resize  (nSpecies);
      pseudoPotFiles.resize (nSpecies);
      for ( species  = group->getGroup  ("species"), iSpecies=0;
            species != NULL;
            species  = species->nextSibling ("species"), iSpecies++ )
      {
         // pseudoPot.species.potential
         pseudoPotFiles(iSpecies) = species->get("potential")->toString();
         // pseudoPot.species.lMax
         lMax(iSpecies) = species->get("lMax")->toInt();
         // pseudoPot.species.lLoc
         lLoc(iSpecies) = species->get("lLoc")->toInt();
         // pseudoPot.species.realSpaceProjectors
         realSpace(iSpecies) = species->contains("realSpaceProjectors")
                             ? species->get("realSpaceProjectors")
                               ->toAttribute ()
                             : false;
         // pseudoPot.species.rGauss
         rGauss(iSpecies) = species->get("rGauss")->toReal();
         // pseudoPot.species.lcaoOrbitals
         intList = species->get("lcaoOrbitals")->toIntList();
         pseudoFocc(iSpecies).resize (lMax(iSpecies)+1);
         for (l=0; l <= lMax(iSpecies); l++)
            pseudoFocc(iSpecies)(l) = intList.contains(l);
         if (species->contains("potType")) potType = species->get("potType")->toString();
         else potType = "FHI";
         // pseudoPot.species.atomicRhoOcc
         if (species->contains("atomicRhoOcc"))  {
            SxSymbol *sym = species->get("atomicRhoOcc");
            int ind, iSpin, nSpin = sym->getRank();
            if (nSpin == 0)  nSpin = 1; // for hydrogen/helium
            foccAtomicRho(iSpecies).resize(lMax(iSpecies)+1);
            for (l=0; l <= lMax(iSpecies); l++)
               foccAtomicRho(iSpecies)(l).resize (nSpin);
            realList  = species->get("atomicRhoOcc")->toList();
            if (realList.getSize() != nSpin*(lMax(iSpecies)+1))  {
               cout << "In species ";
               if (elementName(iSpecies).getSize () > 0)
                  cout << elementName(iSpecies);
               else if (chemName(iSpecies).getSize () > 0)
                  cout << chemName(iSpecies);
               else cout << (iSpecies+1);
               sxprintf (": %d <--> %d\n", (int)realList.getSize(),
                         nSpin*(lMax(iSpecies)+1));
               cout << "Corrupt atomicRhoOcc!" << endl;
               SX_QUIT;
            }
            for (iSpin = 0; iSpin < nSpin; iSpin++)  {
               for (l=0; l <= lMax(iSpecies); l++)  {
                  ind = l + iSpin*(lMax(iSpecies)+1);
                  foccAtomicRho(iSpecies)(l)(iSpin) = realList(ind);
               }
            }
         }
      }

   }  catch (SxException e) {
      e.print();
      SX_EXIT;
   }

   SxArray<SxString> psiFiles(nSpecies); // missing!
   if (potType == "FHI") initPseudoPotentials ();
   else if (potType == "Siesta") initPseudoPotSiesta (psiFiles);
   else {
      cout << "Unknown potential type, S/PHI/nX quits here!" << endl;
      SX_QUIT;
   }
}   

SxPseudoPot::~SxPseudoPot ()
{
   // empty
}


void SxPseudoPot::initPseudoPotentials ()
{
   double zv, radf, psi, pot, rhoCore, lDr;
   int i, j, n, l, lMaxIdx, lLocIdx, iSpecies;
   int iPoint;
   Real8 norm;

   SxList<Real8> radList, psiList, potList;  // :r
   nlcc = false;

   rad.resize (nSpecies);
   pseudoPsi.resize (nSpecies);
   pseudoPot.resize (nSpecies);
   radCore.resize   (nSpecies);
   logDr.resize   (nSpecies);
   for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      radCore(iSpecies).resize(2);
      SxFileParser fp (pseudoPotFiles (iSpecies));

      // --- read valence charge
      fp.topic ("valence charge and max l");
      fp >> zv >> lMaxIdx;
      if (lMaxIdx > lMax(iSpecies)+1)  {
         cout << "WARNING: omitting l > " << lMax(iSpecies) 
              << " from pseudopotential file '" << pseudoPotFiles (iSpecies)
              << "'." << endl;
      } else if (lMaxIdx <= lMax(iSpecies))  {
         cout << "Requested l=0.." << lMax(iSpecies) << " from '"
              << pseudoPotFiles (iSpecies) << "', but it contains only "
              << lMaxIdx
              << " l-channels.\nCheck the pseudopotential metadata!\n";
         SX_QUIT;
      }
      fp.nextLine ();

      // --- skip following 10 lines
      for (i=0; i<10; i++)  fp.nextLine ();

      pseudoPsi(iSpecies).resize (lMax(iSpecies) + 1);
      pseudoPot(iSpecies).resize (lMax(iSpecies) + 1);
      // --- read orbitals (2nd column=rad, 3th column=psi,
      //                    4th column=potential)
      for (l=0; l < lMaxIdx; l++)  {
         // --- number of entries of next orbital
         fp.topic ("number of entries of orbital l=" + l);
         fp >> n >> lDr;
         lDr = log(lDr);
         if (l <= lMax(iSpecies))
            logDr(iSpecies) = lDr; // radial mesh is checked below
         fp.nextLine ();
         fp.topic ("mesh data for l=" + l);
         for (i=0; i<n; i++)  {
            fp >> iPoint >> radf >> psi >> pot;
            radList.append (radf);
            psiList.append (psi/radf);
            potList.append (pot);
         }
         norm = sqrt(( SxVector<double> (radList).cub()
                     * SxVector<double> (psiList).sqr()).integrate(lDr)
                    );
////         for (j=0; j<n; j++)  psiList(j) *= exp (-rad/rc)^2;
         if (fabs (norm) > 1e-10)
            for (j=0; j<n; j++)  psiList(j) /= norm;

         if (l <= lMax(iSpecies))  {
            if (rad(iSpecies).getSize () == 0)
               rad(iSpecies) = radList;
            else if (rad(iSpecies).getSize () != radList.getSize ()
                     || (rad(iSpecies) - SxVector<double>(radList))
                        .normSqr () > 1e-12 * (double)radList.getSize ())  {
               fp.where ();
               cout << endl << "Varying radial meshes are not supported within "
                       "a species." << endl;
               cout << "If necessary, interpolate the offending channel l="
                    << l << endl;
               SX_QUIT;
            }

            pseudoPsi(iSpecies)(l) = psiList;
            // Set up auxData
            SxAuxData &auxData = pseudoPsi(iSpecies)(l).auxData;
            auxData.is = iSpecies;
            auxData.n = 0;
            auxData.l = char(l);
            pseudoPot(iSpecies)(l) = potList;
         }

         radList.removeAll ();
         psiList.removeAll ();
         potList.removeAll ();
      }

      // --- NLCC
      int ir = 0;
      fp.skipWhite ();
      fp.topic ("pseudo-core density");
      while ( !feof (fp.fp) )  {
         // TODO: old/new pseudopotentials (2/4 columns)
         fp >> radf >> rhoCore;
         radCore(iSpecies)(0).append (radf);
         radCore(iSpecies)(1).append (rhoCore);
         if (fabs(radf-rad(iSpecies)(ir)) > 1e-12 * fabs(radf))  {
            fp.where ();
            cout << ": Pseudo-core density must have same mesh as psi's.\n";
            SX_QUIT;
         }
         ir++;
         // skip to end of line
         while (fgetc (fp.fp) != '\n' && !feof(fp.fp)) {}
         fp.skipWhite ();
      }
      if (radCore(iSpecies)(0).getSize() > 0)
         nlcc = true;

   }

   // --- substract local part from pseudopotential
   for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      lLocIdx = lLoc(iSpecies);
      for (l=0; l<=lMax(iSpecies); l++)  {
         if (l != lLocIdx )  {
            pseudoPot(iSpecies)(l) -=  pseudoPot(iSpecies)(lLocIdx);
         }
      } // l
   } // iSpecies
}

void SxPseudoPot::initPseudoPotSiesta (SxArray<SxString> /*psiFiles*/)
{
   SX_EXIT; // not tested
   // TODO: separate pseudopotential reader...

   rad.resize (nSpecies);
   pseudoPsi.resize (nSpecies);
   pseudoPot.resize (nSpecies);
   radCore.resize   (nSpecies);
   logDr.resize   (nSpecies);
   for (int iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      radCore(iSpecies).resize(2);
      SxFileParser fp (pseudoPotFiles (iSpecies));

      // --- read header
      fp.topic ("header");
      SxString element = fp.getItem (),
               XC = fp.getItem (),
               rel = fp.getItem (),
               pseudoCore = fp.getItem ();

      if (!(element = elementName(iSpecies))) {
         cout << "Wrong element! Is " << element << "should be " << elementName(iSpecies) << endl;
         SX_QUIT;
      }
      fp.nextLine (2);
      int nDown, nUp, nGrid;
      double r0, zVal;
      fp >> nDown >> nUp >> nGrid >> r0 >> logDr(iSpecies) >> zVal;
      if (nUp > 0) {
         cout << "Spin in pseudoPotential not yet implemented!" << endl;
         SX_EXIT;
      }

      fp.nextLine ();
      // read radial Basis
      fp.topic ("radial basis");
      rad(iSpecies) = fp.getVector (nGrid);

      // read PseudoPotential
      fp.topic ("pseudo potential");
      pseudoPot(iSpecies).resize(nDown);
      for (int iDown = 0; iDown < nDown; iDown++)  {
         fp.nextLine ();
         int l = fp.getInt ();
         pseudoPot(iSpecies)(l) = fp.getVector (nGrid);
      }

      // read core Charge
      /*
      radCore(iSpecies)(0) = 1.0 * rad(iSpecies);
      radCore(iSpecies)(1).resize(nGrid);
      counter = 0;
      nlcc = true;
      SxVector<double> &rc = radCore(iSpecies);
      for (int iLine = 0; iLine < lines; iLine++,counter+=4)  {
         nVarsRead = fscanf (fp, "%le %le %le %le", &rc(1)(counter), &rc(1)(counter+1), &rc(1)(counter+2), &rc(1)(counter+3));
      }
      switch (lastEntrys)  {
         case 1:
            nVarsRead = fscanf (fp, "%le", &rc(counter));
            break;
         case 2:
            nVarsRead = fscanf (fp, "%le %le", &rc(counter), &rc(counter+1));
            break;
         case 3:
            nVarsRead = fscanf (fp, "%le %le %le", &rc(counter), &rc(counter+1), &rc(counter+2));
            break;
         default:
            break;
      }
      */

      pseudoPsi(iSpecies).resize(nDown);
      // read pseudoPsi missing
      SX_EXIT;

   }

   // --- substract local part from pseudopotential
   for (int iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      int lLocIdx = lLoc(iSpecies);
      for (int l = 0; l <= lMax(iSpecies); l++)  {
         if (l != lLocIdx )  {
            pseudoPot(iSpecies)(l) -=  pseudoPot(iSpecies)(lLocIdx);
         }
      } // l
   } // iSpecies
}

void SxPseudoPot::print () const
{
   int iSpecies;
   Real8 lDr;
   for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      sxprintf ("|    %-15s  ", elementName(iSpecies).ascii());
      cout << "rGauss=" << rGauss(iSpecies)
           << ", zv="   << valenceCharge(iSpecies)
           << ", loc="  << lLoc(iSpecies);
      if (radCore(iSpecies)(1).getSize() > 0)  {
         lDr = logDr(iSpecies);
         cout << ", nCoreEl=" <<
            (
             (SxVector<double> (radCore(iSpecies)(1))
            * SxVector<double> (radCore(iSpecies)(0)).cub()).integrate(lDr)
           );
      }
      if (realSpace(iSpecies)) cout << ", rsp";
      cout << endl;
      cout << "|             Potential:   " << "file:/"
                                            << pseudoPotFiles(iSpecies)
                                            << endl;
   }
}

void SxPseudoPot::reConfinePsi (double sigma) 
{
   int iSpecies, il;
   double norm;


   for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)  {
      for (il = 0; il < pseudoPsi(iSpecies).getSize(); il++)  {
         
               SxVector<double> &psi     = pseudoPsi(iSpecies)(il);
         const SxVector<double> &radFunc = rad(iSpecies);

         // mutiply psi by Gaussian : Exp(-r^2/(2*sigma^2))
         psi *= exp ( (- radFunc.sqr()/(2.0*sigma*sigma)) );

         // --- re-normalize
         norm = sqrt((radFunc.cub() * psi.sqr()).integrate(logDr(iSpecies)));
         if (fabs (norm) > 1e-10) psi /= norm;

      }
   }
   
}

int SxPseudoPot::getLMax () const
{
   int lmax = -1, ll;
   for (int iSpecies = 0; iSpecies < lMax.getSize (); ++iSpecies)
      if (lmax < (ll = lMax(iSpecies))) lmax = ll;
   SX_CHECK (lmax >= 0);
   return lmax;
}

SxArray<SxArray<SxVector<double> > > SxPseudoPot::getPseudoPsi ()
{
   SxArray<SxArray<SxVector<double> > > result (nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      result(is).resize(lMax(is)+1);
      for (int l = 0; l <= lMax(is); l++)   {
         result(is)(l) = pseudoPsi(is)(l);
      }
   }
   return result;
}

const SxVector<double>& SxPseudoPot::getPseudoPsi (int is, int l) const
{
   return pseudoPsi(is)(l);
}

SxMeshG
SxPseudoPot::getAtomRhoG(const SxGBasis &G, int iSpecies, int iSpin) const
{
   SX_CHECK(G.structPtr);
   double volume = G.structPtr->cell.volume;
   SxRadBasis r (rad, logDr);  // |r>
   const SxVector<double> &psRad = rad(iSpecies);
   SxVector<double> radFun (psRad.getSize ());
   radFun.set(0.);
   radFun.setBasis (&r);
   for (int l=0; l <= lMax(iSpecies); l++)  {
      double focc = 0.;
      if (iSpin >= 0) {
         focc =  foccAtomicRho(iSpecies)(l)(iSpin);
      } else {
         SX_LOOP(jSpin) focc += foccAtomicRho(iSpecies)(l)(jSpin);
      }
      radFun += focc * pseudoPsi(iSpecies)(l).sqr();
   }
   //cout << "AO density should be projected Dirac-like" << endl;
   SxMeshG psiAtom  = r.toPWBasis  (psRad,radFun, G, 0, logDr(iSpecies));
   psiAtom *= SxYlm::getYlmNormFactor (0,0) * SQRT_1_4PI // = Y(0,0)
           * FOUR_PI / sqrt(volume); // normalization
   return psiAtom;
}
