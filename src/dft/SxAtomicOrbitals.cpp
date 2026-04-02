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

#include <SxAtomicOrbitals.h>
#include <SxNaturalCubicSpline.h>
#include <SxCubicSpline.h>
#include <SxRadialBasis.h>
#include <SxRegex.h>
#include <SxFileIO.h>
#include <SxTextIO.h>
#include <SxSimpleParser.h>

SxAtomicOrbitals::SxAtomicOrbitals ()
{
   // empty
}

SxAtomicOrbitals::SxAtomicOrbitals (
      const SxArray<SxArray<SxRadBasis::TPsi> > &in,
      const SxConstPtr<SxRadBasis> &radBasisPtrIn)
   : SxPsiSet (PW)
{
   SX_CHECK (radBasisPtrIn.getPtr () != NULL);

   radBasisPtr = radBasisPtrIn;
   muSet     =  in;

   // --- make sure that basis is set
   const SxRadBasis *radPtr = radBasisPtr.getPtr ();
   int iSpecies, nSpecies = (int)in.getSize();
   SX_CHECK (nSpecies <= radBasisPtrIn->radFunc.getSize(),
             nSpecies,   radBasisPtrIn->radFunc.getSize());
   for (iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      int nl = (int)in(iSpecies).getSize();
      for (int l = 0; l < nl; l++)  {
         // check dimensions
         SX_CHECK (in(iSpecies)(l).getSize() ==
                   radBasisPtrIn->radFunc(iSpecies).getSize(),
                   in(iSpecies)(l).getSize(),
                   radBasisPtrIn->radFunc(iSpecies).getSize());
         // set basis
         muSet(iSpecies)(l).setBasis (radPtr);
      }
   }

   createFuncLMap ();

   registerMemoryObservers ();
}

SxAtomicOrbitals::SxAtomicOrbitals (
      const SxArray<SxVector<double> > &in,
      const SxConstPtr<SxRadBasis> &radBasisPtrIn)
   : SxPsiSet (PW)
{
   // --- check dimensions
   int nSpecies = (int)in.getSize();
   muSet.resize(nSpecies);
   SX_CHECK (nSpecies == radBasisPtr->radFunc.getSize(),
             nSpecies, radBasisPtr->radFunc.getSize());
   for (int iSpecies=0; iSpecies < nSpecies; iSpecies++)  {
      int nl = (int)in(iSpecies).getNCols ();
      muSet(iSpecies).resize(nl);
      for (int l=0; l < nl; l++)  {
         /*SX_CHECK (in(iSpecies).colRef(l).getSize() ==
                     radBasisPtr->radFunc(iSpecies).getSize(),
                     in(iSpecies).colRef(l).getSize(),
                     radBasisPtr->.radFunc(iSpecies).getSize());*/
         muSet(iSpecies)(l) = in(iSpecies).colRef(l);
      }
   }
   radBasisPtr = radBasisPtrIn;

   createFuncLMap ();

   registerMemoryObservers ();
}


SxAtomicOrbitals::SxAtomicOrbitals (const SxAtomicOrbitals &in)
   : SxPsiSet(PW)
{
   (*this) = in;
}

SxAtomicOrbitals::SxAtomicOrbitals (SxBinIO &io)
{

   radBasisPtr = SxPtr<SxRadBasis>::create(io);
   try  {
      read(io);
   } catch (SxException e)   {
      e.print ();
      SX_EXIT;
   }
   createFuncLMap ();
}

// rescale radial coordinate of orbital
static
SxVector<double> scaleOrbital (double scale, const SxVecRef<double> &orbital)
{
   SX_CHECK(scale > 1.0);
   const SxVecRef<double> &rad = orbital.getBasis<SxRadBasis> ()
                                 .radFunc(orbital.auxData.is);
   int l = orbital.auxData.l;
   SxVector<double> scaledRad = scale * rad;
   scaledRad(0) = 0.0;
   SxCubicSpline spline(scaledRad, orbital,
                        (l < 2) ? SxCubicSpline::NaturalHermite
                                : SxCubicSpline::Hermite);
   SxVector<double> newOrbital = spline.getY (rad);
   if (l > 0) newOrbital(0) = 0.;
   newOrbital.auxData = orbital.auxData;
   return newOrbital;
}

void SxAtomicOrbitals::setup (const SxSymbolTable *table)
{

   SX_CHECK (table);
   // SxRadBasis to be set up here ...
   SxPtr<SxRadBasis> newRad = SxPtr<SxRadBasis>::create ();
   radBasisPtr = newRad; // ... and kept as const object when we finish here

   // fromPotential cache
   class {
      private:
         SxMap<SxString, SxPtr<SxAtomicOrbitals> > orbs;

         /// the actual read routine
         void read (const SxString &name) {
            SxParser::Table potTable = SxParser ().read (name);
            SxPtr<SxRadBasis> rad;
            if (potTable->containsGroup("pseudoPot"))   {
               SxPseudoPot pot (&*potTable);
               orbs(name) = SxPtr<SxAtomicOrbitals>::create (
                               pot.getPseudoPsi (),
                               SxPtr<SxRadBasis>::create(pot.rad, pot.logDr));
            } else if (potTable->containsGroup("pawPot"))  {
               SxPAWPot pot (&*potTable);
               orbs(name) = SxPtr<SxAtomicOrbitals>::create (pot.getPhiPS (),
                               SxPtr<SxRadBasis>::create(pot.rad, pot.logDr));
            } else {
               cout << "No known potential group found in file '"
                    << name << "'." << endl;
               cout << "SPHInX quits here!" << endl;
               SX_QUIT;
            }
         }
      public:
         /// Get rad basis from a sx input file
         const SxRadBasis& getBasis (const SxString &name)
         {
            if (!orbs.hasKey (name)) read (name);
            return orbs(name)->getRadBasis ();
         }
         /// Get orbitals from a sx input file
         const SxAtomicOrbitals& getOrbs (const SxString &name)
         {
            if (!orbs.hasKey (name)) read (name);
            return *orbs(name);
         }
   } fileCache;

   SYMBOLPARSE (table) {
      int iSpecies=-1;

      // get species
      FOREACH_SYMBOLGROUP("species")  {
         iSpecies++;
         // --- basis
         SYMBOLGROUP ("radBasis")  {
            SYMBOLGROUP("fromFile")  {
               SxString fileName = SYMBOLGET("file");
               int is = SYMBOLGET("is");
               SxRadBasis fileBasis;
               fileBasis.read(fileName);
               newRad->addMesh(fileBasis.radFunc(is), fileBasis.logDr(is));
            }
            SYMBOLGROUP("generate") {
               double rMin = SYMBOLGET("rMin");
               double rMax = SYMBOLGET("rMax");
               int nPoints = SYMBOLGET("nPoints");
               SxRadBasis genBasis;
               genBasis.set(rMin, rMax, nPoints);
               newRad->addMesh(genBasis.radFunc(0),genBasis.logDr(0));
            }
            SYMBOLGROUP("fromPotential")  {
               SxString potFile = SYMBOLGET("file");
               const SxRadBasis &fileBasis = fileCache.getBasis (potFile);
               int is = SYMBOLGET("is");
               newRad->addMesh (fileBasis.radFunc (is), fileBasis.logDr (is));
            }
            if (newRad->radFunc.getSize () != iSpecies + 1)  {
               // should never happen - invalid input file
               cout << "No known radial basis initialization specified!"
                    << endl;
               cout << "SPHInX quits here!" << endl;
               SX_QUIT;
            }
         }
         // --- orbitals
         int iOrbital = -1;
         SxVector<double> newOrbital;
         FOREACH_SYMBOLGROUP("orbital")
         {
            iOrbital++;
            SYMBOLGROUP("fromPotential") {
               SxString potFile = SYMBOLGET("file");
               const SxAtomicOrbitals &fromPot = fileCache.getOrbs (potFile);
               int l   = SYMBOLGET ("l");
               int iot = SYMBOLGET("iot");
               int is  = SYMBOLGET("is");
               newOrbital = fromPot(is,iot);
               // Check Basis
               if (fromPot.getBasis().getPtr() != radBasisPtr.getPtr())  {
                  newOrbital
                     = fromPot.getBasis()
                     ->changeRadBasis(radBasisPtr.getPtr (),newOrbital);
               }
               if (l != newOrbital.auxData.l) {
                  cout << "Inconsistent iot/l combination in fromPotential!"
                       << endl;
                  cout << "l is " << l
                       <<", should be " << int(newOrbital.auxData.l)
                       << endl;
                  SX_QUIT;
               }
               newOrbital.auxData.is = iSpecies;
               if (HAVE_SYMBOL("scale"))
                  newOrbital = scaleOrbital (SYMBOLGET("scale"), newOrbital);
               cout << "Add Orbital from Potential: is iot l " << is
                    << " " << iot << " " << l;
               addOrbital(newOrbital);
            }
            SYMBOLGROUP("fromFile")  {
               SxString fileName = SYMBOLGET("file");
               int l = SYMBOLGET("l");
               int iot = SYMBOLGET("iot");
               int is = SYMBOLGET("is");
               // read Orbitalfile & check Basis
               SxAtomicOrbitals fileOrbitals;
               if (HAVE_SYMBOL("Siesta"))  {
                  fileOrbitals.readSiesta(fileName,radBasisPtr->radFunc(is));
                  newOrbital = fileOrbitals(is,iot);
                  newOrbital.setBasis(&*radBasisPtr);
               } else  {
                  SxRadBasis vecBasis;
                  vecBasis.read(fileName);
                  fileOrbitals.read(fileName);
                  int nSpeciesFile = fileOrbitals.getNSpecies ();
                  if (is >= nSpeciesFile) {
                     cout << "Did not find species "
                          << is << " in quamol file " << fileName
                          << " with " << nSpeciesFile << " species." << endl;
                     SX_QUIT;
                  }
                  int nOrbTypesFile = fileOrbitals.getNOrbTypes(is);
                  if (iot >= nOrbTypesFile) {
                     cout << "Did not find orbitaltype iot="
                          << iot << " in quamol file " << fileName
                          << " with " << nOrbTypesFile << " orbital types."
                          << endl;
                     SX_QUIT;
                  }
                  if (HAVE_SYMBOL("rCut"))  {
                     double rCut = SYMBOLGET("rCut");
                     double rMin = vecBasis.getRMin(is);
                     double rMax = vecBasis.getRMax(is);
                     int nPoints = int ((double)vecBasis.getNPoints (is)
                                        * rCut / rMax);
                     SxRadialBasis cutRad (rMin, rCut, nPoints,
                                           /*realSpace=*/true);
                     newOrbital = vecBasis | ( cutRad | fileOrbitals(is,iot));
                  } else {
                     newOrbital = fileOrbitals(is,iot);
                     newOrbital.setBasis(vecBasis);
                  }

                  if (fabs (vecBasis.radFunc(is)(0)
                           - radBasisPtr->radFunc(iSpecies)(0)) > 1e-6 ||
                      fabs (vecBasis.logDr(is)
                            - radBasisPtr->logDr(iSpecies)) > 1e-6)
                     newOrbital = vecBasis.changeRadBasis(
                           radBasisPtr.getPtr (),
                           newOrbital);
               }
               if (l != newOrbital.auxData.l) {
                  cout << "Inconsistent iot/l combination in fromFile!"
                       << endl;
                  cout << "l is " << l
                       <<", should be " << int(newOrbital.auxData.l)
                       << endl;
                  SX_QUIT;
               }
               newOrbital.auxData.is = iSpecies;
               if (HAVE_SYMBOL("scale"))
                  newOrbital = scaleOrbital (SYMBOLGET("scale"), newOrbital);
               cout << "Add Orbital from File: is iot l " << is
                    << " " << iot
                    << " " << l;
               addOrbital(newOrbital);
            }
            SYMBOLGROUP("fromWaves")  {
               SxString fileName = SYMBOLGET("file");
               int l      = SYMBOLGET("l");
               int iState = SYMBOLGET("iState");
               newOrbital = compressWave(fileName,iState,iSpecies,l);
               newOrbital.auxData.is = iSpecies;
               newOrbital.auxData.ia = -1;
               newOrbital.auxData.n = -1;
               newOrbital.auxData.l = char(l);
               newOrbital.auxData.m = NONE_M;
               cout << "Add Orbital from wavestate: iState l "
                    << iState << " " << l;
               addOrbital(newOrbital);
            }
            SYMBOLGROUP("generateGaussian")  {
               double beta = SYMBOLGET("beta");
               int l = SYMBOLGET("l");
               if (l > 0) beta = beta * sqrt(2./l);
               const SxVector<double> &r = radBasisPtr->radFunc(iSpecies);
               SxVector<double> rPowL = r;
               if (l == 0) rPowL.set(1.0);
               for (int il = 2; il <= l; il++) rPowL *= r;
               newOrbital = rPowL * exp(-beta * r.sqr ());
               newOrbital.auxData.is = iSpecies;
               newOrbital.auxData.ia = -1;
               newOrbital.auxData.n = -1;
               newOrbital.auxData.l = char(l);
               newOrbital.auxData.m = NONE_M;
               newOrbital.setBasis(radBasisPtr.getConstPtr());
               if (l == 0)  {
                  cout << "Generate Gaussian Orbital: is l beta ";
                  cout << iSpecies;
                  cout << " " << l;
                  cout << " " << beta;
               } else  {
                  cout << "Generate Gaussian Orbital: is l rMax ";
                  cout << iSpecies;
                  cout << " " << l;
                  cout << " " << beta / sqrt(2./l);
               }
               addOrbital(newOrbital);
            }
            SYMBOLGROUP("gaussianBasis")  {
               SxVector<double> exps (SYMBOLGET("exponents")->toList());
               int nGauss = (int)exps.getSize ();
               SxVector<double> coeffs (SYMBOLGET("coefficients")->toList());
               int l = SYMBOLGET("l");
               if (nGauss != coeffs.getSize())  {
                  cout << "Wrong Number of coefficients: nGaussians is "
                       << nGauss << ", but have " << coeffs.getSize ()
                       << "coefficients!" << endl;
                  SX_QUIT;
               }
               const SxVector<double> &r = radBasisPtr->radFunc(iSpecies);
               int dim = (int)r.getSize();
               newOrbital.resize (dim);
               newOrbital.set(0.0);
               for (int iGauss = 0; iGauss < nGauss; iGauss++)  {
                  newOrbital
                     += coeffs(iGauss) * exp(-r.sqr () * exps(iGauss));
               }
               newOrbital *= pow(r, double(l));
               newOrbital.auxData.is = iSpecies;
               newOrbital.auxData.ia = -1;
               newOrbital.auxData.n = -1;
               newOrbital.auxData.l = char(l);
               newOrbital.auxData.m = NONE_M;
               newOrbital.setBasis(radBasisPtr.getConstPtr());
               cout << "Generate orbital via gaussian basis: is l nGauss ";
               cout << iSpecies;
               cout << " " << l;
               cout << " " << nGauss;
               addOrbital(newOrbital);
            }
            if (getNOrbTypes (iSpecies) != iOrbital + 1) {
               // should never happen - invalid input file
               cout << "No known Initialization shema found!" << endl;
               cout << "SPHInX quits here!" << endl;
               SX_QUIT;
            }
         }
      }
   }

   createFuncLMap ();
}

void SxAtomicOrbitals::addOrbital (const SxVecRef<double> &orbitalIN)
{
   SX_CHECK(radBasisPtr.getPtr () != NULL);
   
   int iSpecies = orbitalIN.auxData.is;
   int lIN =  orbitalIN.auxData.l;
   int iot = 0;
   if (iSpecies < getNSpecies ()) {
      int nOrbTypes = getNOrbTypes(iSpecies);
      muSet(iSpecies).resize(nOrbTypes + 1,true);
      for (iot = nOrbTypes; iot > 0; iot++)  {
         if (muSet(iSpecies)(iot-1).auxData.l > lIN)  {
            // shift existing orbitals with l > lIN
            muSet(iSpecies)(iot) = std::move(muSet(iSpecies)(iot-1));
            muSet(iSpecies)(iot).auxData.n = char(iot);
         } else {
            break;
         }
      }
#ifndef NDEBUG
      // make sure that orbitals below iot have smaller l's
      for (int jot = 0; jot < iot; jot++)  {
         SX_CHECK (muSet(iSpecies)(jot).auxData.l <= lIN,
                   iSpecies, jot, (int)muSet(iSpecies)(jot).auxData.l, lIN);
      }
#endif
      // shift iot's in funcLMap for all l's larger than current one
      int lmax = (int)funcLMap(iSpecies).getSize () - 1;
      for (int l = lIN + 1; l <= lmax; l++)
         SX_LOOP(ifl) funcLMap(iSpecies)(l)(ifl)++;
      // append new orbital in funcLMap
      if (lIN <= lmax)  {
         // append for existing l
         int nOrbL = (int)funcLMap(iSpecies)(lIN).getSize ();
         funcLMap(iSpecies)(lIN).resize (nOrbL + 1);
         funcLMap(iSpecies)(lIN)(nOrbL) = iot;
      } else {
         // add new l
         funcLMap(iSpecies).resize (lIN + 1, true);
         funcLMap(iSpecies)(lIN).resize (1);
         funcLMap(iSpecies)(lIN)(0) = iot;
      }
   } else {
      // append new species
      muSet.resize (iSpecies + 1, true);
      muSet(iSpecies).resize(1);
      funcLMap.resize (iSpecies + 1, true);
      funcLMap(iSpecies).resize (lIN + 1);
      funcLMap(iSpecies)(lIN).resize (1);
      funcLMap(iSpecies)(lIN)(0) = iot;// = 0
   }
   // insert orbital now
   muSet(iSpecies)(iot) = orbitalIN;
   muSet(iSpecies)(iot).auxData.is = iSpecies;
   muSet(iSpecies)(iot).auxData.ia = -1;
   muSet(iSpecies)(iot).auxData.n = char(iot);
   muSet(iSpecies)(iot).auxData.l = char(lIN);
   muSet(iSpecies)(iot).auxData.m = NONE_M;
   muSet(iSpecies)(iot).setBasis(radBasisPtr.getConstPtr ());
   cout << " as " << iSpecies << " " << iot << " " << lIN << endl;
}

SxAtomicOrbitals::~SxAtomicOrbitals ()
{
   // empty
}

void SxAtomicOrbitals::operator= (const SxAtomicOrbitals &in)
{
   if (&in == this) {
      cout << "a = a Error" << endl;
      SX_EXIT;
   }

   radBasisPtr = in.radBasisPtr;
   muSet = in.muSet;
   funcLMap = in.funcLMap; 
}

void SxAtomicOrbitals::operator+= (const SxAtomicOrbitals &in)
{
   SX_CHECK(radBasisPtr.getPtr() != NULL);
   SX_CHECK(radBasisPtr == in.getRadBasisPtr());
   SX_CHECK(muSet.getSize() == in.muSet.getSize(),
            muSet.getSize(),
            in.muSet.getSize());
   for (int is = 0; is < in.muSet.getSize (); is++)   {
      SX_CHECK(muSet(is).getSize() == in.muSet(is).getSize(),
               muSet(is).getSize(),
               in.muSet(is).getSize());
      for (int l = 0; l < in.muSet(is).getSize (); l++)   {
         muSet(is)(l) += in.muSet(is)(l);
      }
   }
}

SxAtomicOrbitals SxAtomicOrbitals::operator* (double skalar) const
{
   SX_CHECK(radBasisPtr.getPtr() != NULL);
   int nSpecies = (int)muSet.getSize();
   SxArray<SxArray<SxRadBasis::TPsi> > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes(is);
      resMuSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         resMuSet(is)(iot) = skalar * muSet(is)(iot);
      }
   }
   SxAtomicOrbitals result = SxAtomicOrbitals(resMuSet, radBasisPtr);

   return result;
}

SxRadBasis::TPsi & SxAtomicOrbitals::operator() (int is, int iot)
{
   return muSet(is)(iot);
}

const SxRadBasis::TPsi & SxAtomicOrbitals::operator() (int is, int iot) const
{
   return muSet(is)(iot);
}

int SxAtomicOrbitals::getIOT (int is, int n, int l) const
{
   int iot = 0;
   while ((muSet(is)(iot).auxData.n != n)
         || (muSet(is)(iot).auxData.l != l))   {
      iot++;
      if (iot >= muSet(is).getSize()) SX_EXIT;
   }

   return iot;
}


SxRadBasis::TPsi SxAtomicOrbitals::operator() (int is,
                                               int ia, 
                                               int n, 
                                               int l,
                                               int m) const
{
   int iot = getIOT (is,n,l);

   SxRadBasis::TPsi vec ( muSet(is)(iot) );
   vec.setBasis (radBasisPtr.getConstPtr ());

   vec.auxData.is = is; vec.auxData.ia = ia;
   vec.auxData.n  = char(n);  vec.auxData.l  = char(l);
   vec.auxData.m  = char(m);

   return vec;
}

int SxAtomicOrbitals::getNSpecies () const
{
   return (int)muSet.getSize ();
}

int SxAtomicOrbitals::getNOrbTypes (int iSpecies) const
{
   return (int)muSet(iSpecies).getSize ();
}

int SxAtomicOrbitals::getNOrbTypes () const
{
   int result = 0;
   for (int is = 0; is < getNSpecies (); is++)
      result += getNOrbTypes(is);

   return result;
}

void SxAtomicOrbitals::set(double val)
{

   for (int iSpecies = 0; iSpecies < muSet.getSize(); iSpecies++)   {
      for (int iot = 0; iot < muSet(iSpecies).getSize(); iot++)   {
         muSet(iSpecies)(iot).set(val);
      }
   }
}

SxQuantumNumbers SxAtomicOrbitals::getQuantumNumbers (int iSpecies, int idx) const
{
   SX_CHECK(iSpecies < muSet.getSize(), iSpecies, muSet.getSize());
   SX_CHECK(idx < muSet(iSpecies).getSize(), idx, muSet(iSpecies).getSize());
   int n = muSet(iSpecies)(idx).auxData.n;
   int l = muSet(iSpecies)(idx).auxData.l;

   SxQuantumNumbers result (iSpecies,n,l,0);

   return result;
}

SxArray<SxQuantumNumbers> SxAtomicOrbitals::getReducedOrbitalMap () const
{
   SX_CHECK(&muSet);
   int nSpecies = (int)muSet.getSize ();
   SxList<SxQuantumNumbers> list;
   for(int is = 0; is < nSpecies; is++)   {
      int nOrbLocal = (int)muSet(is).getSize();
      for(int iot = 0; iot < nOrbLocal; iot++)   {
         int l = muSet(is)(iot).auxData.l;
         // SxQuantumNumberConstructor forbids uninitialized ia and m
         SxQuantumNumbers qNumbers (is,0,iot,l,0);
         // Therefore unset it here
         qNumbers.iAtom = -1; qNumbers.m = NONE_M; 
         list.append(qNumbers);
      }
   }
   SxArray<SxQuantumNumbers> result = SxArray<SxQuantumNumbers>(list);

   return result;
}


SxArray<SxQuantumNumbers>
SxAtomicOrbitals::getOrbitalMap (const SxAtomicStructure &structure) const
{
   SX_CHECK(&muSet);
   int nSpecies = (int)muSet.getSize ();
   SxList<SxQuantumNumbers> list;
   for(int is = 0; is < nSpecies; is++)   {
      int nAtoms = structure.getNAtoms(is);
      int nOrbLocal = (int)muSet(is).getSize();
      for(int ia = 0; ia < nAtoms; ia++)   {
         for(int iot = 0; iot < nOrbLocal; iot++)   {
            int l = muSet(is)(iot).auxData.l;
            for(int m = -l; m <= l; m++)   {
               list.append(SxQuantumNumbers(is,ia,iot,l,m));
            }
         }
      }
   }

   SxArray<SxQuantumNumbers> result = SxArray<SxQuantumNumbers>(list);

   return result;
}

int SxAtomicOrbitals::getNOrbitals (const SxAtomicStructure &structure) const
{
   SX_CHECK(muSet.getSize() > 0);
   int nSpecies = (int)muSet.getSize ();
   int result = 0;
   for(int is = 0; is < nSpecies; is++)   {
      int nAtoms = structure.getNAtoms(is);
      int nOrbLocal = (int)muSet(is).getSize();
      for(int ia = 0; ia < nAtoms; ia++)   {
         for(int iot = 0; iot < nOrbLocal; iot++)   {
            int l = muSet(is)(iot).auxData.l;
            result += 2*l+1;
         }
      }
   }

   return result;
}

int SxAtomicOrbitals::getOrbitalIdx (int is, int ia, int iot, int l, int m, const SxArray<SxQuantumNumbers> &map) const
{
   for (int iOrbital = 0; iOrbital < map.getSize(); iOrbital++)  {
      if (map(iOrbital).iSpecies == is  &&
          map(iOrbital).iAtom    == ia  &&  
          map(iOrbital).n        == iot &&
          map(iOrbital).l        == l   &&
          map(iOrbital).m        == m) return iOrbital;
   }

   cout << "OrbitalIdx with is = " << is << ", ia = " << ia
        << ", n = " << iot << ", l = " << l << ", m = " << m
        << " not found!" << endl;
   SX_EXIT;
}


void SxAtomicOrbitals::createFuncLMap ()
{
   int nSpecies = (int)muSet.getSize ();
   funcLMap.resize(nSpecies);
   for(int iSpecies = 0; iSpecies < nSpecies; iSpecies++)   {
      int lMax = getLMax (iSpecies);
      funcLMap(iSpecies).resize(lMax + 1);
      SxArray<int> funcPerL(lMax + 1);
      funcPerL.set(0);
      int nOrbTypes = (int)muSet(iSpecies).getSize ();
      for(int iot = 0; iot < nOrbTypes; iot++)   {
         int l = muSet(iSpecies)(iot).auxData.l;
         funcPerL(l)++;
      }
      SxArray<int> currentIFL(lMax + 1);
      currentIFL.set(0);
      for (int l = 0; l <= lMax; l++)  {
         funcLMap(iSpecies)(l).resize(funcPerL(l));
      }
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         int l = muSet(iSpecies)(iot).auxData.l;
         funcLMap(iSpecies)(l)(currentIFL(l)) = iot;
         currentIFL(l)++;
      }
   }
}

SxRadBasis::TPsi SxAtomicOrbitals::getFuncL (int is, int l, int ifl)
{
   int iot = funcLMap(is)(l)(ifl);
   return muSet(is)(iot);
}

const SxRadBasis::TPsi SxAtomicOrbitals::getFuncL (int is, int l, int ifl) const
{
   int iot = funcLMap(is)(l)(ifl);
   return muSet(is)(iot);
}


int SxAtomicOrbitals::getLMax () const
{
   int nSpecies = (int)muSet.getSize ();
   int result = 0;
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = (int)muSet(is).getSize ();
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         if (result <  muSet(is)(iot).auxData.l)
            result = muSet(is)(iot).auxData.l;
      }
   }
   return result;
}

int SxAtomicOrbitals::getLMax (const int iSpecies) const
{
   int nOrbTypes = (int)muSet(iSpecies).getSize ();
   int result = 0;
   for (int iot = 0; iot < nOrbTypes; iot++)   {
         if (result <  muSet(iSpecies)(iot).auxData.l)
            result = muSet(iSpecies)(iot).auxData.l;
   }
   return result;
}

SxArray<SxArray<int> > SxAtomicOrbitals::getFuncPerL () const
{
   int nSpecies = (int)muSet.getSize ();
   SxArray<SxArray<int> > result(nSpecies);
   for(int iSpecies = 0; iSpecies < nSpecies; iSpecies++)   {
      int lMax = getLMax (iSpecies);
      result(iSpecies).resize(lMax + 1);
      for(int l = 0; l <= lMax; l++)   {
         result(iSpecies)(l) = getFuncPerL(iSpecies,l);
      }
   }
   return result;   
}

int SxAtomicOrbitals::getFuncPerL (int is, int l) const
{
   return (int)funcLMap(is)(l).getSize();
}

void SxAtomicOrbitals::registerMemoryObservers ()
{
   TRACK_MEMORY (muSet);
}

void SxAtomicOrbitals::write (SxBinIO &io) const
{
   SX_CHECK(radBasisPtr.getPtr () != NULL);
   const SxRadBasis &radBasis = *radBasisPtr;
   int nSpecies = getNSpecies();

   try {
      //create dimensions
      io.addDimension ("nSpecies", nSpecies);
      for (int is = 0; is < nSpecies; is++)   {
         io.addDimension ("dimOrbTypes-"+SxString(is), getNOrbTypes(is));
         io.addDimension ("dimRad-"+SxString(is),
                          (int)radBasis.radFunc(is).getSize());
      }
      //write data
      SxArray<SxVector<int> > lNumbers (nSpecies);
      for (int is = 0; is < nSpecies; is++)   {
         SxString radialName= "radFunc-" + SxString(is);
         SxString dimRadName = "dimRad-" + SxString(is);
         io.write (radialName, radBasis.radFunc(is), dimRadName);
         SxString logDrName= "logDr-" + SxString(is);
         io.write (logDrName, radBasis.logDr(is));
         lNumbers(is).resize(getNOrbTypes(is));
         for (int iot=0; iot < getNOrbTypes(is); iot++)   {
            radialName= "radial-" + SxString(is) + "-" + SxString(iot);
            io.write (radialName, muSet(is)(iot), dimRadName);
            lNumbers(is)(iot) = muSet(is)(iot).auxData.l;
         }
         SxString lName= "lNumbers-" + SxString(is);
         io.write (lName, lNumbers(is),"dimOrbTypes-"+SxString(is));
      }
   } catch (SxException e)   {
      e.print ();
      SX_EXIT;
   }
}

void SxAtomicOrbitals::write (const SxString &filename) const
{
   SxBinIO io(filename, SxBinIO::BINARY_WRITE_ONLY);
   write(io);
   io.setMode (SxBinIO::WRITE_DATA);
   write(io);
   io.close();
}

void SxAtomicOrbitals::read (SxBinIO &io)
{

   try {
      SxPtr<SxRadBasis> basisPtr = SxPtr<SxRadBasis>::create(io);
      radBasisPtr = basisPtr;
      //get dimensions
      int nSpecies = io.getDimension ("nSpecies");
      muSet.resize(nSpecies);
      for (int is = 0; is < nSpecies; is++)   {
         int nOrbTypes = io.getDimension ("dimOrbTypes-"+SxString(is));
         muSet(is).resize(nOrbTypes);
         int dimRad = io.getDimension ("dimRad-"+SxString(is));
         SxVector<int> lVec (nOrbTypes);
         SxString lName= "lNumbers-" + SxString(is);
         io.read (lName,&lVec,nOrbTypes);
         // TODO further check of SxRadBasis
         for (int iot=0; iot < nOrbTypes; iot++)   {
            muSet(is)(iot).resize(dimRad);
            SxString radialName= "radial-" + SxString(is) + "-" + SxString(iot);
            SxVector<double> &readVec = muSet(is)(iot);
            io.read (radialName,&readVec,dimRad);
            muSet(is)(iot).setBasis (radBasisPtr.getConstPtr());
            muSet(is)(iot).auxData.is = is;
            muSet(is)(iot).auxData.n = char(iot);
            muSet(is)(iot).auxData.l = char(lVec(iot));
            muSet(is)(iot).auxData.m = 0;
         }
         // TODO Inconsistent check --- is in file and is in RadBasis might differ
         /*
         if (dimRad != rBasis.radFunc(is).getSize())   {
            cout << "WARNING: Incompatible grid size in readAtomicOrbitals !" << endl;
            cout << "Grid size in File is " << dimRad << endl;
            cout << "Grid size in rad Basis is " 
                 << rBasis.radFunc(is).getSize() << endl;
            cout << "Interpolate now!" << endl;
            SxRadBasis fileBasis;
            fileBasis.read(io);
            for (int iot=0; iot < nOrbTypes; iot++)   { 
               SxNaturalCubicSpline interpol (toVector(fileBasis.radFunc(is)), 
                     toVector(muSet(is)(iot)));
               muSet(is)(iot).resize(rBasis.radFunc(is).getSize ());
               for (int i = 0; i < rBasis.radFunc(is).getSize (); i++)   {
                  int last = fileBasis.radFunc(is).getSize() - 1;
                  if (fileBasis.radFunc(is)(last) > rBasis.radFunc(is)(i))
                     muSet(is)(iot)(i) = interpol.getVal(rBasis.radFunc(is)(i));
                  else muSet(is)(iot)(i) = 0.0;
               }
            }
         }
         */
      }
   } catch (SxException e)   {
      e.print ();
      SX_EXIT;
   }

   createFuncLMap ();
}

void SxAtomicOrbitals::read (const SxString &file)
{
   try  {
      SxBinIO io (file, SxBinIO::BINARY_READ_ONLY);
      read (io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}
void SxAtomicOrbitals::writeSiesta (const SxArray<SxString> &file, 
                                    const SxArray<SxVector<double> > &basis)
{
   int nSpecies = getNSpecies ();
   SX_CHECK (basis.getSize () == nSpecies);
   SX_CHECK (file.getSize () == nSpecies);

   for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)   {
      SxTextIO io(file(iSpecies));
      int oldL = -1;
      int z = 0;
      for (int iot = 0; iot < getNOrbTypes (iSpecies); iot++)  {
         // orbital header
         int l = muSet(iSpecies)(iot).auxData.l;
         if (l == oldL)  z++;
         else {
            z = 0;
            oldL = l;
         }
         io.printf("%i ? %i ? 0.000000 #orbital l, n, z, is_polarized, population\n",
                   l,z);
         int npts = (int)basis(iSpecies).getSize ();
         double delta = basis(iSpecies)(1) - basis(iSpecies)(0);
         double cutoff = basis(iSpecies)(npts-1);
         io.printf("%i %.12e %.12f # npts, delta, cutoff\n",npts, delta, cutoff);
         SxCubicSpline spline (
               radBasisPtr->radFunc(iSpecies), 
               muSet(iSpecies)(iot), 
               SxCubicSpline::Natural);
         SxVector<double> orbital = spline.getY(basis(iSpecies));
         orbital(npts-1) = 0.0;
         if (l != 0) orbital(0) = 0.0;
         // data
         io.writeXYPlot (basis(iSpecies), orbital);
      }
   }
}


void SxAtomicOrbitals::readSiesta (const SxString &file, 
                                   const SxVecRef<double> &radFunc)
{
   SxString ionFile;
   SxList<SxString> cLine;
   SxString data;
   muSet.resize(1);

   try {
      ionFile = SxFileIO::readBinary (file,-1);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   SxList<SxString> ionFileTok = ionFile.tokenize ('\n');
   int line = 0;
   while (ionFileTok(line).contains("# Lmax for basis, no. of nl orbitals") == 0) 
      if (line >= ionFileTok.getSize()) {
         cout << "Expression not found" << endl;
         SX_EXIT;
      }
   data = ionFileTok(line).left("#").stripWhiteSpace ();
   //int lMax = (data.left(" ").stripWhiteSpace ()).toInt ();
   int nOrbTypes = (data.right(" ").stripWhiteSpace ()).toInt ();
   muSet(0).resize(nOrbTypes);
   for (int iot = 0; iot < nOrbTypes; iot++)  {
      while (ionFileTok(line).contains("#orbital") == 0) {
         line++;
         if (line >= ionFileTok.getSize()) {
            cout << "Expression not found" << endl;
            SX_EXIT;
         }
      }
      cLine = ionFileTok(line).left('#').stripWhiteSpace ().tokenize(' ');
      int l = cLine.first ().toInt ();
      line++;
      cLine = ionFileTok(line).left('#').stripWhiteSpace ().tokenize(' ');
      int nPoints = cLine.first ().toInt ();
      //int cutoff = cLine.last ().toInt ();
      SxVector<double> x (nPoints);
      SxVector<double> y (nPoints);
      line++;
      for (int iPoint = 0; iPoint < nPoints; iPoint++,line++)  {
         x(iPoint) = ionFileTok(line).stripWhiteSpace ().left (' ')
            .stripWhiteSpace ().toDouble();
         y(iPoint) = ionFileTok(line).stripWhiteSpace ().right(' ')
            .stripWhiteSpace ().toDouble();
      }
      SxCubicSpline spline(x,y,
            SxCubicSpline::Natural);
      muSet(0)(iot) = spline.getY(radFunc);
      // Set auxdata
      SxRadBasis::TPsi &mu = muSet(0)(iot);
      mu.auxData.is = -1;
      mu.auxData.ia = -1;
      mu.auxData.n  = char(iot);
      mu.auxData.l  = char(l);
      mu.auxData.m  = NONE_M;
   }

   createFuncLMap ();

   registerMemoryObservers ();
}

void SxAtomicOrbitals::print (const SxString &file) const
{
   int nSpecies = getNSpecies();
   const SxRadBasis &rad = *radBasisPtr;
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes(is);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         SxString outFile = file + SxString(is) + SxString(iot) + ".dat";
         SxTextIO (outFile).writeXYPlot (rad.radFunc(is), muSet(is)(iot));
      }
   }
}

void SxAtomicOrbitals::print () const
{
   SxString file = "AtomicOrbitals";
   print (file);
}

void SxAtomicOrbitals::setBasis (SxConstPtr<SxRadBasis> radBasisPtrIn)
{
   radBasisPtr = radBasisPtrIn;
   for (int is = 0; is < muSet.getSize (); ++is)  {
      for (int iot = 0; iot < muSet(is).getSize (); ++iot)  {
         muSet(is)(iot).setBasis (radBasisPtr.getConstPtr ());
      }
   }
}

void SxAtomicOrbitals::readOrbital (FILE *fp, int is, int n, int l, int iot,
                                    bool ignoreComments)
{
   SX_CHECK (fp);
   SX_CHECK (is >= 0 && (   (is <= getNSpecies () && iot == -1)
             || (is < getNSpecies ())),
             is, getNSpecies (), iot);
   SX_CHECK (l >= 0, l);
   SX_CHECK ((iot >= 0 && iot < muSet(is).getSize ()) || (iot == -1),
             iot, muSet(is).getSize ());
   
   // --- read the file
   int N = 0;
   char c;
   SxStack<double> rFile, psiFile;
   while (!feof(fp))  {
      if (ignoreComments)  {
         while (fscanf (fp, " #%c", &c) == 1)  {
            // comment line
            while (c != '\n')  {
               if (feof(fp)) break;
               c = (char)fgetc (fp);
            }
         }
      }
      double r, f;
      if (fscanf(fp, "%lf %lf", &r, &f) == 2)  {
         rFile << r;
         psiFile << f;
         N++;
      } else {
         break;
      }
   }
   if (N == 0)  {
      cout << "Could not read orbital from text file!" << endl;
      SX_EXIT;
   }

   // --- transform data to vectors
   SxVector<double> r(rFile), psi(psiFile);

   // --- change radial grid?
   if (radBasisPtr)  {
      if (N != radBasisPtr->radFunc(is).getSize ()
          || (radBasisPtr->radFunc(is) - r).normSqr () > 1e-10)
      {
         // --- spline interpolation
         SxNaturalCubicSpline spline;
         if (l == 0)
            spline = SxNaturalCubicSpline (r, r * r * psi);
         else
            spline = SxNaturalCubicSpline (r, psi);
         // linear extrapolation below r0
         double rmax = r(r.getSize () - 1), rmin = r(0);
         double dpsi0 = (psi(1) - psi(0)) / (r(1) - r(0)),
                psi0 = psi(0);

         // --- get psi on rBasis radial grid
         psi = SxVector<double>((*radBasisPtr)(is));
         for (int ir = 0; ir < radBasisPtr->radFunc(is).getSize (); ++ir)  {
            double rr = radBasisPtr->radFunc(is)(ir);
            if (rr > rmax)  {
               psi(ir) = 0.;
            } else if (rr < rmin)  {
               psi(ir) = psi0 + (rr - rmin) * dpsi0;
            } else {
               psi(ir) = spline.getVal (rr);
               if (l == 0) psi(ir) /= rr * rr;
            }
         }
         SX_VALIDATE_VECTOR (psi);
      }
   }

   // --- metadata
   if (iot != -1 && muSet(is)(iot).getSize () > 0)
      psi.auxData = muSet(is)(iot).auxData;
   psi.auxData.is = is;
   psi.auxData.n = char(n);
   psi.auxData.l = char(l);

   // --- store away data in right place
   if (iot == -1)
      addOrbital (psi);
   else
      muSet(is)(iot) = psi;
}

void SxAtomicOrbitals::readOrbitals (FILE *fp, int is)
{
   char buffer[10240];
   int n, l;
   SxRegex lMatch("l *= *([0-9]+)");
   SxRegex nMatch("n *= *([0-9]+)");
   for ( ; !feof(fp); )  {
      buffer[0] = 0;
      char *res = fgets (buffer, 10240, fp);
      if (!res) {
         if (feof(fp)) break;
         cout << "Error while reading file" << endl;
         SX_QUIT;
      }
      SxString line(buffer);
      line = line.stripWhiteSpace ();
      if (line.getSize () <= 1) continue;
      SxList<SxString> match;
      n = muSet.getSize () >= is ? 0 : (int)muSet(is).getSize ();
      try {
         match = lMatch.match (line);
         if (match.getSize () == 2)  {
            cout << line;
            cout << match(0) << endl;
            l = match(1).toInt ();
         } else {
            cout << match << endl;
            cout << "Missing l in '" << line << "'." << endl;
            SX_EXIT;
         }
         match = nMatch.match (line);
         if (match.getSize () == 2) n = match(1).toInt ();
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
      readOrbital (fp, is, n, l, -1, false);
   }
}

double SxAtomicOrbitals::getNormSqr(int is, int iot) const
{
   const SxVector<double> &vec = muSet(is)(iot);
   double logDr = radBasisPtr->logDr(is);
   const SxVector<double> &r = radBasisPtr->radFunc(is);
   double result = (vec * vec * r * r * r).integrate(logDr);
   return result;
}

double SxAtomicOrbitals::getNormSqrSum() const
{
   double result = 0.0;
   int nSpecies = (int)muSet.getSize();
   for (int is = 0; is < nSpecies; is++)  {
      int nOrbTypes = (int)muSet(is).getSize ();
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         result += getNormSqr(is,iot);
      }
   }

   return result;
}

double SxAtomicOrbitals::dot(const SxAtomicOrbitals &in) const
{
   double result = 0.0;
   int nSpecies = getNSpecies ();
   SX_CHECK(nSpecies == in.getNSpecies());

   for (int is = 0; is < nSpecies; is++)  {
      SX_CHECK (getNOrbTypes(is) == in.getNOrbTypes (is));
      int nOrbTypes = (int)muSet(is).getSize ();
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         result += dot(in,is,iot,is,iot);
      }
   }

   return result;
}

double SxAtomicOrbitals::dot(const SxAtomicOrbitals &in, int is, int iot, int js, int jot) const
{
   const SxVector<double> &vec1 = muSet(is)(iot);
   const SxVector<double> &vec2 = in.muSet(js)(jot);

   double logDr = radBasisPtr->logDr(is);
   const SxVector<double> &r = radBasisPtr->radFunc(is);
   SX_CHECK (logDr == in.radBasisPtr->logDr(is));
   SX_CHECK ((r(0) - in.radBasisPtr->radFunc(is)(0)) < 1e-6);
   double result = (vec1 * vec2 * r * r * r).integrate(logDr);
   return result;
}

void SxAtomicOrbitals::normalize ()
{
   int nSpecies = (int)muSet.getSize();
   for (int is = 0; is < nSpecies; is++)  {
      int nOrbTypes = (int)muSet(is).getSize ();
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         double norm = getNormSqr(is,iot);
         if (norm > 1e-6) muSet(is)(iot) /= sqrt(norm);
         else {
            cout << "WARNING: Try to normalize zero norm orbital." << endl;
            cout << "Norm is " << norm << endl;
            muSet(is)(iot).set(0.0);
         }
      }
   }
}

SxArray<SxArray<SxVector<double> > > SxAtomicOrbitals::orthogonalize ()
{
   SX_CHECK(muSet.getSize() == funcLMap.getSize(),
            muSet.getSize(), funcLMap.getSize());
   // Gram-Schmidt via Cholesky
   SxAtomicOrbitals org = 1.0 * *this;
   int nSpecies = (int)muSet.getSize ();
   SxArray<SxArray<SxVector<double> > > result(nSpecies);
   for (int is = 0; is < nSpecies; is++)  {
      int lMax = getLMax(is);
      result(is).resize(lMax + 1);
      for (int l = 0; l <= lMax; l++)  {
         int nFL = getFuncPerL (is,l);
         if (nFL == 0) continue;
         SxVector<double> S = getOverlap(is,l);
         //result(is)(l) = std::move(S).choleskyDecomposition ().transpose (). inverse ();
         result(is)(l) = std::move(S).choleskyDecomposition (UpperRight). inverse ();
         for (int ifl = 0; ifl < nFL; ifl++)  {
            int iot = funcLMap(is)(l)(ifl);
            muSet(is)(iot).set(0.0);
            for (int jfl = 0; jfl < nFL; jfl++)  {
               muSet(is)(iot) += result(is)(l)(jfl,ifl) * org.getFuncL(is,l,jfl);
            }
         }
      }
   }
   return result;
}

SxArray<SxArray<SxVector<double> > > SxAtomicOrbitals::getOverlap () const
{
   int nSpecies = getNSpecies ();
   SxArray<SxArray<SxVector<double> > > result(nSpecies);
   for (int is = 0; is < nSpecies; is++)  {
      int lMax = getLMax(is);
      result(is).resize(lMax + 1);
      for (int l = 0; l <= lMax; l++)  {
         result(is)(l) = getOverlap (is, l);
      }
   }
   return result;
}

SxVector<double> SxAtomicOrbitals::getOverlap (int is, int l) const
{
   const SxVector<double> &r = radBasisPtr->radFunc(is);
   double logDr = radBasisPtr->logDr(is);
   int nFL = getFuncPerL (is,l);
   SxVector<double> result(nFL,nFL);
   for (int ifl = 0; ifl < nFL; ifl++)  {
      for (int jfl = ifl; jfl < nFL; jfl++)  {
         result(ifl,jfl) = result(jfl,ifl) = 
            (getFuncL(is,l,ifl) * getFuncL(is,l,jfl)*r*r*r).integrate(logDr);
      }
   }
   return result;
}

void SxAtomicOrbitals::refine (SxArray<int> &factor)
{
   SX_CHECK(radBasisPtr.getPtr () != NULL);
   const SxRadBasis &radBasis = *radBasisPtr;

   int nSpecies = getNSpecies ();
   SX_CHECK(factor.getSize() == nSpecies, factor.getSize(), nSpecies);
   SxArray<SxVector<double> > rad(nSpecies);
   SxArray<double> logDr(nSpecies);
   for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      SX_CHECK(factor(iSpecies) > 0, factor(iSpecies));
      if (factor(iSpecies) == 1) {
         cout << "Mesh for species " << (iSpecies + 1) 
              << " not refined." << endl;
         rad(iSpecies) = radBasis.radFunc(iSpecies);
         logDr(iSpecies) = radBasis.logDr(iSpecies);
      }
      else  {
         // --- interpolate phi on finer mesh
         cout << "Refining mesh for species " << (iSpecies + 1) 
              << " by factor " << factor(iSpecies) << "." << endl;
         double logDrOrig = radBasis.logDr(iSpecies);
         double r0 = radBasis.radFunc(iSpecies)(0);
         double r0New = r0;

         // --- define finer mesh
         int nrOrig = (int)radBasis.radFunc(iSpecies).getSize ();
         logDr(iSpecies) = logDrOrig / factor(iSpecies);
         int nExtra = (r0New < 0.) ? 0 : (int)(log(r0/r0New) / logDr(iSpecies));
         int nFine = (nrOrig-1)*factor(iSpecies) + 1 + nExtra;
         rad(iSpecies).resize(nFine);

         for (int i = 0; i < nFine; ++i)  {
            rad(iSpecies)(i) = r0 * exp(i * logDr(iSpecies));
         }

         // interpolate
         for (int iot = 0; iot < getNOrbTypes (iSpecies); iot++)  {
            // update rad,psi,logDr in psPot
            muSet(iSpecies)(iot) 
               = interpolateRad (muSet(iSpecies)(iot),
                     radBasis.radFunc(iSpecies)(0),
                     radBasis.logDr(iSpecies),
                     rad(iSpecies));
         }
      }
   }
   radBasisPtr = SxConstPtr<SxRadBasis>::create(rad,logDr);
   setBasis (radBasisPtr);
}

SxVector<double> SxAtomicOrbitals::compressWave(SxString &file, 
            int iState, 
            int iSpecies, 
            int l)
{
   SxPW waves;
   SxAtomicStructure structure;
   try {
      SxBinIO io (file, SxBinIO::BINARY_READ_ONLY);
      waves = SxPW (file, SxPW::InMemory);
      structure.read (io);
      io.close ();
   }
   catch (SxException e)  {
      e.print ();
      SX_QUIT;
   }

   SxStack<double> xData;
   SxStack<double> yData;
   SxVector<double> yRad;
   SxVector<double> result;

   SxGkBasis &GkBasis = *waves.getGkBasisPtr ();
   GkBasis.changeTau(structure);
   double naNorm = 1.0/(structure.getNAtoms(iSpecies));
   double mNorm = 1.0/(2.0*l+1.0);
   double spinNorm = 1.0/double(waves.getNSpin ());

   for (int ik = 0; ik < waves.getNk(); ik++)  {
      yRad.resize(GkBasis(ik).g2.getSize());
      yRad.set(0.0);
      for (int iSpin = 0; iSpin < waves.getNSpin (); iSpin++)  {
         for (int iAtom = 0; iAtom < structure.getNAtoms(iSpecies); iAtom++)  {
            for (int m = -l; m <= l; m++)  {
               SxVector<double> wave = sqrt(waves(iState,iSpin,ik).absSqr ());
               SxVector<SxComplex16> radial = wave
                  * GkBasis(ik).getPhaseFactors(iSpecies,iAtom).conj () //shift
                  * GkBasis(ik).getYlm(l,m) // project l
                  * SxYlm::getYlmNormFactor(l,m)
                  * sqrt(2.0 * structure.cell.volume/PI) 
                  * spinNorm
                  * naNorm 
                  * mNorm;
               yRad += radial.real ();
            }
         }
      }
      for (int i = 0; i < GkBasis(ik).g2.getSize(); i++)  {
         xData.push(sqrt(GkBasis(ik).g2(i)));
         yData.push(yRad(i));
      }
   }

   SxVector<double> xToFit(xData);
   SxVector<double> yToFit(yData);

   SxCubicSpline fit;
   SxRadialBasis radGBasis (0.0, GkBasis.getMaxGk(), 100, /*realSpace*/ false,
                            SxRadialBasis::Linear);
   if (l == 0) // ?(May have CUSP, no Mirror)?
      fit  = SxCubicSpline (
            xToFit,
            yToFit,
            radGBasis.getRadFunc(),
            SxCubicSpline::Natural,
            SxCubicSpline::MirrorPlane);
   else if (l & 1) // (l odd)
      fit  = SxCubicSpline (
            xToFit,
            yToFit,
            radGBasis.getRadFunc(),
            SxCubicSpline::Natural,
            SxCubicSpline::MirrorPoint);
   else // (l even and not zero)
      fit  = SxCubicSpline (
            xToFit,
            yToFit,
            radGBasis.getRadFunc(),
            SxCubicSpline::Natural,
            SxCubicSpline::MirrorPlane);
   SxVector<double> resultG = fit.getYFit ();
   resultG.setBasis(&radGBasis);
   resultG.auxData.is = iSpecies;
   resultG.auxData.ia = -1;
   resultG.auxData.n = -1;
   resultG.auxData.l = char(l);
   resultG.auxData.m = NONE_M;

   result = *radBasisPtr | resultG;

   result.auxData.is = iSpecies;
   result.auxData.ia = -1;
   result.auxData.n = -1;
   result.auxData.l = char(l);
   result.auxData.m = NONE_M;

   // ensure localization
   SxVector<double> rCut (radBasisPtr->radFunc(iSpecies).getSize ());
   rCut.set(10.0);
   result *= 1.0 / (1.0 + exp((radBasisPtr->radFunc(iSpecies) - rCut) / 1.0));

   result.normalize ();

   return result;
}

#ifdef USE_HDF5
void SxAtomicOrbitals::writeHDF5(const SxString &name)
{
   SxHDF5 file (name, SxHDF5::BINARY_CREATE);
   SxArray<SxArray<int> > funcPerL = getFuncPerL ();
   for (int is = 0; is < funcPerL.getSize (); is++)  {
      SxString isName = "Species-" + SxString(is);
      file.createGroup(isName);
      file.enterGroup(isName);
      SxString basisName = "radBasis";
      SxVector<double> radBasis = toVector(radBasisPtr->radFunc(is));
      file.writeVector(basisName, radBasis);
      for (int l = 0; l < funcPerL(is).getSize(); l++)  {
         SxString lName = "Angularmomentum-" + SxString(l);
         file.createGroup(lName);
         file.enterGroup(lName);
         for (int il = 0; il < funcPerL(is)(l); il++)  {
            SxString radName = "Radial-" + SxString(il);
            SxVector<double> radial = toVector(getFuncL(is,l,il));
            file.writeVector(radName, radial);
         }
         file.leaveGroup();
      }
      file.leaveGroup();
   }
}
#endif
