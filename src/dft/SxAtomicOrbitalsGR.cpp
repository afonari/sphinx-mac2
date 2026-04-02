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

#include <SxAtomicOrbitalsGR.h>
#include <SxCubicSpline.h>
#include <SxTextIO.h>
SxAtomicOrbitalsGR::SxAtomicOrbitalsGR ()
{
   // empty
}

SxAtomicOrbitalsGR::SxAtomicOrbitalsGR (
      const SxArray<SxArray<SxVector<double> > > &in,
      const SxConstPtr<SxRadialBasis> &radialBasisPtrIn,
      bool splineRepIn)
{
   // --- check dimensions
#  ifndef NDEBUG
   for (int is = 0; is < in.getSize(); is++)  {
      int nOrbTypes = (int)in(is).getSize();
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         SX_CHECK(in(is)(iot).auxData.is >= 0, in(is)(iot).auxData.is);
         SX_CHECK(in(is)(iot).auxData.n >= 0, (int)in(is)(iot).auxData.n);
         SX_CHECK(in(is)(iot).auxData.l >= 0, (int)in(is)(iot).auxData.l);
         if (splineRepIn) {
            SX_CHECK (in(is)(iot).getSize() == 4 * radialBasisPtrIn->getNElements(),
                      in(is)(iot).getSize(), 4 * radialBasisPtrIn->getNElements());
         }
         else  {   
            SX_CHECK (in(is)(iot).getSize() == radialBasisPtrIn->getNElements(),
                      in(is)(iot).getSize(), radialBasisPtrIn->getNElements());
         }
      }
   }
#  endif /* NDEBUG */

   radialBasisPtr = radialBasisPtrIn;
   muSet = in;
   splineRep = splineRepIn;

   createFuncLMap ();
}

SxAtomicOrbitalsGR::SxAtomicOrbitalsGR (const SxAtomicOrbitalsGR &in)
{
   (*this) = in;
}

SxAtomicOrbitalsGR::~SxAtomicOrbitalsGR ()
{
   // empty
}

void SxAtomicOrbitalsGR::operator= (const SxAtomicOrbitalsGR &in)
{
   radialBasisPtr = in.radialBasisPtr;
   muSet = in.muSet;
   splineRep = in.splineRep;
   funcLMap = in.funcLMap;
}

void SxAtomicOrbitalsGR::operator+= (const SxAtomicOrbitalsGR &in)
{
   SX_CHECK(splineRep == in.splineRep);
   SX_LOOP2(is,iot) {
      SX_CHECK(muSet(is)(iot).getSize() == in.muSet(is)(iot).getSize(),
               muSet(is)(iot).getSize(), in.muSet(is)(iot).getSize());
      muSet(is)(iot) += in.muSet(is)(iot);
   }
}

void SxAtomicOrbitalsGR::operator-= (const SxAtomicOrbitalsGR &in)
{
   SX_CHECK(splineRep == in.splineRep);
   SX_LOOP2(is,iot)  {
      SX_CHECK(muSet(is)(iot).getSize() == in.muSet(is)(iot).getSize(),
               muSet(is)(iot).getSize(), in.muSet(is)(iot).getSize());
      muSet(is)(iot) -= in.muSet(is)(iot);
   }
}

SxAtomicOrbitalsGR SxAtomicOrbitalsGR::operator+ (const SxAtomicOrbitalsGR &in) const
{
   SX_CHECK(splineRep == in.splineRep);
   SX_CHECK(radialBasisPtr->getNElements() == in.radialBasisPtr->getNElements());
   int nSpecies = (int)muSet.getSize();
   SX_CHECK(nSpecies == in.getNSpecies(), nSpecies, in.getNSpecies());
   SxArray<SxArray<SxVector<double> > > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes (is);
      SX_CHECK(nOrbTypes == in.getNOrbTypes(is),
               nOrbTypes, in.getNOrbTypes(is));
      resMuSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         resMuSet(is)(iot) = muSet(is)(iot) + in.muSet(is)(iot);
         resMuSet(is)(iot).auxData.is = in.muSet(is)(iot).auxData.is;
         resMuSet(is)(iot).auxData.n  = in.muSet(is)(iot).auxData.n;
         resMuSet(is)(iot).auxData.l  = in.muSet(is)(iot).auxData.l;
      }
   }
   SxAtomicOrbitalsGR result = SxAtomicOrbitalsGR(resMuSet, radialBasisPtr, splineRep);

   return result;
}

void SxAtomicOrbitalsGR::sumMPI() const
{
   ssize_t nSpecies = muSet.getSize();
   SxArray<SxArray<SxVector<double> > > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      ssize_t nOrbTypes = getNOrbTypes (is);
      resMuSet(is).resize(nOrbTypes);
      for (ssize_t iot = 0; iot < nOrbTypes; iot++) 
         SxLoopMPI::sum (muSet(is)(iot));
   }
}

SxAtomicOrbitalsGR SxAtomicOrbitalsGR::operator- (const SxAtomicOrbitalsGR &in) const
{
   SX_CHECK(splineRep == in.splineRep);
   SX_CHECK(radialBasisPtr->getNElements() == in.radialBasisPtr->getNElements());
   int nSpecies = (int)muSet.getSize();
   SX_CHECK(nSpecies == in.getNSpecies(), nSpecies, in.getNSpecies());
   SxArray<SxArray<SxVector<double> > > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes (is);
      SX_CHECK(nOrbTypes == in.getNOrbTypes(is),
               nOrbTypes, in.getNOrbTypes(is));
      resMuSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         resMuSet(is)(iot) = muSet(is)(iot) - in.muSet(is)(iot);
         resMuSet(is)(iot).auxData.is = in.muSet(is)(iot).auxData.is;
         resMuSet(is)(iot).auxData.n  = in.muSet(is)(iot).auxData.n;
         resMuSet(is)(iot).auxData.l  = in.muSet(is)(iot).auxData.l;
      }
   }
   SxAtomicOrbitalsGR result = SxAtomicOrbitalsGR(resMuSet, radialBasisPtr, splineRep);

   return result;
}

SxAtomicOrbitalsGR SxAtomicOrbitalsGR::operator* (double skalar) const
{
   SX_CHECK(radialBasisPtr->getNElements() > 0);
   int nSpecies = (int)muSet.getSize();
   SxArray<SxArray<SxVector<double> > > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes(is);
      resMuSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         resMuSet(is)(iot) = muSet(is)(iot) * skalar;
         resMuSet(is)(iot).auxData.is = muSet(is)(iot).auxData.is;
         resMuSet(is)(iot).auxData.n  = muSet(is)(iot).auxData.n;
         resMuSet(is)(iot).auxData.l  = muSet(is)(iot).auxData.l;
      }
   }
   SxAtomicOrbitalsGR result = SxAtomicOrbitalsGR(resMuSet, radialBasisPtr, splineRep);

   return result;
}

SxAtomicOrbitalsGR SxAtomicOrbitalsGR::operator/ (double skalar) const
{
   SX_CHECK(radialBasisPtr->getNElements() > 0);
   int nSpecies = (int)muSet.getSize();
   SxArray<SxArray<SxVector<double> > > resMuSet(nSpecies);
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = getNOrbTypes(is);
      resMuSet(is).resize(nOrbTypes);
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         resMuSet(is)(iot) = muSet(is)(iot) / skalar;
         resMuSet(is)(iot).auxData.is = muSet(is)(iot).auxData.is;
         resMuSet(is)(iot).auxData.n  = muSet(is)(iot).auxData.n;
         resMuSet(is)(iot).auxData.l  = muSet(is)(iot).auxData.l;
      }
   }
   SxAtomicOrbitalsGR result = SxAtomicOrbitalsGR(resMuSet, radialBasisPtr, splineRep);

   return result;
}

SxVector<double> & SxAtomicOrbitalsGR::operator() (int is, int iot)
{
   return muSet(is)(iot);
}

const SxVector<double> & SxAtomicOrbitalsGR::operator() (int is, int iot) const
{
   return muSet(is)(iot);
}

SxVector<double> SxAtomicOrbitalsGR::operator() (int is, int ia, int iot, int l, int m) const
{
   SxVector<double> result = muSet(is)(iot);
   // Set auxData
   result.setBasis(*radialBasisPtr);
   result.auxData.is = is;
   result.auxData.ia = ia;
   result.auxData.n  = char(iot);
   result.auxData.l  = char(l);
   result.auxData.m = char(m);

   return result;
}

int SxAtomicOrbitalsGR::getIOT (int is, int n, int l) const
{
   int iot = 0;
   while ((muSet(is)(iot).auxData.n != n)
         || (muSet(is)(iot).auxData.l != l))   {
      iot++;
      if (iot >= muSet(is).getSize()) SX_EXIT;
   }

   return iot;
}

int SxAtomicOrbitalsGR::getOrbitalIdx (int is, int ia, int iot, int l, int m, SxArray<SxQuantumNumbers> map) const
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

int SxAtomicOrbitalsGR::getNSpecies () const
{
   return (int)muSet.getSize ();
}

int SxAtomicOrbitalsGR::getNOrbTypes (int iSpecies) const
{
   return (int)muSet(iSpecies).getSize ();
}

int SxAtomicOrbitalsGR::getNOrbTypes () const
{
   int result = 0;
   for(int is = 0; is < getNSpecies (); is++)   
      result += getNOrbTypes(is);
      
   return result;
}

void SxAtomicOrbitalsGR::set(double val)
{

   for (int iSpecies = 0; iSpecies < muSet.getSize(); iSpecies++)   {
      for (int iot = 0; iot < muSet(iSpecies).getSize(); iot++)   {
         muSet(iSpecies)(iot).set(val);
      }
   }
}

SxArray<SxQuantumNumbers> SxAtomicOrbitalsGR::getReducedOrbitalMap () const
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
SxAtomicOrbitalsGR::getOrbitalMap (const SxAtomicStructure &structure) const
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

int SxAtomicOrbitalsGR::getNOrbitals (const SxAtomicStructure &structure) const
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

void SxAtomicOrbitalsGR::createFuncLMap ()
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

SxRadialBasis::TPsi &SxAtomicOrbitalsGR::getFuncL (int is, int l, int ifl)
{
   int iot = funcLMap(is)(l)(ifl);
   return muSet(is)(iot);
}

const SxRadialBasis::TPsi &SxAtomicOrbitalsGR::getFuncL (int is, int l, int ifl) const
{
   int iot = funcLMap(is)(l)(ifl);
   return muSet(is)(iot);
}

int SxAtomicOrbitalsGR::getLMax () const
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

int SxAtomicOrbitalsGR::getLMax (const int iSpecies) const
{
   int nOrbTypes = (int)muSet(iSpecies).getSize ();
   int result = 0;
   for (int iot = 0; iot < nOrbTypes; iot++)   {
         if (result <  muSet(iSpecies)(iot).auxData.l)
            result = muSet(iSpecies)(iot).auxData.l;
   }
   return result;
}

int SxAtomicOrbitalsGR::getFuncPerL (int is, int l) const
{
   return (int)funcLMap(is)(l).getSize();
}

void SxAtomicOrbitalsGR::print (const SxString &fileIn) const
{
   int nSpecies = (int)muSet.getSize ();
   for(int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = (int)muSet(is).getSize();
      for(int iot = 0; iot < nOrbTypes; iot++)   {
         SxString file = fileIn + SxString(is) + SxString(iot) + ".dat";
         SxVecRef<double> vec;
         if (splineRep) vec = toVec(is,iot);
         else           vec = muSet(is)(iot);
         SxTextIO (file).writeXYPlot(radialBasisPtr->getRadFunc(), vec);
      }
   }
}

void SxAtomicOrbitalsGR::setBasis (const SxConstPtr<SxRadialBasis> &radialBasisPtrIn)
{
   radialBasisPtr = radialBasisPtrIn;
   for (int is = 0; is < muSet.getSize (); ++is)  {
      for (int iot = 0; iot < muSet(is).getSize (); ++iot)  {
         muSet(is)(iot).setBasis (radialBasisPtr.getConstPtr ());
      }
   }
}


double SxAtomicOrbitalsGR::getNormSqr(int is, int iot) const
{
   SxVecRef<double> vec;
   if (splineRep) vec = toVec(is,iot);
   else           vec = muSet(is)(iot);
   double result = radialBasisPtr->integrate(vec * vec);
   return result;
}

double SxAtomicOrbitalsGR::getNormSqrSum() const
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
double SxAtomicOrbitalsGR::dot (const SxAtomicOrbitalsGR &orbitalsIn) const
{
   double result = 0.0;
   int nSpecies = getNSpecies ();
   SX_CHECK(nSpecies == orbitalsIn.getNSpecies());

   for (int is = 0; is < nSpecies; is++)  {
      SX_CHECK (getNOrbTypes(is) == orbitalsIn.getNOrbTypes (is));
      int nOrbTypes = (int)muSet(is).getSize ();
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         result += dot(orbitalsIn,is,iot,is,iot);
      }
   }

   return result;
}

double SxAtomicOrbitalsGR::dot (const SxAtomicOrbitalsGR &orbitalsIn, int is, int iot, int js, int jot) const
{
   SX_CHECK(radialBasisPtr == orbitalsIn.radialBasisPtr);
   SxVecRef<double> vec1;
   if (splineRep) vec1 = toVec(is,iot);
   else           vec1 = muSet(is)(iot);
   SxVecRef<double> vec2;
   if (orbitalsIn.splineRep) vec2 = orbitalsIn.toVec(js,jot);
   else           vec2 = orbitalsIn.muSet(js)(jot);
   double result = radialBasisPtr->integrate(vec1 * vec2);
   return result;
}

double SxAtomicOrbitalsGR::sum (const SxAtomicOrbitalsGR &orbitalsIn) const
{
   double result = 0.0;
   int nSpecies = getNSpecies ();
   SX_CHECK(nSpecies == orbitalsIn.getNSpecies());

   for (int is = 0; is < nSpecies; is++)  {
      SX_CHECK (getNOrbTypes(is) == orbitalsIn.getNOrbTypes (is));
      int nOrbTypes = (int)muSet(is).getSize ();
      for (int iot = 0; iot < nOrbTypes; iot++)  {
         result += sum(orbitalsIn,is,iot,is,iot);
      }
   }

   return result;
}

double SxAtomicOrbitalsGR::sum (const SxAtomicOrbitalsGR &orbitalsIn, int is, int iot, int js, int jot) const
{
   SX_CHECK(radialBasisPtr == orbitalsIn.radialBasisPtr);
   SxVecRef<double> vec1;
   if (splineRep) vec1 = toVec(is,iot);
   else           vec1 = muSet(is)(iot);
   SxVecRef<double> vec2;
   if (orbitalsIn.splineRep) vec2 = orbitalsIn.toVec(js,jot);
   else           vec2 = orbitalsIn.muSet(js)(jot);
   double result = (vec1 * vec2).sum ();
   return result;
}

void SxAtomicOrbitalsGR::toSpline () const
{
   SX_CHECK (!splineRep);
   if (!splineRep)  {
      const SxVecRef<double> &radGFunc = radialBasisPtr->getRadFunc();
      int nSpecies = getNSpecies ();
      for (int is = 0; is < nSpecies; is++)   {
         int nOrbtypes = getNOrbTypes (is);
         for (int iot = 0; iot < nOrbtypes; iot++)   {
            SX_CHECK(is == muSet(is)(iot).auxData.is,
                     is, muSet(is)(iot).auxData.is);
            int ia = muSet(is)(iot).auxData.ia;
            int n = muSet(is)(iot).auxData.n;
            int l = muSet(is)(iot).auxData.l;
            int m = muSet(is)(iot).auxData.m;
            SxCubicSpline spline;
            if (l < 2) 
               spline = SxCubicSpline (
                     radGFunc, 
                     muSet(is)(iot), 
                     SxCubicSpline::NaturalHermite);
            else
               spline = SxCubicSpline (
                     radGFunc, 
                     muSet(is)(iot), 
                     SxCubicSpline::Hermite);
            muSet(is)(iot) = spline.getSpline ();
            muSet(is)(iot).setBasis(&*radialBasisPtr);
            muSet(is)(iot).auxData.is = is;
            muSet(is)(iot).auxData.ia = ia;
            muSet(is)(iot).auxData.n = char(n);
            muSet(is)(iot).auxData.l = char(l);
            muSet(is)(iot).auxData.m = char(m);
         }
      }
      splineRep = true;
   }
}

void SxAtomicOrbitalsGR::toVec () const
{
   SX_CHECK (splineRep);
   if(splineRep)  {
      const SxVecRef<double> &radGFunc = radialBasisPtr->getRadFunc();
      int nSpecies = getNSpecies ();
      for (int is = 0; is < nSpecies; is++)   {
         int nOrbtypes = getNOrbTypes (is);
         for (int iot = 0; iot < nOrbtypes; iot++)   {
            SX_CHECK(is == muSet(is)(iot).auxData.is,
                     is, muSet(is)(iot).auxData.is);
            int ia = muSet(is)(iot).auxData.ia;
            int n = muSet(is)(iot).auxData.n;
            int l = muSet(is)(iot).auxData.l;
            int m = muSet(is)(iot).auxData.m;
            SxCubicSpline spline (radGFunc, muSet(is)(iot));
            muSet(is)(iot) = spline.getY (radGFunc);
            muSet(is)(iot).setBasis(&*radialBasisPtr);
            muSet(is)(iot).auxData.is = is;
            muSet(is)(iot).auxData.ia = ia;
            muSet(is)(iot).auxData.n = char(n);
            muSet(is)(iot).auxData.l = char(l);
            muSet(is)(iot).auxData.m = char(m);
         }
      }
      splineRep = false;
   } 
}

SxVector<double> SxAtomicOrbitalsGR::toVec (int is, int iot) const
{
   SX_CHECK (splineRep);
   const SxVecRef<double> &radFunc = radialBasisPtr->getRadFunc();
   SxCubicSpline spline (radFunc, muSet(is)(iot));
   SxVector<double> result = spline.getY (radFunc);
   return result;
}

void SxAtomicOrbitalsGR::normalize ()
{
   int nSpecies = (int)muSet.getSize();
   for (int is = 0; is < nSpecies; is++)   {
      int nOrbTypes = (int)muSet(is).getSize ();
      for (int iot = 0; iot < nOrbTypes; iot++)   {
         double norm = getNormSqr(is,iot);
         if (norm > 1e-12) muSet(is)(iot) /= sqrt(norm);
         else {
            cout << "WARNING: Try to normalize zero norm orbital." << endl;
            cout << "Norm is " << norm << endl;
            cout << "Orbital is set to zero." << endl;
            muSet(is)(iot).set(0.0);
         }
      }
   }
}

SxArray<SxArray<SxVector<double> > > SxAtomicOrbitalsGR::orthogonalize ()
{
   SX_CHECK(muSet.getSize() == funcLMap.getSize(),
            muSet.getSize(), funcLMap.getSize());
   // Gram-Schmidt via Cholesky
   SxAtomicOrbitalsGR org = 1.0 * *this;
   int nSpecies = (int)muSet.getSize ();
   SxArray<SxArray<SxVector<double> > > result(nSpecies);
   for (int is = 0; is < nSpecies; is++)  {
      int lMax = getLMax(is);
      result(is).resize(lMax + 1);
      for (int l = 0; l <= lMax; l++)  {
         int nFL = getFuncPerL (is,l);
         if (nFL == 0) continue;
         SxVector<double> S = getOverlap(is,l);
         result(is)(l) = std::move(S).choleskyDecomposition ().transpose (). inverse ();
         for (int ifl = 0; ifl < nFL; ifl++)  {
            int iot = funcLMap(is)(l)(ifl);
            muSet(is)(iot).set(0.0);
            for (int jfl = 0; jfl < nFL; jfl++)  {
               muSet(is)(iot) += result(is)(l)(jfl,ifl) * org.getFuncL(is,l,jfl);
            }
         }
      }
   }

   // return CholeskyDecompositions
   return result;
}

SxArray<SxArray<SxVector<double> > > SxAtomicOrbitalsGR::getOverlap () const
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

SxVector<double>  SxAtomicOrbitalsGR::getOverlap (int is, int l) const
{
   int nFL = getFuncPerL (is,l);
   SxVector<double> result(nFL,nFL);
   for (int ifl = 0; ifl < nFL; ifl++)  {
      int iot = funcLMap(is)(l)(ifl);
      SxVecRef<double> vec1;
      if (splineRep) vec1 = toVec(is,iot);
      else           vec1 = muSet(is)(iot);
      for (int jfl = ifl; jfl < nFL; jfl++)  {
         int jot = funcLMap(is)(l)(jfl);
         SxVecRef<double> vec2;
         if (splineRep) vec2 = toVec(is,jot);
         else           vec2 = muSet(is)(jot);
         result(ifl,jfl) = result(jfl,ifl) = 
            radialBasisPtr->integrate(vec1 * vec2);
      }
   }
   return result;
}

void SxAtomicOrbitalsGR::orthogonalizeOn(SxAtomicOrbitalsGR &basis)
{
   //Gram-Schmidt ortho
   bool switchThis = false;
   bool switchBasis = false;
   int nSpecies = (int)muSet.getSize ();
   SX_CHECK(nSpecies == basis.muSet.getSize());
   
   if (isSpline()) { 
      toVec ();
      switchThis = true;   
   }

   if (basis.isSpline())  {
      basis.toVec ();
      switchBasis = true;
   }

   for (int is = 0; is < nSpecies; is++)  {
      int lMax = getLMax(is);
      SX_CHECK(lMax <= basis.getLMax(is));
      for (int l = 0; l <= lMax; l++)  {
         int nFL = getFuncPerL(is,l);
         int nFLbasis = basis.getFuncPerL(is,l);
         for (int ifl = 0; ifl < nFL; ifl++)  {
            int iot = funcLMap(is)(l)(ifl);
            SxVector<double> vec = muSet(is)(iot);
            for (int jfl = 0; jfl < nFLbasis; jfl++)  {
               int jot = basis.funcLMap(is)(l)(jfl);
               muSet(is)(iot) -= radialBasisPtr->integrate(vec * basis.muSet(is)(jot))
                               / basis.getNormSqr(is,jot) * basis.muSet(is)(jot);
            }
         }
      }
   }

   if (switchThis) toSpline ();
   if (switchBasis) basis.toSpline ();
}

void SxAtomicOrbitalsGR::rotate(const SxArray<SxArray<SxVector<double> > > &rotMat)
{
   SX_CHECK(muSet.getSize() == funcLMap.getSize(),
            muSet.getSize(), funcLMap.getSize());
   // Gram-Schmidt via Cholesky
   SxAtomicOrbitalsGR org = 1.0 * *this;
   int nSpecies = (int)muSet.getSize ();
   SX_CHECK (nSpecies == rotMat.getSize(),nSpecies, rotMat.getSize());
   for (int is = 0; is < nSpecies; is++)  {
      int lMax = getLMax(is);
      SX_CHECK (lMax + 1 == rotMat(is).getSize(), lMax + 1, rotMat(is).getSize());
      for (int l = 0; l <= lMax; l++)  {
         int nFL = getFuncPerL (is,l);
         SX_CHECK (nFL == rotMat(is)(l).getNCols (), nFL, rotMat(is)(l).getNCols ());
         for (int ifl = 0; ifl < nFL; ifl++)  {
            int iot = funcLMap(is)(l)(ifl);
            muSet(is)(iot).set(0.0);
            for (int jfl = 0; jfl < nFL; jfl++)  {
               muSet(is)(iot) += rotMat(is)(l)(jfl,ifl) * org.getFuncL(is,l,jfl);
            }
         }
      }
   }
}

