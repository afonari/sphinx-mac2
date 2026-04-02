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


/* Implementation of n band k dot p perturbation theory.
*/

#include <SxKdotP.h>
#include <SxGBasis.h>
#include <SxFSAction.h>
#include <SxFileIO.h>
#include <SxFileParser.h>
#include <SxIO.h>
#include <SxEigensystem.h>
#include <sobol.hpp>

namespace Timer {
   enum KPNxNHamTimer {
      sndDer,
      fstDer,
      NxNham,
      init,
      evalTree,
      EvalTreeBS,
      EvalByteCode,
      EvalByteCodeBS,
      derivatives,
      kk2s,
      extField,
      complexNr,
      doubleNr,
      known,
      constant,
      FitBandsTotal,
      FitBands,
      FitBandsDiag
   };
}

SX_REGISTER_TIMERS (Timer::KPNxNHamTimer)
{
   using namespace Timer;
   regTimer (sndDer,         "second derivative");
   regTimer (fstDer,         "first derivative");
   regTimer (NxNham,         "Hamiltonian routine");
   regTimer (init,           "Hamiltonian init");
   regTimer (evalTree,       "Evaluate tree");
   regTimer (EvalTreeBS,     "Evaluate tree BS");
   regTimer (EvalByteCode,   "Evaluate byte code");
   regTimer (EvalByteCodeBS, "byte code BS");
   regTimer (kk2s,           "kk2s");
   regTimer (extField,       "ext");
   regTimer (complexNr,      "complex");
   regTimer (doubleNr,       "double");
   regTimer (known,          "known");
   regTimer (constant,       "constant");
   regTimer (FitBandsTotal,  "fit bands total");
   regTimer (FitBands,       "fit bands");
   regTimer (FitBandsDiag,   "fit bands: diag");
}

#define UNIQUE_PAUSE(x)  UNIQUE_TIMER_NAME(x)
#define SX_CLOCK_PAUSE(ID) SxClockPause UNIQUE_PAUSE(__LINE__) (SxTimer::getTimerId (ID))
class SxClockPause {
   int id;
   public:
      inline SxClockPause (int idIn) : id(idIn)  {
         if (SxTimer::getGlobalTimer ().isIdle (id))
            id = -1;
         else
            SxTimer::getGlobalTimer ().stop (id);
      }
      inline ~SxClockPause ()
      {
         if (id >= 0)  {
            SxTimer &timer = SxTimer::getGlobalTimer ();
            timer.start (id);
            timer.nCalls(id)--;
         }
      }
};

SxKdotP::SxKdotP (const SxPWSet &waves,
                        const RhoR &rhoRIn)
   : SxHamiltonian (),
     SxRho (rhoRIn),
     wavesPtr (&waves)
// SxKdotP
{
   SX_CHECK (rBasisPtr);
   sdCache.resize(6);
   fdCache.resize(3);
   
   gBasisPtr = &(rBasisPtr->getGBasis ());
   
   D.resize(1);
   
}

PrecEnergy SxKdotP::getETrial()
{
   return eTrial;
}

/*
// TODO: This is an ugly routine. Does 2 different things when called, distinguished by
// the counter variable. It is required for optical spectra that rely on the Hamiltonian
// but not used for the electronic structure calculation itself at all.
// Problem is that all communication with the Hamiltonian from anywhere else
// (e.g. a tool) is done via the HamPtr, so whatsoever is required needs a virtual
// function in SxHamiltonian.h. Maybe it is reasonable to allow for some arbitrary
// communication channel in SxHamiltonian that can be used to transfer various
// data between a Hamiltonian and other tools & add-ons?
void SxKdotP::formDerivatives()
{
   SX_EXIT; // CF 2021-03-24: this won't work with bytecode implementation
   if (kDerivFirstRound)  {
      accurateInterfaces = false;
      cout << "Switch to simplified interfaces for optical spectra..." << endl;
      kDerivFirstRound = false;
   } else
   {
      cout << "form k-derivative of Hamiltonian" << endl;
      derivatives = true;
      // --- replace with templates to ensure no parameters are substituted twice
      if (inFile.contains("k2") > 0) inFile = inFile.substitute("k2", "###K2###");
      if (inFile.contains("kx2") > 0) inFile = inFile.substitute("kx2", "###KX2###");
      if (inFile.contains("ky2") > 0) inFile = inFile.substitute("ky2", "###KY2###");
      if (inFile.contains("kz2") > 0) inFile = inFile.substitute("kz2", "###KZ2###");
      if (inFile.contains("kxy") > 0) inFile = inFile.substitute("kxy", "###KXY###");
      if (inFile.contains("kxz") > 0) inFile = inFile.substitute("kxz", "###KXZ###");
      if (inFile.contains("kyz") > 0) inFile = inFile.substitute("kyz", "###KYZ###");
      if (inFile.contains("kx") > 0) inFile = inFile.substitute("kx", "###KX###");
      if (inFile.contains("ky") > 0) inFile = inFile.substitute("ky", "###KY###");
      if (inFile.contains("kz") > 0) inFile = inFile.substitute("kz", "###KZ###");
      if (inFile.contains("k") > 0) inFile = inFile.substitute("k", "###K###");
      // --- replace with k-derivatives, multiplied by light polarisation a
      if (inFile.contains("###K2###")) inFile = inFile.substitute("###K2###", "(2*kx*ax+2*ky*ay+2*kz*az)");
      if (inFile.contains("###KX2###")) inFile = inFile.substitute("###KX2###", "2*kx*ax");
      if (inFile.contains("###KY2###")) inFile = inFile.substitute("###KY2###", "2*ky*ay");
      if (inFile.contains("###KZ2###")) inFile = inFile.substitute("###KZ2###", "2*kz*az");
      if (inFile.contains("###KXY###")) inFile = inFile.substitute("###KXY###", "(ky*ax+kx*ay)");
      if (inFile.contains("###KXZ###")) inFile = inFile.substitute("###KXZ###", "(kz*ax+kx*az)");
      if (inFile.contains("###KYZ###")) inFile = inFile.substitute("###KYZ###", "(kz*ay+ky*az)");
      if (inFile.contains("###K###")) inFile = inFile.substitute("###K###", "(ax+ay+az)");
      if (inFile.contains("###KX###")) inFile = inFile.substitute("###KX###", "(ax)");
      if (inFile.contains("###KY###")) inFile = inFile.substitute("###KY###", "(ay)");
      if (inFile.contains("###KZ###")) inFile = inFile.substitute("###KZ###", "(az)");
      cout << "build updated tree: " << endl;
      SxString hamString = ((inFile.right("Hamiltonian")).right("=[")).left("];");
      SxList<SxString> cols = hamString.tokenize(']');
      int i, j;
      for (i = 0; i < nComp; i++)  {
         cols(i) = cols(i).right("[");
         SxList<SxString> rows = cols(i).tokenize(',');
         int nRows = (int)rows.getSize();
         // --- split columns in single elements
         for (j = 0; j < nRows; j++)  {
            expression(i)(j).resize(0);
            // --- evaluate single element
            buildTree(i, j, rows(j), inFile);
         }
      }
      cout << "success." << endl;
   }
}
*/

PsiG SxKdotP::operator* (const SxVecRef<PrecCoeffG> &psiG)
{
   // --- test: show <psiI|H|psiJ>
   SX_CHECK (psiG.getBasisPtr());

   const SxGBasis &G = psiG.getBasis<SxGBasis> ();

   if (psiG.getNCols () > 1)  {
      PsiG res(psiG.getNRows (), psiG.getNCols ());
      res.auxData = psiG.auxData;
      SX_LOOP (iState)
         res.colRef (iState) = G | hamNxN (psiG.colRef (iState));
      return res;
   } // else
   return (G | hamNxN (psiG));
}

void SxKdotP::writeMeshASCII (const SxString &name, const PsiR &data) const
{
   
   const SxRBasis &R = *rBasisPtr;
   SxMesh3D mesh = R.getMesh();
   SxString fileName = name + ".dat";
   ofstream fileN ((fileName).ascii (), fstream::out | fstream::trunc);
   for (int x = 0; x < mesh(0); x++)  {
      for (int y = 0; y < mesh(1); y++)  {
         for (int z = 0; z < mesh(2); z++)  {
            int idx = (int)mesh.getMeshIdx(x, y, z, SxMesh3D::Positive);
            fileN << data(idx).re << endl;
         }
      }
   }
   cout << name << " written to " << fileName << endl;
}

SxVector<SxComplex16>
SxKdotP::hMatrixBS(const SxVector3<double> &kVec,
                   const SxVector<double> &strain,
                   double vPol, double vExtern, double vCharge,
                   const SxVecRef<double> &params,
                   SxStack<SxComplex16> &stack) const
{
   SxComplex16 e;

   SxVector<SxComplex16> hamiltonian(nComp, nComp);
   SX_CLOCK(Timer::EvalByteCodeBS);
   for (int jComp = 0; jComp < nComp; jComp++)  {
      for (int iComp = 0; iComp < nComp; iComp++)  {
         e = evaluateByteCode(byteCodeHam(iComp, jComp), kVec,
                              strain, vPol, vExtern, vCharge,
                              params, stack);
         hamiltonian(iComp, jComp) = e;
      }
   }

   return hamiltonian;
}

// --- provide potential landscape of all bands throughout cell
void SxKdotP::potentialLandscape (SxString outFile, SxVector3<double> kvals)
{

   int iBand, iPar;
   // clean up output file
   SX_OUTPUT_FILE (outFile);
   const SxRBasis &R = *rBasisPtr;
   SxMesh3D mesh = R.getMesh();
   
   int xMax = mesh(0);
   int yMax = mesh(1);
   int zMax = mesh(2);
   SxVector3<int> r;
//   zero3(0) = 0;
//   zero3(1) = 0;
//   zero3(2) = 0;

   SxVector<double> params(nParam);
   SxVector<double> strain(6);
   SxStack<SxComplex16> workStack;
   SxVector<SxComplex16> hMat;
   SxVector<double> eigVals;

   // --- setup parameters for bytecode evaluation
//   ssize_t coordIdx = rBasisPtr->fft3d.mesh.getMeshIdx(r, SxMesh3D::Positive);
   SxFileIO outio;
   SX_OUTPUT_FILE (outFile); 
   try  {
      outio.open (outFile, "a");
   } catch (const SxException &e) {
      e.print ();
      SX_EXIT;
   }
   int x,y,z;
   if ((plLayer(0) == -1) && (plLayer(1) == -1) && (plLayer(2) == -1))  {
      cout  << "Full 3D potential landscape requested. This may take a while..." << endl;
   }
   if ((plLayer(0) > -1) || (plLayer(1) > -1) || (plLayer(2) > -1))  {
      cout << "Potential landscape required for " << plLayer << endl;
   }
   for (x = 0; x < xMax; x++)  {
      if ((plLayer(0) == -1) || (plLayer(0) == x))  { // compute and plot for all (-1) or only one layer
      r(0) = x;
      for (y = 0; y < yMax; y++)  {
         if ((plLayer(1) == -1) || (plLayer(1) == y))  {
         r(1) = y;
         for (z = 0; z < zMax; z++)  {
            if ((plLayer(2) == -1) || (plLayer(2) == z))  {
            ssize_t coordIdx = rBasisPtr->fft3d.mesh.getMeshIdx(r, SxMesh3D::Positive);
            r(2) = z;
            // parameter values
            for (iPar = 0; iPar < nParam; ++iPar)
               params(iPar) = parameters(iPar)(coordIdx);
            // strain
            if (eIJ.getSize () > 0)
               SX_LOOP(i) strain(i) = eIJ(i)(coordIdx);
            // potentials
            double vPol    = (vP.getSize () > 0  ) ? vP  (coordIdx) : 0.,
                   vExtern = (vExt.getSize () > 0) ? vExt(coordIdx) : 0.,
                   vCharge = (chargePotential.getSize () > 0)
                           ? chargePotential(coordIdx) : 0.;
         
            hMat = hMatrixBS(kvals, strain, vPol, vExtern, vCharge, params,
                   workStack);
            eigVals = symEigenvalues (hMat);
            for (iBand = 0; iBand < nComp; iBand++)  
               (ostream&)outio << "   " << (eigVals(iBand) * HA2EV);
            (ostream&)outio << endl;
            
         }}
      }}
   }}
   outio.close ();
   cout << "+------- Potential Landscape calculation finished. -------" << endl;
}
// --- provide bandstructure between given k's
void SxKdotP::showBandstructure (SxString outFile, SxVector3<int> r)
{

   int iStep, iBand;
   // clean up output file
   SX_OUTPUT_FILE (outFile);

   // --- setup parameters for bytecode evaluation
   ssize_t coordIdx = rBasisPtr->fft3d.mesh.getMeshIdx(r, SxMesh3D::Positive);
   // parameter values
   SxVector<double> params(nParam);
   for (int iPar = 0; iPar < nParam; ++iPar)
      params(iPar) = parameters(iPar)(coordIdx);
   // strain
   SxVector<double> strain(6);
   if (eIJ.getSize () > 0)
      SX_LOOP(i) strain(i) = eIJ(i)(coordIdx);
   // potentials
   double vPol    = (vP.getSize () > 0  ) ? vP  (coordIdx) : 0.,
          vExtern = (vExt.getSize () > 0) ? vExt(coordIdx) : 0.,
          vCharge = (chargePotential.getSize () > 0)
                  ? chargePotential(coordIdx) : 0.;
         
   SxFileIO outio;
   SX_OUTPUT_FILE (outFile); 
   try  {
      outio.open (outFile, "a");
   } catch (const SxException &e) {
      e.print ();
      SX_EXIT;
   }
   SxStack<SxComplex16> workStack;
   for (iStep = 0; iStep < bsPts.getSize(); iStep++) {
      SxVector<SxComplex16> hMat
         = hMatrixBS(bsPts(iStep), strain, vPol, vExtern, vCharge, params,
                     workStack);
      (ostream&)outio << bsPts(iStep)(0) << "  " << bsPts(iStep)(1) << "  "
                      << bsPts(iStep)(2);
      SxVector<double> eigVals = symEigenvalues (hMat);
      for (iBand = 0; iBand < nComp; iBand++)
         (ostream&)outio << "   " << (eigVals(iBand) * HA2EV);
      (ostream&)outio << endl;
   }
   outio.close ();
   cout << "+------- Band structure calculation finished. -------" << endl;
}

// --- the Hamiltonian itself
PsiR SxKdotP::hamNxN (const SxVecRef<PrecCoeffG> &psiG)
{
   SX_CLOCK(Timer::NxNham);
   SX_CHECK (psiG.getBasisPtr());
   const SxGBasis &G = psiG.getBasis<SxGBasis> ();
   const SxRBasis &R = *rBasisPtr;
   PsiR psiRFull = (R | psiG);
   int iComp, jComp, i;
   PsiR elem;

   // initialize preconditioner and potential from input charges
   if (D.getSize() == 1)  {
      SX_CLOCK(Timer::init);
      gZero.resize(G.ng);
      gOne.resize(G.ng);
      gZero.set (0.);
      gOne.set (1.);
      // --- resize fdC, to be used by hamElem(i,j,prec)
      for (i = 0; i < 3; i++)  {
         fdCache(i) = gZero;
      }

      // second derivatives in preconditioner
      D.resize(psiG.getSize());
      D.auxData = psiG.auxData;
      for (i = 0; i < 3; i++)  {
         sdCache(i).unref ();
         sdCache(i) = G.gVec.colRef(i).sqr ();
         sdCache(i+3).unref ();
         sdCache(i+3) = gZero;
      }
      for (iComp = 0; iComp < nComp; iComp++)  {
         D.compRef (iComp) = evaluateTreePrecond(iComp,iComp, 0).abs();
      }
   }   // timings: 0.0%

   PsiR result (psiRFull.getSize());
   result.auxData = psiRFull.auxData;
   result.set (0.);
   iMult = 0;
   for (iComp = 0; iComp < nComp; iComp++)  {
      psiR.unref ();
      psiR = psiRFull.compRef(iComp);
      psiGcomp.unref ();
      psiGcomp = psiG.compRef(iComp);
      for (i = 0; i < 3; i++)
         fdCache(i).unref ();
      for (i = 0; i < 6; i++)
         sdCache(i).unref ();


      for (jComp = 0; jComp < nComp; jComp++)  {
         // SX_CLOCK(Timer::evalTree);
         // result.compRef (jComp) += evaluateTree(jComp,iComp, 0);
         SX_CLOCK(Timer::EvalByteCode);
         result.compRef (jComp) += evaluateByteCode(byteCodeHam(jComp,iComp));
      }
   }
   return result;
}

const PsiR& SxKdotP::fdC(int i)
{
   SX_CHECK (i >= 0 && i < 3);
   if (fdCache(i).getSize () == 0)  {
      SX_CLOCK_PAUSE(Timer::evalTree);
      SX_CLOCK_PAUSE(Timer::EvalByteCode);
      SX_CLOCK_PAUSE(Timer::kk2s);
      SX_CLOCK(Timer::fstDer);
      const SxGBasis &G = psiGcomp.getBasis<SxGBasis> ();
      // -I d/dx_i = -I * (-I G_i) = -G_i
      fdCache(i) = *rBasisPtr | (-wgt(i) * G.gVec.colRef(i) * psiGcomp);
   }
   return fdCache(i);
}

const PsiRef& SxKdotP::sdC(int i, int j)
{
   if (i == j)  {
      if (sdCache(i).getSize () == 0)  {
         SX_CLOCK_PAUSE(Timer::evalTree);
         SX_CLOCK_PAUSE(Timer::EvalByteCode);
         SX_CLOCK_PAUSE(Timer::kk2s);
         SX_CLOCK(Timer::sndDer);
         const PsiG &G = psiGcomp.getBasis<SxGBasis> ().gVec;
         SxVector<double> Gii = G.colRef (i).sqr ();
         Gii *= sqr(wgt(i));
         sdCache(i) = *rBasisPtr | (Gii * psiGcomp);
      }
      return sdCache(i);
   }
   int iCache = -1;
   switch (3*i + j) {
      case 1: // 0,1 or 1,0
      case 3: iCache = 3; break;
      case 2: // 0,2 or 2,0
      case 6: iCache = 4; break;
      case 5: // 1,2 or 2,1
      case 7: iCache = 5; break;
      default: SX_EXIT;
   }
   if (sdCache(iCache).getSize () == 0)  {
      SX_CLOCK_PAUSE(Timer::evalTree);
      SX_CLOCK_PAUSE(Timer::EvalByteCode);
      SX_CLOCK_PAUSE(Timer::kk2s);
      SX_CLOCK(Timer::sndDer);
      const PsiG &G = psiGcomp.getBasis<SxGBasis> ().gVec;
      SxVector<double> Gij = G.colRef (i) * G.colRef (j);
      Gij *= wgt(i) * wgt(j);
      sdCache(iCache) = *rBasisPtr | (Gij * psiGcomp);
   }
   return sdCache(iCache);
}

const SxSymbolTable *
SxKdotP::getHamiltonianGroup (const SxSymbolTable *table)
{
   SX_CHECK (table);
   SxSymbolTable *hGroup = NULL;
   try  { hGroup = table->getGroup("kpHamiltonian"); }
   catch (const SxException&) { /*EXIT;*/ }
   return hGroup;
}


SxVecRef<double> SxKdotP::parameters(int iParam)
{
   SxMeshR p;
   int iMat;
   if (!speedOpt)  {
      p = matParam(0, iParam) * materials(0);
      if (bowParam(0, iParam) != 0)
            p -= bowParam(0, iParam)
               * (1. - materials(0)) * materials(0);
      for (iMat = 1; iMat < nMat; iMat++)  {
         p = p + matParam(iMat, iParam) * materials(iMat);
         if (bowParam(iMat, iParam) != 0)
            p -= bowParam(iMat, iParam)
               * (1. - materials(iMat)) * materials(iMat);
      }
      return p;
   }
   else
      return pMem(iParam);
}

int SxKdotP::firstOutsideBracket(SxString expr, char op)
{
   int pos = -1;
   int bLevel = 0;
   int i;

   for (i = 0; i < expr.getSize(); i++)  {
      if (expr(i) == '(') bLevel++;
      if (expr(i) == ')') bLevel--;
      if ((expr(i) == op) && (bLevel == 0) && (pos == -1))
         pos = i;
   }

   return pos;
}

SxString SxKdotP::containsUnknown(SxList<SxString> inList)
{
   SxString unknownElem = "###";
   int iElem;
   for (iElem = 0; iElem < inList.getSize(); iElem++)  {
      SxString elem = inList(iElem);
      if (
            ( // --- is one of the known operators
                (elem == "k2") ||
                (elem == "kx2") ||
                (elem == "ky2") ||
                (elem == "kz2") ||

                (elem == "kxy") ||
                (elem == "kxz") ||
                (elem == "kyz") ||

                (elem == "k") ||
                (elem == "kx") ||
                (elem == "ky") ||
                (elem == "kz") ||

                (elem == "e") ||
                (elem == "eXX") ||
                (elem == "eYY") ||
                (elem == "eZZ") ||

                (elem == "eXY") ||
                (elem == "eXZ") ||
                (elem == "eYZ") ||

                (elem == "Vp") ||
                (elem == "Vchg") ||
                (elem == "Vext")

         )  || // --- or known parameter
            (paramNames.contains(elem))
            || // --- or complex number
            (elem == "i")
            ||
            (
             (elem.substitute("i","").substitute("I","")
              ).isDouble()
             )
            || // --- or function Sqrt{}
            (
             (elem.trim().substitute("_{", "", 1)
               .substitute("}", "", 1)
              ).isDouble()
            )
            || // --- or function Sin{}
            (
             (elem.trim().substitute("!{", "", 1)
               .substitute("}", "", 1)
              ).isDouble()
            )
            || // --- or function Cos{}
            (
             (elem.trim().substitute("?{", "", 1)
               .substitute("}", "", 1)
              ).isDouble()
            )
            || // --- or function Tan{}
            (
             (elem.trim().substitute(".{", "", 1)
               .substitute("}", "", 1)
              ).isDouble()
            )
      )  {
      } else  {
         unknownElem = elem;
      }
   }   
  return unknownElem;
}

SxString SxKdotP::replaceUnknown(SxString elem, SxString ham)
{
   elem = elem.substitute(" ", "");
   SxString elemOut = elem;
   /* Blanks are required to make sure that unknown
    strings are not replaced as parts of other strings.
    Blanks have to be removed afterwards again.*/
   elemOut = elemOut.substitute("+", " + ");
   elemOut = elemOut.substitute("-", " - ");
   elemOut = elemOut.substitute("*", " * ");
   elemOut = elemOut.substitute("/", " / ");
   elemOut = elemOut.substitute("(", " ( ");
   elemOut = elemOut.substitute(")", " ) ");

   // --- split element in single contributions
   elem = elem.substitute("+", "|");
   elem = elem.substitute("-", "|");
   elem = elem.substitute("*", "|");
   elem = elem.substitute("/", "|");
   elem = elem.substitute("(", "|");
   elem = elem.substitute(")", "|");

   // --- set up list of single contributions
   SxList<SxString> split = elem.tokenize('|');
   
   SxString unknownElem = containsUnknown(split);

   if (unknownElem != "###")  { // --- unknown element detected
      SxString replacement;
      SxString search = unknownElem + "=";
      int srch = (int)ham.find(search);
      if ( (srch > -1 ) && (srch < ham.getSize()) )  {
         // --- check wether all terms A = ... contain semicolon
         replacement = "(" + ham.right(unknownElem + "=").left(";") + ")";
         if (replacement.contains("=") > 0) {
             cout << "Error: '" << unknownElem << "=' lacks semicolon."
                  << endl;
             SX_EXIT;
         }
         elemOut = elemOut.substitute(unknownElem, replacement);
      }
      else  {
      cout << "Error: element " << unknownElem
         << " is not found anywhere. EXITING."<< endl;
      cout << paramNames << endl;
      SX_EXIT;
      }
      elemOut = replaceUnknown(elemOut.trim(), ham);

   }
   return elemOut;
}

int SxKdotP::resolveExpr(SxString expr, int col, int row)
{
   int pos = -1;
   // --- check for +
   SxString left, right;
   int l = 0, r  = 0, rtn = 0;
   pos = firstOutsideBracket(expr, '+');
   if (pos > -1)  {
      left = expr.subString(0, pos - 1);
      if (left.trim() == "") left = "0";
      right = expr.subString(pos + 1);
      expression(col)(row) << "+";
      rtn = (int)expression(col)(row).getSize() - 1;
      leftPtr(col)(row).resize(rtn + 1);
      rightPtr(col)(row).resize(rtn + 1);
      l = resolveExpr(left, col, row);
      r = resolveExpr(right, col, row);
      leftPtr(col)(row)(rtn) = l;
      rightPtr(col)(row)(rtn) = r;
   } 
   // --- check for -
   if (pos == -1) {
      pos = firstOutsideBracket(expr, '-');
      if (pos > -1)  {
         left = expr.subString(0, pos - 1);
         if (left.trim() == "") left = "0";
         right = expr.subString(pos + 1);
         expression(col)(row) << "-";
         rtn = (int)expression(col)(row).getSize() - 1;
         leftPtr(col)(row).resize(rtn + 1);
         rightPtr(col)(row).resize(rtn + 1);
         l = resolveExpr(left, col, row);
         r = resolveExpr(right, col, row);
         leftPtr(col)(row)(rtn) = l;
         rightPtr(col)(row)(rtn) = r;
      } 
   }
   // --- check for *
   if (pos == -1) {
      pos = firstOutsideBracket(expr, '*');
      if (pos > -1)  {
         left = expr.subString(0, pos - 1);
         right = expr.subString(pos + 1);
         expression(col)(row) << "*";
         rtn = (int)expression(col)(row).getSize() - 1;
         leftPtr(col)(row).resize(rtn + 1);
         rightPtr(col)(row).resize(rtn + 1);
         l = resolveExpr(left, col, row);
         r = resolveExpr(right, col, row);
         leftPtr(col)(row)(rtn) = l;
         rightPtr(col)(row)(rtn) = r;
      }
   }
   // --- check for /
   if (pos == -1) {
      pos = firstOutsideBracket(expr, '/');
      if (pos > -1)  {
         left = expr.subString(0, pos - 1);
         right = expr.subString(pos + 1);
         expression(col)(row) << "/";
         rtn = (int)expression(col)(row).getSize() - 1;
         leftPtr(col)(row).resize(rtn + 1);
         rightPtr(col)(row).resize(rtn + 1);
         l = resolveExpr(left, col, row);
         r = resolveExpr(right, col, row);
         leftPtr(col)(row)(rtn) = l;
         rightPtr(col)(row)(rtn) = r;
      } 
   }
   // if none of the above is found: search for brackets
   if (pos == -1)  {
      int lBrck = -1, rBrck = -1; // left bracket, right bracket's position
      int length = (int)expr.getSize();
      int i;
      for (i = 0; i < length; i++)  {
         if ((lBrck == -1) && (expr(i) == '('))
            lBrck = i;
      }
      for (i = length-1; i > 0; i--)  {
         if ((rBrck == -1) && (expr(i) == ')'))
            rBrck = i;
      }
      if ((lBrck > -1) && (rBrck > -1))  {
         expression(col)(row) << "()";
         rtn = (int)expression(col)(row).getSize() - 1;
         leftPtr(col)(row).resize(rtn + 1);
         rightPtr(col)(row).resize(rtn + 1);

         leftPtr(col)(row)(rtn) = rtn + 1;
         rightPtr(col)(row)(rtn) = -1;
         resolveExpr(expr.subString(lBrck + 1, rBrck - 1), col, row);
         pos = 0;
      } 
   }
   // --- if still no element identified: expr is a leaf
   if (pos == -1) {
      expression(col)(row) << expr.trim();
      rtn = (int)expression(col)(row).getSize() - 1;
      leftPtr(col)(row).resize(rtn + 1);
      rightPtr(col)(row).resize(rtn + 1);
      leftPtr(col)(row)(rtn) = -1;
      rightPtr(col)(row)(rtn) = -1;
   }
   return rtn;
}

SxComplex16 SxKdotP::evaluateTreeNumber (int col, int row, int pos) const
{
   const SxString &expr = expression(col)(row)(pos);
   if (expr == "()")  {
      return evaluateTreeNumber(col, row, leftPtr(col)(row)(pos));
   } else if (expr == "*")  {
      return  evaluateTreeNumber(col, row, leftPtr (col)(row)(pos))
            * evaluateTreeNumber(col, row, rightPtr(col)(row)(pos));
   } else if (expr == "/")  {
      return  evaluateTreeNumber(col, row, leftPtr (col)(row)(pos))
            / evaluateTreeNumber(col, row, rightPtr(col)(row)(pos));
   } else if (expr == "+")  {
      return  evaluateTreeNumber(col, row, leftPtr (col)(row)(pos))
            + evaluateTreeNumber(col, row, rightPtr(col)(row)(pos));
   } else if (expr == "-")  {
      return  evaluateTreeNumber(col, row, leftPtr (col)(row)(pos))
            - evaluateTreeNumber(col, row, rightPtr(col)(row)(pos));
   }
   // --- function
   else if (expr.contains("_{") > 0)  {
      if ((expr.right("_{").left("}")).isDouble())  {
         return sqrt((expr.right("_{").left("}")).toDouble());
      }
      SX_EXIT;
   } else if (expr.contains("!{") > 0)  {
      if ((expr.right("!{").left("}")).isDouble())  {
         return sin((expr.right("!{").left("}")).toDouble());
      }
      SX_EXIT;
   } else if (expr.contains("?{") > 0)  {
      if ((expr.right("?{").left("}")).isDouble())  {
         return cos((expr.right("?{").left("}")).toDouble());
      }
      SX_EXIT;
   } else if (expr.contains(".{") > 0)  {
      if ((expr.right(".{").left("}")).isDouble())  {
         return tan((expr.right(".{").left("}")).toDouble());
      }
      SX_EXIT;
   }
   // --- complex factor
   else if ((expr.contains("i") > 0) && (expr.left("i").isDouble()))  {
      return I * expr.left("i").toDouble();
   } else if (expr == "i")  {
      return I;
   } else if (expr.isDouble()) {
      return expr.toDouble();
   }

   SX_EXIT;
   return 0.;
}

// Virtual machine for the byte code
// - stack of vals
// topmost value corresponds to right leaf of the binary tree
// secondtopmost value corresponds to left leaf of the binary tree

enum ByteCode : char {
   /// Add two topmost vals
   Plus,
   // Subtract topmost val (right) from secondtopmost val (left)
   Minus,
   /// Multiply two topmost vals
   Multiply,
   // Divide secondtopmost val (left) by topmost val (right)
   Divide,
   /// operator (char: opKey)
   Operator,
   /// Multiplication with operator
   MultiplyOp,
   /// Push strain field on stack (char: straincomp) 0-5; 6 => hydrostatic pressure 0 + 1 + 2
   Strain,
   /// Push polarization potential on stack
   PolarizationPotential,
   /// Push external potential on stack
   ExtPotential,
   /// Push charge potential on stack
   ChargePotential,
   /// Add a constant (2x double: re, im)
   PlusVal,
   /// Subtract from a constant (2x double: re, im)
   MinusVal,
   /// Multiply with a complex number (2x double: re, im)
   MultiplyVal,
   /// Set topmost to a value (2x double: re, im)
   LoadVal,
   /// Load parameter field (int: iParam)
   LoadParameter,
   /// Multiply topmost val with psi-component
   MultiplyPsi
};

// --- auxiliary routine: put arbitray type on char stack byte by byte
template <class T>
inline void push (SxStack<char> &code, const T &val)
{
   const char *bytes = reinterpret_cast<const char*>(&val);
   for (size_t i = 0; i < sizeof(T); i++)
      code << bytes[i];
}

// --- auxiliary routine: get arbitray type from char array,
//     advancing the idx by the necessary amount
//     note: idx must point to the byte BEFORE the data starts
//           (=>idx still points to the bytecode that defines the type)
template <class T>
inline T getVal (const SxArray<char> &code, int *idx)
{
   T res = *reinterpret_cast<const T*>(code.elements + *idx + 1);
   (*idx) += (int)sizeof(T);
   return res;
}

// This generates the byte code from the tree (recursively); the byte code
// is added to code (character stack)
void SxKdotP::evaluateTreeByteCode(int col, int row, int pos,
                                   SxStack<char> &code) const
{
   if (pos == 0 && expression(col)(row)(0) == "0") return;
   SxString expr = expression(col)(row)(pos);
   bool noOperator = false;
   if (expr.contains("!"))  {
      noOperator = true;
      expr = expr.left("!");
   }
   // --- shortcut: if tree evaluates to a numerical constant, put the result
   if (isRealNumber(pos, col, row))  {
      SxComplex16 val = evaluateTreeNumber (col, row, pos);
      code << LoadVal;
      push (code, val);
   }
   // --- evaluate +, -, *, /
   else if (expr == "+")  {
      if (isRealNumber(rightPtr(col)(row)(pos), col, row))  {
         SxComplex16 val = evaluateTreeNumber (col, row, rightPtr(col)(row)(pos));
         if (isRealNumber(leftPtr(col)(row)(pos), col, row))  {
            // val + val
            val += evaluateTreeNumber (col, row, leftPtr(col)(row)(pos));
            code << LoadVal;
            push (code, val);
         } else {
            // vec + val
            evaluateTreeByteCode(col, row, leftPtr(col)(row)(pos), code);
            code << PlusVal;
            push (code, val);
         }
      } else if (isRealNumber(leftPtr(col)(row)(pos), col, row))  {
         SxComplex16 val = evaluateTreeNumber (col, row, leftPtr(col)(row)(pos));
         // val + vec
         evaluateTreeByteCode(col, row, rightPtr(col)(row)(pos), code);
         if (val.abs () > 1e-14)  {
            // add if non-zero
            code << PlusVal;
            push (code, val);
         }
      } else {
         // vec + vec
         evaluateTreeByteCode(col, row, leftPtr(col)(row)(pos), code);
         evaluateTreeByteCode(col, row, rightPtr(col)(row)(pos), code);
         code << Plus;
      }
   } else if (expr == "-")  {
      if (isRealNumber(rightPtr(col)(row)(pos), col, row))  {
         // minus is performed here:
         SxComplex16 val = -evaluateTreeNumber (col, row, rightPtr(col)(row)(pos));
         if (isRealNumber(leftPtr(col)(row)(pos), col, row))  {
            // val - val (minus is done above)
            val += evaluateTreeNumber (col, row, leftPtr(col)(row)(pos));
            code << LoadVal;
            push (code, val);
         } else {
            // vec - val
            evaluateTreeByteCode(col, row, leftPtr(col)(row)(pos), code);
            code << PlusVal; // minus is done above
            push (code, val);
         }
      } else if (isRealNumber(leftPtr(col)(row)(pos), col, row))  {
         SxComplex16 val = evaluateTreeNumber (col, row, leftPtr(col)(row)(pos));
         // val - vec
         evaluateTreeByteCode(col, row, rightPtr(col)(row)(pos), code);
         code << MinusVal;
         push (code, val);
      } else {
         // vec - vec
         evaluateTreeByteCode(col, row, leftPtr(col)(row)(pos), code);
         evaluateTreeByteCode(col, row, rightPtr(col)(row)(pos), code);
         code << Minus;
      }
   } else if (expr == "*") {
      int posL = leftPtr (col)(row)(pos);
      int posR = rightPtr(col)(row)(pos);
      bool leftIsNumber  = isRealNumber(posL, col, row);
      bool rightIsNumber = isRealNumber(posR, col, row);
      bool leftBracketOp = expression(col)(row)(posL) == "()"
                         && containsOperator(posL, col, row);
      bool rightBracketOp = expression(col)(row)(posR) == "()"
                         && containsOperator(posR, col, row);
      if (   ( leftBracketOp && ! rightIsNumber)
          || (rightBracketOp && !  leftIsNumber))
      {
         cout << "Error in element: " << col << ", " << row << ": ";
         printElement(pos,col,row);
         cout << endl;
         cout << "Multiplication with bracket that contains operator "
              << "detected.\nThis construct is not allowed in "
              << "byteCode mode. EXITING." << endl;
         SX_EXIT;
      }

      if (expression(col)(row)(posR).contains("k"))
      {
         // swap order: expr * op = op * expr
         //cout << "Assuming operator commutes in row=" << row
         //     << " col=" << col << endl;
         int oldL = posL;
         posL = posR;
         posR = oldL;
      }
      if (isRealNumber(posR, col, row))  {
         SxComplex16 val = evaluateTreeNumber (col, row, posR);
         if (isRealNumber(posL, col, row))  {
            // val * val
            val *= evaluateTreeNumber (col, row, posL);
            code << LoadVal;
            push (code, val);
         } else {
            // vec * val
            evaluateTreeByteCode(col, row, posL, code);
            code << MultiplyVal;
            push (code, val);
         }
      } else if (isRealNumber(posL, col, row))  {
         SxComplex16 val = evaluateTreeNumber (col, row, posL);
         // val * vec
         evaluateTreeByteCode(col, row, posR, code);
         code << MultiplyVal;
         push (code, val);
      } else {
         // vec * vec
         evaluateTreeByteCode(col, row, posL, code);
         evaluateTreeByteCode(col, row, posR, code);
         // dangerous: should be replaced by proper check
         // simple parameters could contain the letter 'k'
         if (expression(col)(row)(posL).contains ("k"))
            code << MultiplyOp;
         else
            code << Multiply;
      }
   } else if (expr == "/") {
      if (isRealNumber(rightPtr(col)(row)(pos), col, row))  {
         // division is done here:
         SxComplex16 val = 1. / evaluateTreeNumber (col, row, rightPtr(col)(row)(pos));
         if (isRealNumber(leftPtr(col)(row)(pos), col, row))  {
            // val / val
            val *= evaluateTreeNumber (col, row, leftPtr(col)(row)(pos));
            code << LoadVal;
            push (code, val);
         } else {
            // vec / val
            evaluateTreeByteCode(col, row, leftPtr(col)(row)(pos), code);
            code << MultiplyVal;
            push (code, val);
         }
      } else if (isRealNumber(leftPtr(col)(row)(pos), col, row))  {
         // val / vec
         SxComplex16 val = evaluateTreeNumber (col, row, leftPtr(col)(row)(pos));
         code << LoadVal;
         push (code, val);
         evaluateTreeByteCode(col, row, rightPtr(col)(row)(pos), code);
         code << Divide;
      } else {
         // vec / vec
         evaluateTreeByteCode(col, row, leftPtr(col)(row)(pos), code);
         evaluateTreeByteCode(col, row, rightPtr(col)(row)(pos), code);
         code << Divide;
      }
   } else if (expr == "()")  {
      // --- evaluate ()
      evaluateTreeByteCode(col, row, leftPtr(col)(row)(pos), code);
   }

   // --- evaluate leafs

   // --- k^2 - operators, also mixed op's
   else if (expr.contains ("k"))  { // dangerous, see above
      int op = resolveOpKey (expr);
      code << Operator << char(op);
   }

   // --- strains
   else if (expr == "e") code << Strain << 6;
   else if (expr == "eXX") code << Strain << 0;
   else if (expr == "eYY") code << Strain << 1;
   else if (expr == "eZZ") code << Strain << 2;
   else if (expr == "eXY") code << Strain << 3;
   else if (expr == "eXZ") code << Strain << 4;
   else if (expr == "eYZ") code << Strain << 5;

   // --- external or polarization potential
   else if (expr == "Vp") code << PolarizationPotential;
   else if (expr == "Vext") code << ExtPotential;
   else if (expr == "Vchg") code << ChargePotential;

   // --- function
   else if (expr.contains("_{") > 0)  {
      if ((expr.right("_{").left("}")).isDouble())  {
         double val = sqrt((expr.right("_{").left("}")).toDouble());
         code << LoadVal;
         push (code, val);
         push (code, double(0.));
      } else {
         SX_EXIT;
      }
   }
   else if (expr.contains("!{") > 0)  {
      if ((expr.right("!{").left("}")).isDouble())  {
         double val = sin((expr.right("!{").left("}")).toDouble());
         code << LoadVal;
         push (code, val);
         push (code, double(0.));
      } else {
         SX_EXIT;
      }
   }
   else if (expr.contains("?{") > 0)  {
      if ((expr.right("?{").left("}")).isDouble())  {
         double val = cos((expr.right("?{").left("}")).toDouble());
         code << LoadVal;
         push (code, val);
         push (code, double(0.));
      } else {
         SX_EXIT;
      }
   }
   else if (expr.contains(".{") > 0)  {
      if ((expr.right(".{").left("}")).isDouble())  {
         double val = tan((expr.right(".{").left("}")).toDouble());
         code << LoadVal;
         push (code, val);
         push (code, double(0.));
      } else {
         SX_EXIT;
      }
   }
   // --- complex factor
   else if ((expr.contains("i") > 0) && (expr.left("i").isDouble()))  {
      double val = expr.left("i").toDouble();
      code << LoadVal;
      push (code, double(0.));
      push (code, val);
   }
   else if (expr == "i")  {
      code << LoadVal;
      push (code, double(0.));
      push (code, double(1.));
   }
   else if (expr.isDouble())  {
      double val = expr.toDouble();
      code << LoadVal;
      push (code, val);
      push (code, double(0.));
   }

   // --- known parameter
   else  {
      int pIdx;
      for (pIdx = 0; pIdx < paramNames.getSize(); pIdx++)
         if (expr == paramNames(pIdx)) {
            code << LoadParameter;
            push (code, pIdx);
         }
   }
   if (noOperator)
      code << MultiplyPsi;
}

SxArray2<SxArray<char> > SxKdotP::generateByteCode () const
{
   SxArray2<SxArray<char> > res(nComp, nComp);
   for (int col = 0; col < nComp; col++)  {
      for (int row = 0; row < nComp; row++)  {
         SxStack<char> code;
         evaluateTreeByteCode (col, row, 0, code);
         if (!code.isEmpty ())
            res(col, row) = code;
         //cout << "Evaluation program for col=" << col << " row=" << row << endl;
         //printByteCode (res(col, row));
      }
   }
   return res;
}

// --- auxiliary routine: find step to next bytecode command, keeping
//     track of the stack size
//     This routine allows to scan through bytecode without doing anything.
//     return value: number of bytes for the next byte code
//                   this is larger than 1 if the byte code has parameters
inline int nextCode (char code, int *stackSize)
{
   SX_CHECK (stackSize);
   switch (code)  {
      case Plus:
      case Minus:
      case Multiply:
      case Divide:
         (*stackSize)--; return 1;

      case Operator:
         (*stackSize)++; return 2;
      case MultiplyOp:
         (*stackSize)--; return 1;
      case Strain:
         (*stackSize)++; return 2;

      case PolarizationPotential:
      case ExtPotential:
      case ChargePotential:
         (*stackSize)++; return 1;
      case PlusVal:
      case MinusVal:
      case MultiplyVal:
         return 1 + int(sizeof(SxComplex16));
      case LoadVal:
         (*stackSize)++;
         return 1 + int(sizeof(SxComplex16));
      case LoadParameter:
         (*stackSize)++; return 1 + sizeof(int);
      case MultiplyPsi:
         return 1;
      default:
         SX_EXIT;
   }
   return 1;
}


void SxKdotP::printByteCode (const SxArray<char> &code) const
{
   if (code.getSize () == 0) { cout << "No op." << endl; return; }
   int stackSize = 0;
   for (int i = 0; i < code.getSize (); i++)  {
      switch (code(i))  {
         case Plus:
            cout << "Plus (" << --stackSize << ")" << endl; break;
         case Minus:
            cout << "Minus (" << --stackSize << ")" << endl; break;
         case Multiply:
            cout << "Multiply (" << --stackSize << ")" << endl; break;
         case Divide:
            cout << "Divide (" << --stackSize << ")" << endl; break;
         case Operator:
            {
               char opCode = code(++i);
               switch (opCode)  {
                  case 0:  cout << "k"; break;
                  case 1:  cout << "kx"; break;
                  case 2:  cout << "ky"; break;
                  case 3:  cout << "kz"; break;
                  case 4:  cout << "k2"; break;
                  case 5:  cout << "kx2"; break;
                  case 6:  cout << "ky2"; break;
                  case 7:  cout << "kz2"; break;
                  case 8:  cout << "kxy"; break;
                  case 9:  cout << "kxz"; break;
                  case 10: cout << "kyz"; break;
                  default: SX_EXIT;
               }
               cout << " (" << ++stackSize << ")" << endl;
               break;
            }
         case MultiplyOp:
            cout << "Multiply with operator (" << --stackSize << ")" << endl;
            break;
         case Strain:
            switch (code(++i))  {
               case 0: cout << "strain xx"; break;
               case 1: cout << "strain yy"; break;
               case 2: cout << "strain zz"; break;
               case 3: cout << "strain xy"; break;
               case 4: cout << "strain xz"; break;
               case 5: cout << "strain yz"; break;
               case 6: cout << "pressure"; break;
            }
            cout << " (" << ++stackSize << ")" << endl;
            break;
         case PolarizationPotential:
            cout << "Vp (" << ++stackSize << ")" << endl; break;
         case ExtPotential:
            cout << "Vext (" << ++stackSize << ")" << endl; break;
         case ChargePotential:
            cout << "Vchg (" << ++stackSize << ")" << endl; break;
         case PlusVal:
            cout << "Add " << getVal<SxComplex16> (code, &i) << ". ("
                 << stackSize << ")" << endl;
            break;
         case MinusVal:
            cout << "Subtract from " << getVal<SxComplex16> (code, &i) << ". ("
                 << stackSize << ")" << endl;
            break;
         case MultiplyVal:
            cout << "Multiply with " << getVal<SxComplex16> (code, &i) << ". ("
                 << stackSize << ")" << endl;
            break;
         case LoadVal:
            cout << "Load " << getVal<SxComplex16> (code, &i) << ". ("
                 << ++stackSize << ")" << endl;
            break;
         case LoadParameter:
            cout << "Parameter '" << paramNames(getVal<int> (code, &i))
                 << "' (" << ++stackSize << ")" << endl;
            break;
         case MultiplyPsi:
            cout << "Multiply with psi (" << stackSize << ")" << endl;
            break;
         default:
            SX_EXIT;
      }
   }
   if (stackSize != 1)  {
      cout << "Broken code: final stacksize is not 1" << endl;
      SX_EXIT;
   }
}

// --- auxiliary routine
// create a changeable temporary reference from a const PsiRef
inline PsiRef newRef (const PsiRef &in)
{
   return const_cast<PsiRef&> (in);
}

PsiRef SxKdotP::evaluateByteCode (const SxArray<char> &code)
{
   if (code.getSize () == 0) {
      return const_cast<PsiR&>(zero);
   }
   SxStack<PsiRef> stack;
   for (int i = 0; i < code.getSize (); i++)  {
      switch (code(i))  {
         case Plus:
            {
               PsiRef a = stack.pop (), b = stack.pop ();
               stack.push (a + b);
               break;
            }
         case Minus:
            {
               PsiRef a = stack.pop (), b = stack.pop ();
               stack.push (b - a);
               break;
            }
         case Multiply:
            {
               PsiRef a = stack.pop (), b = stack.pop ();
               stack.push (a * b);
               break;
            }
         case Divide:
            {
               PsiRef a = stack.pop (), b = stack.pop ();
               stack.push (b / a);
               break;
            }
         case Operator:
            {
               int opCode = code(++i);
               int stackSize = 1;
               int j;
               for (j = i + 1; j < code.getSize (); )  {
                  if (code(j) == MultiplyOp) break;
                  if (stackSize == 0) break;
                  j += nextCode(code(j), &stackSize);
               }
               if (j < code.getSize () && code(j) == MultiplyOp)
               {
                  if (par(iMult).getSize () > 0)  {
                     // short-cut: use precomputed par(iMult)
                     stack.push (returnAccurate(iMult));
                     i = j;
                     iMult++;
                  } else {
                     opKey(iMult) = opCode;
                     // and now continue until we arrive at MultiplyOp below...
                  }
               } else {
                  // --- apply operator
                  switch (opCode)  {
                     case 0:  stack << (fdC(0) + fdC(1) + fdC(2)); break; // k
                     case 1:  stack << newRef(fdC(0)); break;   // kx
                     case 2:  stack << newRef(fdC(1)); break;   // ky
                     case 3:  stack << newRef(fdC(2)); break;   // kz
                     case 4:  stack << (sdC(0,0) + sdC(1,1) + sdC(2,2)); // k2
                              break;
                     case 5:  stack << newRef(sdC(0,0)); break; // kx2
                     case 6:  stack << newRef(sdC(1,1)); break; // ky2
                     case 7:  stack << newRef(sdC(2,2)); break; // kz2
                     case 8:  stack << newRef(sdC(0,1)); break; // kxy
                     case 9:  stack << newRef(sdC(0,2)); break; // kxz
                     case 10: stack << newRef(sdC(1,2)); break; // kyz
                     default: SX_EXIT;
                  }
               }
               break;
            }
         case MultiplyOp:
            // cache parameter factor
            par(iMult) = stack.pop ();
            // compute the necessary derivatives
            determineDerivatives (iMult);
            // and calculate the accurate product
            // 0.5 * [ op | (par * psi) + par * (op | psi) ]
            stack.push (returnAccurate(iMult));
            iMult++;
            break;
         case Strain:
            {
               int idx = code (++i);
               if (idx == 6)
                  stack << (eIJ(0) + eIJ(1) + eIJ(2)); // pressure
               else
                  stack << eIJ(idx);
               break;
            }
         case PolarizationPotential:
            stack << vP; break;
         case ExtPotential:
            stack << vExt; break;
         case ChargePotential:
            stack << chargePotential; break;
         case PlusVal:
            stack.top () += getVal<SxComplex16> (code, &i);
            break;
         case MinusVal:
            stack << (getVal<SxComplex16> (code, &i) - stack.pop ());
            break;
         case MultiplyVal:
            stack.top () *= getVal<SxComplex16> (code, &i);
            break;
         case LoadVal:
            {
               PsiR valVec (*rBasisPtr);
               valVec.set ( getVal<SxComplex16> (code, &i) );
               stack << (PsiRef&&)(valVec);
               break;
            }
         case LoadParameter:
            stack << parameters(getVal<int> (code, &i) );
            break;
         case MultiplyPsi:
            stack.top () *= psiR;
            break;
         default:
            SX_EXIT;
      }
   }
   SX_CHECK (stack.getSize () == 1, stack.getSize ());
   return stack.pop ();
}

SxComplex16 SxKdotP::evaluateByteCode (const SxArray<char> &code,
                                       const SxVector3<double> &kVec,
                                       const SxVector<double> &strain,
                                       double vPol, double vExtern, double vCharge,
                                       const SxVecRef<double> &params,
                                       SxStack<SxComplex16> &stack) const
{
   SX_CHECK (strain.getSize () == 6, strain.getSize ());
   SX_CHECK (params.getSize () == nParam, params.getSize (), nParam);
   if (code.getSize () == 0) {
      return 0.;
   }
   // This routine seems to be really fast. It may be necessary to remove the timer
   // to ensure that the timer does not become the bottleneck.
   //SX_CLOCK(Timer::EvalByteCodeBS);
   for (int i = 0; i < code.getSize (); i++)  {
      switch (code(i))  {
         case Plus:
            {
               SxComplex16 a = stack.pop (), b = stack.pop ();
               stack.push (a + b);
               break;
            }
         case Minus:
            {
               SxComplex16 a = stack.pop (), b = stack.pop ();
               stack.push (b - a);
               break;
            }
         case Multiply:
         case MultiplyOp: // same as normal multiply
            {
               SxComplex16 a = stack.pop (), b = stack.pop ();
               stack.push (a * b);
               break;
            }
         case Divide:
            {
               SxComplex16 a = stack.pop (), b = stack.pop ();
               stack.push (b / a);
               SX_CHECK_NUM (stack.top ());
               break;
            }
         case Operator:
            // --- apply operator
            switch (code(++i))  {
               case 0:  stack << kVec.sum (); break; // k
               case 1:  stack << kVec(0); break;   // kx
               case 2:  stack << kVec(1); break;   // ky
               case 3:  stack << kVec(2); break;   // kz
               case 4:  stack << kVec.normSqr ();
                        break;
               case 5:  stack << kVec(0) * kVec(0); break; // kx2
               case 6:  stack << kVec(1) * kVec(1); break; // ky2
               case 7:  stack << kVec(2) * kVec(2); break; // kz2
               case 8:  stack << kVec(0) * kVec(1); break; // kxy
               case 9:  stack << kVec(0) * kVec(2); break; // kxz
               case 10: stack << kVec(1) * kVec(2); break; // kyz
               default: SX_EXIT;
            }
            break;
         case Strain:
            {
               int idx = code (++i);
               if (idx == 6)
                  stack << (strain(0) + strain(1) + strain(2)); // pressure
               else
                  stack << strain(idx);
               SX_CHECK_NUM (stack.top ());
               break;
            }
         case PolarizationPotential:
            stack << vPol;
            SX_CHECK_NUM (stack.top ());
            break;
         case ExtPotential:
            stack << vExtern;
            SX_CHECK_NUM (stack.top ());
            break;
         case ChargePotential:
            stack << vCharge;
            SX_CHECK_NUM (stack.top ());
            break;
         case PlusVal:
            stack.top () += getVal<SxComplex16> (code, &i);
            break;
         case MinusVal:
            stack << (getVal<SxComplex16> (code, &i) - stack.pop ());
            break;
         case MultiplyVal:
            stack.top () *= getVal<SxComplex16> (code, &i);
            break;
         case LoadVal:
            stack << getVal<SxComplex16> (code, &i);
            break;
         case LoadParameter:
            stack << params(getVal<int> (code, &i) );
            SX_CHECK_NUM (stack.top ());
            break;
         case MultiplyPsi:
            // no op
            break;
         default:
            SX_EXIT;
      }
   }
   SX_CHECK (stack.getSize () == 1, stack.getSize ());
   return stack.pop ();
}

int SxKdotP::resolveOpKey(const SxString &opString)
{
   if (opString == "k") return 0;
   else if (opString == "kx") return 1;
   else if (opString == "ky") return 2;
   else if (opString == "kz") return 3;
   else if (opString == "k2") return 4;
   else if (opString == "kx2") return 5;
   else if (opString == "ky2") return 6;
   else if (opString == "kz2") return 7;
   else if (opString == "kxy") return 8;
   else if (opString == "kxz") return 9;
   else if (opString == "kyz") return 10;
   else return -1;
}

PsiR SxKdotP::returnAccurate(int i)
{
   PsiR rtn;
   SX_CLOCK(Timer::kk2s);
   switch (opKey(i)) {
      case 0: // k = kx + ky + kz
              rtn = par(i) * (fdC(0) + fdC(1) + fdC(2));
              rtn.plus_assign_ax (0.5, kIPar(i) * psiR);
              break;
      case 1: // kx
              rtn = par(i) * fdC(0);
              rtn.plus_assign_ax (0.5, kIPar(i) * psiR);
              break;
      case 2: // ky
              rtn = par(i) * fdC(1);
              rtn.plus_assign_ax (0.5, kIPar(i) * psiR);
              break;
      case 3: // kz
              rtn = par(i) * fdC(2);
              rtn.plus_assign_ax (0.5, kIPar(i) * psiR);
              break;
      case 4: // k2 = kx2 + ky2 + kz2
              rtn = par(i) * (sdC(0,0) + sdC(1,1) + sdC(2,2))
                  + kIPar(i) * fdC(0)
                  + kJPar(i) * fdC(1)
                  + kKPar(i) * fdC(2);
              break;
      case 5: // kx2
              rtn = par(i) * sdC(0,0) + kIPar(i) * fdC(0);
              break;
      case 6: // ky2
              rtn = par(i) * sdC(1,1) + kIPar(i) * fdC(1);
              break;
      case 7: // kz2
              rtn = par(i) * sdC(2,2) + kIPar(i) * fdC(2);
              break;
      case 8: // kxy = kx*ky
              rtn = par(i) * sdC(0,1);
              rtn.plus_assign_ax (0.5, kIPar(i) * fdC(1) + kJPar(i) * fdC(0));
              break;
      case 9: // kxz = kx * kz
              rtn = par(i) * sdC(0,2);
              rtn.plus_assign_ax (0.5, kIPar(i) * fdC(2) + kJPar(i) * fdC(0));
              break;
      case 10: // kyz = ky * kz
              rtn = par(i) * sdC(1,2);
              rtn.plus_assign_ax (0.5, kIPar(i) * fdC(2) + kJPar(i) * fdC(1));
              break;
      default:
              SX_EXIT;
   }
   return rtn;
}

void SxKdotP::determineDerivatives(int i)
{
   const SxGBasis &G = *gBasisPtr;
   const SxRBasis &R = *rBasisPtr;
   PsiG parG = G | par(i);
   // ensure that par(i) is Fourier filtered
   // this is relevant for keeping the Hamiltonian Hermitean
   // but not sufficient if mesh is too small
   par(i) = R | parG;

   // first derivatives of par
   SxArray<PsiR> dPar(3);
   // -I d/dr = -I * (-I G) = -G
   parG *= -1; // common prefactor
   for (int iDir = 0; iDir < 3; iDir++) {
      dPar(iDir) = R | (G.gVec.colRef(iDir) * parG);
      dPar(iDir) *= wgt(iDir);
   }

   if (opKey(i) == 0)  {// k
      kIPar(i) = dPar(0) + dPar(1) + dPar(2);
      kJPar(i) = zero;
   }
   if (opKey(i) == 1)  {// kx
      kIPar(i) = dPar(0);
      kJPar(i) = zero;
   }
   if (opKey(i) == 2)  {// ky
      kIPar(i) = dPar(1);
      kJPar(i) = zero;
   }
   if (opKey(i) == 3)  {// kz
      kIPar(i) = dPar(2);
      kJPar(i) = zero;
   }
   if (opKey(i) == 4)  {// k2
      kIPar(i) = dPar(0);
      kJPar(i) = dPar(1);
      kKPar(i) = dPar(2);
   }
   if (opKey(i) == 5)  {// kx2
      kIPar(i) = dPar(0);
      kJPar(i) = kIPar(i);
   }
   if (opKey(i) == 6)  {// ky2
      kIPar(i) = dPar(1);
      kJPar(i) = kIPar(i);
   }
   if (opKey(i) == 7)  {// kz2
      kIPar(i) = dPar(2);
      kJPar(i) = kIPar(i);
   }
   if (opKey(i) == 8)  {// kxy
      kIPar(i) = dPar(0);
      kJPar(i) = dPar(1);
   }
   if (opKey(i) == 9)  {// kxz
      kIPar(i) = dPar(0);
      kJPar(i) = dPar(2);
   }
   if (opKey(i) == 10)  {// kyz
      kIPar(i) = dPar(1);
      kJPar(i) = dPar(2);
   }
}

void SxKdotP::printElement(int pos, int col, int row) const
{
   SxString expr = expression(col)(row)(pos);
   if (expr.contains ("!"))  {
      cout << '!';
      expr = expr.left ("!");
   }
   if (expr == "()")  {
      cout << "(";
      printElement(leftPtr(col)(row)(pos), col, row);
      cout << ")";
   }
   else if (expr == "*")  {
      printElement(leftPtr(col)(row)(pos), col, row);
      cout << " * ";
      printElement(rightPtr(col)(row)(pos), col, row);
   }
   else if (expr == "/")  {
      printElement(leftPtr(col)(row)(pos), col, row);
      cout << " / ";
      printElement(rightPtr(col)(row)(pos), col, row);
   }
   else if (expr == "+")  {
      printElement(leftPtr(col)(row)(pos), col, row);
      cout << " + ";
      printElement(rightPtr(col)(row)(pos), col, row);
   }
   else if (expr == "-")  {
      printElement(leftPtr(col)(row)(pos), col, row);
      cout << " - ";
      printElement(rightPtr(col)(row)(pos), col, row);
   } else {
     cout << expr;
   }

}

// --- if used for preconditioner
PsiRef SxKdotP::evaluateTreePrecond(int col, int row, int pos)
{
   PsiRef rtn;
   if (col != row) {
      // zero preconditioner element in reciprocal space
      return gZero;
   }
   int pIdx;
   // --- in case no operator present in a summand, multiply with psiR
   bool noOperator = false;
   SxString expr = expression(col)(row)(pos);
   if (expr.contains("!"))  {
      noOperator = true;
      expr = expr.left("!");
   }
   // --- evaluate +, -, *, /
   if (expr == "+")
      rtn = evaluateTreePrecond(col, row, leftPtr(col)(row)(pos))
         + evaluateTreePrecond(col, row, rightPtr(col)(row)(pos));
   else if (expr == "-")
      rtn = evaluateTreePrecond(col, row, leftPtr(col)(row)(pos))
         - evaluateTreePrecond(col, row, rightPtr(col)(row)(pos));
   else if (expr == "*")
      rtn = evaluateTreePrecond(col, row, leftPtr(col)(row)(pos))
         * evaluateTreePrecond(col, row, rightPtr(col)(row)(pos));
   else if (expr == "/")
      rtn = evaluateTreePrecond(col, row, leftPtr(col)(row)(pos))
         / evaluateTreePrecond(col, row, rightPtr(col)(row)(pos));
   // --- evaluate ()
   else if (expr == "()") {
        rtn = evaluateTreePrecond(col, row, leftPtr(col)(row)(pos));
      }
   // --- evaluate leafs

   // --- k^2 - operators, also mixed op's
   else if (expr == "k2") rtn = sdC(0,0) + sdC(1,1) + sdC(2,2) ;
   else if (expr == "kx2") rtn = sdC(0,0);
   else if (expr == "ky2") rtn = sdC(1,1);
   else if (expr == "kz2") rtn = sdC(2,2);
   else if (expr == "kxy") rtn = sdC(0,1);
   else if (expr == "kxz") rtn = sdC(0,2);
   else if (expr == "kyz") rtn = sdC(1,2);

   // --- linear k - operators

   else if (expr == "k") rtn = fdC(0) + fdC(1) + fdC(2);
   else if (expr == "kx") rtn = fdC(0);
   else if (expr == "ky") rtn = fdC(1);
   else if (expr == "kz") rtn = fdC(2);

   // --- strains - no contribution to preconditioner

   else if (expr == "e") rtn = gZero;
   else if (expr == "eXX") rtn = gZero;
   else if (expr == "eYY") rtn = gZero;
   else if (expr == "eZZ") rtn = gZero;
   else if (expr == "eXY") rtn = gZero;
   else if (expr == "eXZ") rtn = gZero;
   else if (expr == "eYZ") rtn = gZero;

   // --- additional potentials - no contribution to preconditioner
   else if (expr == "Vp") rtn = gZero;
   else if (expr == "Vext") rtn = gZero;
   else if (expr == "Vchg") rtn = gZero;

   // --- complex factor
   else if (expr.contains("_{") > 0)  {
      if ((expr.right("_{").left("}")).isDouble())  {
         rtn = sqrt((expr.right("_{").left("}")).toDouble()) * gOne;
      }
   }
   else if (expr.contains("!{") > 0)  {
      if ((expr.right("!{").left("}")).isDouble())  {
         rtn = sin((expr.right("!{").left("}")).toDouble()) * gOne;
      }
   }
   else if (expr.contains("?{") > 0)  {
      if ((expr.right("?{").left("}")).isDouble())  {
         rtn = cos((expr.right("?{").left("}")).toDouble()) * gOne;
      }
   }
   else if (expr.contains(".{") > 0)  {
      if ((expr.right(".{").left("}")).isDouble())  {
         rtn = tan((expr.right(".{").left("}")).toDouble()) * gOne;
      }
   }
   else if ((expr.contains("i") > 0) && (expr.left("i").isDouble()))  {
      rtn = I * expr.left("i").toDouble() * gOne;
   }
   else if (expr == "i")
      rtn = I * gOne;
   else if (expr.isDouble()) rtn = expr.toDouble() * gOne;

   // --- known parameter
   else
      for (pIdx = 0; pIdx < paramNames.getSize(); pIdx++) {
         if (expr == paramNames(pIdx))
            rtn = matParam(precMaterial,pIdx) * gOne;
      }

   // --- in preconditioner: if no operator is present, element = 0
   if (noOperator)  {
      rtn.unref ();
      rtn = gZero;
   }
   if (rtn.getSize () == 0)
      rtn = gZero;

   SX_CHECK (rtn.getSize () == gZero.getSize (),
             rtn.getSize (), gZero.getSize ());
   return rtn;
}

bool SxKdotP::isRealNumber(int pos, int col, int row) const
{
bool rtn = false;
   SxString expr = expression(col)(row)(pos);
   if (expr == "()")  {
      rtn = isRealNumber(leftPtr(col)(row)(pos), col, row);
   }
   else if (expr == "*")  {
      rtn = isRealNumber(leftPtr(col)(row)(pos), col, row) 
      && isRealNumber(rightPtr(col)(row)(pos), col, row);
   }
   else if (expr == "/")  {
      rtn = isRealNumber(leftPtr(col)(row)(pos), col, row) 
      && isRealNumber(rightPtr(col)(row)(pos), col, row);
   }
   else if (expr == "+")  {
      rtn = isRealNumber(leftPtr(col)(row)(pos), col, row) 
      && isRealNumber(rightPtr(col)(row)(pos), col, row);
   }
   else if (expr == "-")  {
      rtn = isRealNumber(leftPtr(col)(row)(pos), col, row) 
      && isRealNumber(rightPtr(col)(row)(pos), col, row);
   } else if (
     ((expr.substitute("i","")).substitute("_{",""))
      .substitute("}","").isDouble())  {
      rtn = true;
   }

return rtn;
}

bool SxKdotP::containsOperator(int pos, int col, int row) const
{
   bool rtn = false;
   SxString expr = expression(col)(row)(pos);
//   cout << "expr: " << expr << endl;
   if (expr == "()")  {
      rtn = containsOperator(leftPtr(col)(row)(pos), col, row);
   }
   if (expr == "*")  {
      rtn = containsOperator(leftPtr(col)(row)(pos), col, row) 
      || containsOperator(rightPtr(col)(row)(pos), col, row);
   }
   if (expr == "/")  {
      rtn = containsOperator(leftPtr(col)(row)(pos), col, row) 
      || containsOperator(rightPtr(col)(row)(pos), col, row);
   }
   if (expr == "+")  {
      rtn = containsOperator(leftPtr(col)(row)(pos), col, row) 
      || containsOperator(rightPtr(col)(row)(pos), col, row);
   }
   if (expr == "-")  {
      rtn = containsOperator(leftPtr(col)(row)(pos), col, row) 
      || containsOperator(rightPtr(col)(row)(pos), col, row);
   }
   if ((expr == "k")
   || (expr == "kx")
   || (expr == "ky")
   || (expr == "kz")
   || (expr == "k2")
   || (expr == "kx2")
   || (expr == "ky2")
   || (expr == "kz2")
   || (expr == "kxy")
   || (expr == "kxz")
   || (expr == "kyz"))  {
      rtn = true;
   };
   return rtn;
}

int SxKdotP::whatIsElement(int pos, int col, int row)
{
   int out = -1; // 0: operator, 1: potential, 2: not clear

   /* bracket:
      If element in bracket is operator, mark bracket with "?"
      and remove "?" from element in bracket to prevent double
      operator countings.
      If element in bracket is a potential, mark bracket with "!"
      and remove "!" from element in bracket to prevent double
      multiplication with psi.
    */
   if (expression(col)(row)(pos) == "()")  {
      out = whatIsElement(leftPtr(col)(row)(pos), col, row);
      if ((out == 0))  {
         expression(col)(row)(pos) += "?";
         expression(col)(row)(leftPtr(col)(row)(pos)) =
         expression(col)(row)(leftPtr(col)(row)(pos))
            .substitute("?", "", 1);
/*    // DEBUG 11.5.2015
         expression(col)(row)(rightPtr(col)(row)(pos)) =
         expression(col)(row)(rightPtr(col)(row)(pos))
            .substitute("?", "", 1);*/
         out = 0;
      }
      else  {
         if (expression(col)(row)(leftPtr(col)(row)(pos)).contains("!"))  {
            expression(col)(row)(pos) += "!";
            expression(col)(row)(leftPtr(col)(row)(pos)) =
               expression(col)(row)(leftPtr(col)(row)(pos))
               .substitute ("!", "", 1);
            out = 1;
         }
         else  {
            out = 0;
         }
      }
   }
   
   /* multiplication, division:
      If one of the elements is an operator, mark expression with "?"
      and remove "?" and "!" from left and right branches.
      If left and right branches are potentials, mark expression with "!"
      and remove "!" from left and right branches.
      */
   else if ((expression(col)(row)(pos) == "*")
         || (expression(col)(row)(pos) == "/"))  {
      // --- multiplication yields operator
      if ((whatIsElement(leftPtr(col)(row)(pos), col, row) == 0)
         || (whatIsElement(rightPtr(col)(row)(pos), col, row) == 0))  {
         out = 0;
         expression(col)(row)(pos) += "?";
         expression(col)(row)(pos) = expression(col)(row)(pos)
            .substitute("!", "");
         expression(col)(row)(leftPtr(col)(row)(pos)) = 
         (expression(col)(row)(leftPtr(col)(row)(pos))
            .substitute ("?", "", 1)).substitute("!", "", 1);
         expression(col)(row)(rightPtr(col)(row)(pos)) = 
         (expression(col)(row)(rightPtr(col)(row)(pos))
            .substitute ("?", "", 1)).substitute("!", "", 1);
      }
      // --- multiplication yields potential
      else if ((whatIsElement(leftPtr(col)(row)(pos), col, row) == 1)
         || (whatIsElement(rightPtr(col)(row)(pos), col, row) == 1))  {
         out = 2;
         if((whatIsElement(leftPtr(col)(row)(pos), col, row) == 1)
               && (whatIsElement(rightPtr(col)(row)(pos), col, row) == 1))  {
            out = 1;
            expression(col)(row)(pos) += "!";
         }
         expression(col)(row)(leftPtr(col)(row)(pos)) = 
         expression(col)(row)(leftPtr(col)(row)(pos))
          .substitute ("!", "", 1);
         expression(col)(row)(rightPtr(col)(row)(pos)) = 
         expression(col)(row)(rightPtr(col)(row)(pos))
          .substitute ("!", "", 1);
      }
      else out = 2;
   }
   /*
      addition, subtraction:
      If left and right branches are operators, mark expression with "?"
      and remove "?" from left and right branches.
      If left and right branches are potentials, mark expression with "!"
      and remove "!" from left and right branches.
      If left and right branches are not of the same kind, do nothing
      and leave left and right branch marks unaffected.
      */
   else if ((expression(col)(row)(pos) == "+")
         || (expression(col)(row)(pos) == "-"))  {
      if ((whatIsElement(leftPtr(col)(row)(pos), col, row) == 0)
         && (whatIsElement(rightPtr(col)(row)(pos), col, row) == 0))  {
         out = 0;
         expression(col)(row)(pos) += "?";
         expression(col)(row)(leftPtr(col)(row)(pos)) = 
         expression(col)(row)(leftPtr(col)(row)(pos))
         .substitute ("?", "", 1);
         expression(col)(row)(rightPtr(col)(row)(pos)) = 
         expression(col)(row)(rightPtr(col)(row)(pos))
         .substitute ("?", "", 1);
      }
      else if ((whatIsElement(leftPtr(col)(row)(pos), col, row) == 1)
         && (whatIsElement(rightPtr(col)(row)(pos), col, row) == 1)) {
         if (
               (expression(col)(row)(leftPtr(col)(row)(pos)).contains("!"))
               && 
               (expression(col)(row)(rightPtr(col)(row)(pos)).contains("!"))
            ) {
            out = 1;
            expression(col)(row)(pos) += "!";
            expression(col)(row)(leftPtr(col)(row)(pos)) = 
               expression(col)(row)(leftPtr(col)(row)(pos))
               .substitute ("!", "", 1);
            expression(col)(row)(rightPtr(col)(row)(pos)) = 
               expression(col)(row)(rightPtr(col)(row)(pos))
               .substitute ("!", "", 1);
      }}
      else {
         out = 2;
         expression(col)(row)(pos) = expression(col)(row)(pos)
            .substitute("!", "", 1)
            .substitute("?", "", 1);
      }

   }
   else if (
         (expression(col)(row)(pos) == "k2") ||
         (expression(col)(row)(pos) == "kx2") ||
         (expression(col)(row)(pos) == "ky2") ||
         (expression(col)(row)(pos) == "kz2") ||
         (expression(col)(row)(pos) == "kxy") ||
         (expression(col)(row)(pos) == "kxz") ||
         (expression(col)(row)(pos) == "kyz") ||
         (expression(col)(row)(pos) == "k") ||
         (expression(col)(row)(pos) == "kx") ||
         (expression(col)(row)(pos) == "ky") ||
         (expression(col)(row)(pos) == "kz") 
            )  {
      out = 0;
      expression(col)(row)(pos) += "?";

   }
      else (out = 1);
   if ((out == 1) && (!expression(col)(row)(pos).contains("!"))
         && (!expression(col)(row)(pos).contains("?"))
         && (pos > 0))  {
      expression(col)(row)(pos) += "!";

   }

   if (expression(col)(row)(pos).contains("?")) out = 0;
   if (expression(col)(row)(pos).contains("!")) out = 1;
   return out;

}

void SxKdotP::removeOpMarks(int pos, int col, int row)
{

   expression(col)(row)(pos) = expression(col)(row)(pos)
      .substitute("?", "");
   if (leftPtr(col)(row)(pos) > 0)
      removeOpMarks(leftPtr(col)(row)(pos), col, row);
   if (rightPtr(col)(row)(pos) > 0)
      removeOpMarks(rightPtr(col)(row)(pos), col, row);
}

void SxKdotP::correctOperators(int col, int row)
{
   whatIsElement(0, col, row);
   removeOpMarks(0, col, row);
}

void SxKdotP::countOperatorMultiplications(SxString in)
{
   int length = (int)in.getSize();
   for (int i = 0; i < length-1; i++)
      if ((in(i) == '*') && (in(i+1) == 'k'))
   nOpMult++;
}

void SxKdotP::buildTree (int col, int row, SxString elem, SxString ham)
{
   // --- substitute all unknown variables in Hamiltonian element, if possible
   SxString elemFull = replaceUnknown(elem, ham);
   // --- only in accurateInterface mode: check for operator in brackets
   if (accurateInterfaces)  {
      countOperatorMultiplications(elemFull.substitute(" ", ""));
   }
   // --- routine that actually sets up the tree
   resolveExpr(elemFull, col, row);
   /* Problem: operators kx, ky, kz automatically apply to Psi.
    Constants need multiplication with Psi if not somehow multiplied
    with an operator.*/
   correctOperators(col,row);
   // --- print out tree, UNCOMMENT FOR DEBUGGING PURPOSES
/*   cout << "List of expressions, tree " << col << ", " << row << ": " << endl;
   for (int i = 0; i < expression(col)(row).getSize(); i++)  {
      cout << i << ":\t" << expression(col)(row)(i);
      if (leftPtr(col)(row)(i) > -1)
         cout << "\t" << leftPtr(col)(row)(i);
      if (rightPtr(col)(row)(i) > -1)
         cout << "\t" << rightPtr(col)(row)(i);
     cout <<endl;
   }*/
//   SX_EXIT;
}

void SxKdotP::writeUpdatedParameters ()
{
   cout << "+-----------------------------------------------------------------------------" << endl;
   cout << "| Writing updated parameter files..." << endl;
   cout << "| Material: " << matNames(nMat) << endl;
   SxString outfile = "" + matNames(nMat) + ".updated.sx";
   SX_OUTPUT_FILE (outfile); 
   try  {
      SxFileIO outio;
      outio.open (outfile, "a");
      SxString outLine;
      outLine = "parameterSet {\n"; outio.write(outLine);
      outLine = "   material=\"" + matNames(nMat) + "\";\n\n";
      outio.write(outLine);
      for (int iParam = 0; iParam < nParam; iParam++)  {
         outLine = "   parameter { name = \"" + paramNames(iParam)
                 + "\"; value = " << SxString(matParam(nMat, iParam))
                 + ";}\n";
         outio.write(outLine);
      }
      outLine = "\n}"; outio.write(outLine);
      outio.close ();
   } catch (const SxException &e) {
      e.print ();
      SX_EXIT;
   }

}

void SxKdotP::read (const SxSymbolTable *table)
{

   accurateInterfaces = true;
   derivatives = false;
   kDerivFirstRound = true;
   // --- read the material parameter files ---
   nMat = 0;
   nParam = 0;
   speedOpt = false;
   moreIO   = true;
   precMaterial = 0;
   SxSymbolTable *paramSet = NULL, *param = NULL, *materialMap = NULL,
                 *hamiltonian = NULL, *strain = NULL, *bandstructure = NULL,
                 *path = NULL, *extCharge = NULL, *fitData = NULL,
                 *plandscape = NULL;

   // --- loop over all materials parameter sets
   SxStack<double> paramVals, bowingVals, nonzeroFitParam, paramRange;
   SxMap<SxString,int> paramMap; // each parameter gets an integer index
   SxMap<SxString, int>::Iterator it;
   int iEpsR = -1;
   kPriority.resize(0);
   SxArray<double> bandGap(2);
   bool onlyFit = false;
   for (paramSet  = table->getGroup("parameterSet");
         paramSet != NULL;
         paramSet  = paramSet->nextSibling ("parameterSet"))

   {
      // --- first material defines parameter map name->integer
      if (nMat == 0)  {
         for (param  = paramSet->getGroup("parameter");
               param != NULL;
               param  = param->nextSibling ("parameter"))
         {
            SxString name = param->get("name")->toString ();
            if ((name == "k2") 
               || (name == "kx2")
               || (name == "ky2")
               || (name == "kz2")
               || (name == "kxy")
               || (name == "kxz")
               || (name == "kyz")
               || (name == "k")
               || (name == "kx")
               || (name == "ky")
               || (name == "kz")
               || (name == "kz")
               || (name == "e")
               || (name == "eXX")
               || (name == "eYY")
               || (name == "eZZ")
               || (name == "eXY")
               || (name == "eXZ")
               || (name == "eYZ")
               || (name == "Vp")
               || (name == "Vchg")
               || (name == "Vext"))  {
               cout << "parameter name is reserved keyword " << name
                    << "! EXITING." << endl;
               SX_EXIT;
            }
       
            paramMap(name) = nParam;
            paramNames << name;
            if (name == "epsilon") iEpsR = nParam;
            nParam++;
         }
      }

      SxString material = paramSet->get("material")->toString ();
      matNames << material;
      if (paramSet->contains("useForPreconditioner"))
         precMaterial = nMat;

      // --- browse parameter set for all parameter names
      for (it = paramMap.begin(); it != paramMap.end(); it++)  {
         bool found = false; // parameter found in material parameter group?
         param = paramSet->getGroup("parameter");
         while ((!found) && (param != NULL)) {
            SxString name = param->get("name")->toString();
            if (it.getKey() == name)
            {
               double val = param->get("value")->toReal();
               paramVals << val;
               if (param->contains("range"))  {
                  paramRange << param->get("range")->toReal();
               }
               else paramRange << 0;
               if (param->contains("nonzero"))  {
                  nonzeroFitParam << 1.0;
               }
               else  {
                  nonzeroFitParam << 0.0;
               }

               found = true;
               double bowing = 0.;
               if (param->contains("bowing"))  {
                  cout << "bowing employed" << endl;
                  bowing = param->get("bowing")->toReal();
               }
               bowingVals << bowing;

            }
            param = param->nextSibling ("parameter");
         }
         if (!found)  {
            cout << "Error: Parameter " << it.getKey() <<
               " was not found in the " << material << " parameter file."
               << endl << "EXITING." << endl;
            SX_EXIT;
         }
      }
      // --- is parameter fitting required?
      
      if (paramSet->containsGroup("fit"))  {
         fitParams << nMat;
         fitData = paramSet->getGroup("fit");
         
         SxString fileName = fitData->get("file")->toString ();
         fitFileNames << fileName;
         int fitSteps = fitData->get("nSteps")->toInt ();
         nFitSteps << fitSteps;
         int fitPoints = fitData->get("nPoints")->toInt ();
         nFitPoints << fitPoints;
         if (fitData->contains("bandPriority") ) {
            bandPriority = SxVector<double>(fitData->get("bandPriority")->toList());
            bandPrioritiesList << bandPriority;
         }
         if (fitData->contains("kPriority") ) {
            kPriority = SxVector<double>(fitData->get("kPriority")->toList());
            double kPrioWeight = 10.0;
            double kPrioFWHM = 10.0;
            if (fitData->contains("kPriorityWeight") ) {
               kPrioWeight = fitData->get("kPriorityWeight")->toReal();
            }
            if (fitData->contains("kPriorityFWHM") ) {
               kPrioFWHM = fitData->get("kPriorityFWHM")->toReal();
            }

            kPriorityWeights << kPrioWeight;
            kPriorityFWHMs << kPrioFWHM;
            kPrioritiesList << kPriority;
         }
         bandGap(0) = 0;
         bandGap(1) = 0;
         double insideGapPrice = 0; // default
         int BZSampling = 0; // default
         if (fitData->contains("bandGap") ) {
            bandGap = fitData->get("bandGap")->toList();
            insideGapPrice = 1;
            BZSampling = 5;
         }
         eVBList << bandGap(0);
         eCBList << bandGap(1);

         if (fitData->contains("insideGapPrice"))  {
            insideGapPrice = fitData->get("insideGapPrice")->toReal();
         }
         if (fitData->contains("onlyFit")) onlyFit = true;
         gapPriceList << insideGapPrice;
         if (fitData->contains("BZSampling"))  {
            BZSampling = fitData->get("BZSampling")->toInt();
         }
         BZSamplingList << BZSampling;


      }

      nMat++;
     }
   
   // --- read external charge file names and prefactors //TODO:NEW FEATURE
   int nCharges = 0;
   short dataType = -1;
   SxMeshR chargeTmp;

   SxList<double> chargePrefactors;
   if (table->containsGroup("charge"))  {
      for (extCharge  = table->getGroup("charge");
            extCharge != NULL;
            extCharge  = extCharge->nextSibling ("charge"))

      {
         SxString filename = extCharge->get("file")->toString ();
         double chgPrefactor = extCharge->get("prefactor")->toReal ();
         chargeFiles << filename;
         chargePrefactors << chgPrefactor;
         nCharges++;
      }
   }

   matParam = SxVector<double> (paramVals);
   matParam.reshape (nParam, nMat);
   matParam = matParam.transpose (); // make it nMat x nParam
   bowParam = SxVector<double> (bowingVals);
   bowParam.reshape (nParam, nMat);
   bowParam = bowParam.transpose (); // make it nMat x nParam
   fitParam = SxVector<double> (paramRange);
   fitParam.reshape (nParam, nMat);
   fitParam = fitParam.transpose (); // make it nMat x nParam
   nonzeros = SxVector<double> (nonzeroFitParam);
   nonzeros.reshape (nParam, nMat);
   nonzeros = nonzeros.transpose (); // make it nMat x nParam

   cout << "####### Input materials and parameters ########" << endl;
   cout << "Mat.\t  |  ";
   for (int iMat = 0; iMat < nMat; iMat++)
      cout << matNames(iMat) << "\t| ";
   cout << endl;
   for (int iParam = 0; iParam < nParam; iParam++)  {
      cout << paramNames(iParam) << " \t| ";
      for (int iMat = 0; iMat < nMat; iMat++)  {
         cout << matParam(iMat, iParam) << " \t| ";
      }
      cout << "\n";
   }

   // --- read the material map file ---
   const SxRBasis &R = *rBasisPtr;
   SxString inputASCII = "";
   SxString inputBinary = "";
   SxMesh3D mesh = R.getMesh();
   
   int xMax = mesh(0);
   int yMax = mesh(1);
   int zMax = mesh(2);
   int i, j;
   const int meshSize = xMax * yMax * zMax;

   materialMap  = table->getGroup("materialMap");
   
   if (materialMap->contains("asciiFile"))  {
      inputASCII = materialMap->get("asciiFile")->toString();
      dataType = 1;
   } else if (materialMap->contains("binaryFile"))  {
      inputBinary = materialMap->get("binaryFile")->toString();
      dataType = 0;
   }

   hamiltonian  = table->getGroup("kpHamiltonian");
   if (hamiltonian->contains("simplifiedInterfaces"))
      accurateInterfaces = false;
   // --- optimize for speed?
   if (hamiltonian->contains("speedOpt"))
      speedOpt = true;

   if (hamiltonian->contains("lessIO"))
      moreIO = false;

   materials.resize(nMat);
   for (int iMat = 0; iMat < nMat; iMat++)  {
      materials(iMat).resize(meshSize);
      materials(iMat).setBasis (rBasisPtr);
   }
   zero.resize(meshSize);
   one.resize(meshSize);
   zero.set (0.);
   one.set (1.);
   int ok = 0; // nr of mesh points that have a total composition of 1
   SxMatrix3<double> cell = R.cell;
   if (dataType == 0)  {
      SxBinIO io (inputBinary, SxBinIO::BINARY_READ_ONLY);
      materials = io.readMesh (&cell, &mesh);
      io.close();
   } else if (dataType == 1)  {
      int line = 0;
      SxFileParser fp (inputASCII.ascii());
      SxVector3<int> pos;
      for (int x = 0; x < xMax; x++)  {
         pos(0) = x;
         for (int y = 0; y < yMax; y++)  {
            pos(1) = y;
            for (int z = 0; z < zMax; z++)  {
               pos(2) = z;
               int meshIdx = (int)mesh.getMeshIdx(pos, SxMesh3D::Positive);
               line++;
               double checksum = 0.;
               for (int iMat = 0; iMat < nMat; iMat++)  {
                  double value = fp.getDouble ();
                  materials(iMat)(meshIdx) = value;
                  checksum += value;
               }
               if (fabs(checksum - 1.) < 1.e-5) ok++;
               else if (fabs(checksum - 1.) > 1.e-5)
                  cout << "checksum error at (" << x << ", " << y << ", "
                  << z << "), diff = " << fabs(1. - checksum) << ", line: "
                  << line << endl;
            }
         }
      }
   }
   // --- check if size of map fits mesh dimensions
   if (materials(0).getSize() != meshSize)  {
      cout << "Mesh size does not fit map size! EXITING." << endl;
      SX_EXIT;
   }
   if (speedOpt)  {
      cout << "optimized for speed -> higher memory consumption" << endl;
      pMem.resize(nParam);
      speedOpt = false; // allows to use SxMeshR SxKdotP::parameters(iParam)
      for (int iParam = 0; iParam < nParam; iParam++)  {
         pMem(iParam) = parameters(iParam);
      }
      speedOpt = true;
   }
   else  {
      cout << "optimized for memory -> higher computational time" << endl;
   }

   int iComp;

   // --- read energy where to search for states
   eTrial = hamiltonian->get("eTrial")->toReal();
   int nBands = hamiltonian->get("nBands")->toInt();
   cout << "eTrial = " << eTrial << " Hartree, nBands = " << nBands << endl;

   cout << "open Hamiltonian file ";
   SxString hamFile = hamiltonian->get("hamFile")->toString();
   cout << hamFile << "... ";
   FILE *file = sxfopen(hamFile, "r");
   cout << "success." << endl;
   cout << "Build Hamiltonian tree...";

   SxStack<char> fileContent; // file content without white space characters
   while (!feof(file))  {
      int c = fgetc(file);
      if (c == EOF) break;
      char ch = (char)c;
      // check for whitespace
      if ((ch == '\n') || (ch == '\t') || (ch == ' ')) continue;
      // --- check for possible /* ... */ comment
      if (ch == '/')  {
         c = fgetc(file);
         if (c == '*')  {
            // --- remove /* ... */ comments
            ch = 0;
            while (fscanf (file, "%*[^*]*%c", &ch) == 1)  {
               if (ch == '/') break;
               ungetc (ch, file);
            }
            if (ch != '/')  {
               cout << "Comment not closed in Hamiltonian! Exiting." << endl;
               SX_QUIT;
            }
            continue;
         } else if (c != EOF)  {
            ungetc (c, file);
         }
      }
      fileContent << ch;
   }
   fclose (file);
   SxString inFile;
   inFile.resize (fileContent.getSize ());
   fileContent.exportStack (inFile.elements, inFile.getSize ());
   cout << "success." << endl;

   inFile = inFile.substitute("Sqrt", "_");
   inFile = inFile.substitute("Sin", "!");
   inFile = inFile.substitute("Cos", "?");
   inFile = inFile.substitute("Tan", ".");
   // --- allow treatment of A - B - C terms
   inFile = inFile.substitute("-", "+(-1.)*"); 
   inFile = inFile.substitute("=+", "="); 
   inFile = inFile.substitute(",+", ","); 
   inFile = inFile.substitute("[+", "[");
   // --- check wether parameter names are overwritten
   for (int pIdx = 0; pIdx < paramNames.getSize(); pIdx++)  {
      if (inFile.contains(paramNames(pIdx) + "=") > 0) {
         cout << "Error: Parameter " << paramNames(pIdx)
              << " must not be overwritten.\nPlease rename the parameter "
              << "in either the parameter file or the Hamiltonian." << endl;
         SX_QUIT;
      }
   }
   // ------------------------------------------------

   SxString hamString = ((inFile.right("Hamiltonian")).right("=[")).left("];");
   if (nBands != int(sqrt(1. * (hamString.contains(",") + 1))))  {
      cout << "Number of bands is not consistent with provided Hamiltonian. "
           << endl << "nBands = " << nBands << endl
           << "Elements in Hamiltonian: " << hamString.contains(",") + 1
           << ", should be " << nBands * nBands << ". EXITING." << endl;
      SX_EXIT;
   }
   
   // --- substitute keywords by symbols
   SxList<SxString> cols = hamString.tokenize(']');
   int nCols = (int)cols.getSize();
   // --- split Hamiltonian in columns
   nComp = nCols;
   cout << nComp << "-band Hamiltonian used here..." << endl;
   // --- set up matrix with Hamiltonian elements
   expression.resize(nComp);
   leftPtr.resize(nComp);
   rightPtr.resize(nComp);
   for (iComp = 0; iComp < nComp; iComp++)  {
      expression(iComp).resize(nComp);
      leftPtr(iComp).resize(nComp);
      rightPtr(iComp).resize(nComp);
   }
   nOpMult = 0;
   cout << "cols: " << nCols << endl;
   for (i = 0; i < nCols; i++)  {
      cols(i) = cols(i).right("[");
      SxList<SxString> rows = cols(i).tokenize(',');
      int nRows = (int)rows.getSize();
      if ( i == 0) cout << "rows: " << nRows << endl;
      if (nRows != nCols) {
         cout << "Hamiltonian is no square matrix. EXITING." << endl;
         cout << "nRows: " << nRows << endl;
         cout << "nCols: " << nCols << endl;
         SX_EXIT;
      }
      // --- split columns in single elements
      for (j = 0; j < nRows; j++)  {
         // --- evaluate single element
         buildTree(i, j, rows(j), inFile);
      }
   }
   cout << "tree set up." << endl;
   if (accurateInterfaces)  {
      cout << "Resize arrays for operator multiplications" << endl;
      opKey.resize(nOpMult); // int that represents operator:
      //0: k, 1: kx, 2: ky, 3: kz, 4: k2, 5: kx2, 6: ky2, 7: kz2
      par.resize(nOpMult);
      kIPar.resize(nOpMult);
      kJPar.resize(nOpMult);
      kKPar.resize(nOpMult);
   }
   // generate byte code
   byteCodeHam = generateByteCode ();

   // --- input for bandstructure plot
   bool onlyBS = false;
   bool onlyPL = false;
   int outputPar = 0;
   outPar = ""; // empty
   if (hamiltonian->contains("outputParameter"))  {
      outPar = hamiltonian->get("outputParameter")->toString();
      for (int pIdx = 0; pIdx < paramNames.getSize(); pIdx++)
         if (outPar == paramNames(pIdx)) {
            outputPar = pIdx;
         }
   }
   outMesh = parameters(outputPar);

   // --- read polarization potential.
   SxArray<SxMeshR> vPolR, vExtR, strainR, chargeR;
   SxString polFile;
   vP = zero;
   if (hamiltonian->contains("polarization"))  {
      cout << "read polarization potential..." << endl;

      dataType = -1;
      polFile = hamiltonian->get("polarization")->toString();
      if (polFile.contains(".sxb")) dataType = 0;
      else if (polFile.contains(".dat")) dataType = 1;
      else {
         cout << polFile
            << " has unknown data format. Please provide .dat or.sxb file."
            << " Exiting." << endl;
         SX_EXIT;
      };

      if (dataType == 0)  {
         SxBinIO io (polFile, SxBinIO::BINARY_READ_ONLY);
         vPolR.resize (1);
         vPolR = io.readMesh (&cell, &mesh);
         vP = (vPolR(0)*(1./HA2EV));
         io.close();

      } else if (dataType == 1)  {
         // --- read ascii polarization file
         SxFileParser fp(polFile);
         vP.resize(meshSize);
         for (int x = 0; x < xMax; x++)  {
            for (int y = 0; y < yMax; y++)  {
               for (int z = 0; z < zMax; z++)  {
                  ssize_t meshIdx = mesh.getMeshIdx(x,y,z, SxMesh3D::Positive);
                  vP(meshIdx) = fp.getDouble () / HA2EV;
               }
            }
         }
      }
      // --- check if mesh size reflects provided data
      if (meshSize != vP.getSize())  {
         cout <<
            "Polarisation potential size does not match mesh size! EXITING."
            << endl;
         SX_EXIT;
      }
   }

   // --- read external potential.
   SxString extFile;
   vExt = zero;
   if (hamiltonian->contains("extPotential"))  {
      dataType = -1;
      extFile = hamiltonian->get("extPotential")->toString();
      if (extFile.contains(".sxb")) dataType = 0;
      else if (extFile.contains(".dat")) dataType = 1;
      else {
         cout << extFile
            << " has unknown data format. Please provide .dat or.sxb file."
            << " Exiting." << endl;
         SX_EXIT;
      }
      if (dataType == 0)  {
         SxBinIO io (extFile, SxBinIO::BINARY_READ_ONLY);
         vExtR.resize (1);
         vExtR = io.readMesh (&cell, &mesh);
         vExt = (vExtR(0)*(1./HA2EV));
         io.close();
      } else if (dataType == 1)  {
         // --- read ascii polarization file
         SxFileParser fp (extFile);
         SxVector3<int> pos;
         vExt.resize(meshSize);
         for (int x = 0; x < xMax; x++)  {
            pos(0) = x;
            for (int y = 0; y < yMax; y++)  {
               pos(1) = y;
               for (int z = 0; z < zMax; z++)  {
                  pos(2) = z;
                  int meshIdx = (int)mesh.getMeshIdx(pos, SxMesh3D::Positive);
                  vExt(meshIdx) = fp.getDouble () / HA2EV;
               }
            }
         }
      }
      if (meshSize != vExt.getSize())  {
         cout << "External potential size does not match mesh size! EXITING."
              << endl;
         SX_EXIT;
      }
   }

   SxVecRef<double> epsilonR;
   if (iEpsR > -1)
      epsilonR = parameters(iEpsR);

   // --- read strain fields.
   if (hamiltonian->containsGroup("strain"))  {
      cout << "read strain fields..." << endl;
      strain = hamiltonian->getGroup("strain");
      SxArray<SxString> strainInput; // 3 diagonal, 3 off-diagonal strains
      dataType = -1; // 0 - binary, 1 - ascii
      strainInput.resize(6);
      eIJ.resize(6);
      strainInput(0) = strain->get("eXX")->toString();
      strainInput(1) = strain->get("eYY")->toString();
      strainInput(2) = strain->get("eZZ")->toString();
      strainInput(3) = strain->get("eXY")->toString();
      strainInput(4) = strain->get("eXZ")->toString();
      strainInput(5) = strain->get("eYZ")->toString();
      for (int idx = 0; idx < 6; idx++)  {
         if (strainInput(idx).contains(".sxb")) dataType = 0;
         else if (strainInput(idx).contains(".dat")) dataType = 1;
         else {
            cout << strainInput(idx)
               << " has unknown data format. Please provide .dat or.sxb file."
               << " Exiting." << endl;
            SX_EXIT;
         };
         if (dataType == 0)  {
         // --- read binary strain file
            SxBinIO io (strainInput(idx), SxBinIO::BINARY_READ_ONLY);
            strainR.resize (1);
            strainR = io.readMesh (&cell, &mesh);
            eIJ(idx)= strainR(0);
            cout << "open " << strainInput(idx) << endl;
            io.close();

         } else if (dataType == 1)  {
         // --- read ascii strain file
            SxFileParser fp(strainInput(idx));
            SxVector3<int> pos;
            eIJ(idx).resize(meshSize);
            for (int x = 0; x < xMax; x++)  {
               pos(0) = x;
               for (int y = 0; y < yMax; y++)  {
                  pos(1) = y;
                  for (int z = 0; z < zMax; z++)  {
                     pos(2) = z;
                     int meshIdx = (int)mesh.getMeshIdx(pos, SxMesh3D::Positive);
                     eIJ(idx)(meshIdx) = fp.getDouble ();
                  }
               }
            } // end x,y,z loops
         }
      
         if (meshSize != eIJ(idx).getSize())  {
            cout << "Strain field size does not match mesh size! EXITING."
                 << endl;
            SX_EXIT;
         }
  
      } // end loop over strains 
   }
   // --- read external charges
   SxMeshR totalCharge;
   for (int idx = 0; idx < nCharges; idx++)  {
      if (chargeFiles(idx).contains(".sxb")) dataType = 0;
         else if (chargeFiles(idx).contains(".dat")) dataType = 1;
         else {
            cout << chargeFiles(idx)
               << " has unknown data format. Please provide .dat or.sxb file."
               << " Exiting." << endl;
            SX_EXIT;
         };
         if (dataType == 0)  {
         // --- read binary external charge file
            SxBinIO io (chargeFiles(idx), SxBinIO::BINARY_READ_ONLY);
            chargeR.resize(1);
            chargeR = io.readMesh (&cell, &mesh);
            chargeTmp = (chargeR(0));
            cout << "open " << chargeFiles(idx) << endl;
            io.close();

         } else if (dataType == 1)  {
         // --- read ascii external charge file
            SxFileParser fp(chargeFiles(idx));
            SxVector3<int> pos;
            chargeTmp.resize(meshSize);
            for (int x = 0; x < xMax; x++)  {
               pos(0) = x;
               for (int y = 0; y < yMax; y++)  {
                  pos(1) = y;
                  for (int z = 0; z < zMax; z++)  {
                     pos(2) = z;
                     int meshIdx = (int)mesh.getMeshIdx(pos, SxMesh3D::Positive);
                     chargeTmp(meshIdx) = fp.getDouble ();
                  }
               }
            } // end x,y,z loops
         }
      
         if (meshSize != chargeTmp.getSize())  {
            cout << "external charge map size does not match mesh size! EXITING."
                 << endl;
            SX_EXIT;
         }
         if (iEpsR < 0)  {
            cout << "Error: no permittivity parameters given. "
                 << "Parameter `epsilon` required in all material files, "
                 << "alternatively, remove charges. EXITING." << endl;
            SX_EXIT;
         }
         if (idx == 0)  { totalCharge  = chargeTmp * chargePrefactors(idx); }
         else           { totalCharge += chargeTmp * chargePrefactors(idx); }
      }
   
   cout << "found the following charge density files: " << endl;
   for (i = 0; i < nCharges; i++)  {cout << chargeFiles(i) << " * " << chargePrefactors(i) << endl;}
   cout << "total charge in system: " << totalCharge.sum() << endl;

   // --- scaled dimension ?
   wgt(0) = 1.;
   wgt(1) = 1.;
   wgt(2) = 1.;
   if (hamiltonian->contains("weight"))  {
      wgt = SxVector3<double>(hamiltonian->get("weight")->toList());
      for (i = 0; i < 3; i ++)
         wgt(i) = 1./wgt(i);
      }
   cout << "weighting k-vectors with: " << wgt << endl;

   // --- initialize chargePotential
   if (nCharges > 0 && chargePotential.getSize () == 0)  {
      const SxGBasis &G = *gBasisPtr;
      SxVector<double> scaledG2 (G);
      for (int ig = 0; ig < G.ng; ig++)  {
         scaledG2(ig) = sqr(G.gVec(ig,0) * wgt(0)) +
                        sqr(G.gVec(ig,1) * wgt(1)) +
                        sqr(G.gVec(ig,2) * wgt(2));
      }

      cout << "compute potential from external charges" << endl;
      double dVol = R.cell.volume / (double)mesh.product ();
      SxVector<PrecCoeffG> vG;
      if (totalCharge.getSize() > 0)  {
         PsiR chgR = totalCharge / epsilonR;
         // prefactors: 4pi      - prefactor of Coulomb potential
         //             wgt(...) - dimension scaling
         //             dVol     - switch from charge to charge density
         vG = (G | chgR) * FOUR_PI * (wgt(0) * wgt(1) * wgt(2) / dVol);
         vG(0) = 0;
         vG(SxIdx(1,G.ng-1)) /= scaledG2(SxIdx(1,G.ng-1));
         chargePotential = SxMeshR(R | vG);
         cout << "Charge Potential: " << chargePotential.minval() * HA2EV
              << " eV to " << chargePotential.maxval() * HA2EV << " eV."
              << endl << "Potential difference: "
              << (chargePotential.maxval() - chargePotential.minval() ) * HA2EV
              << " eV." << endl;
         cout << "write charge potential " << endl;

         ofstream vChgFile;
         fstream::openmode mode = fstream::out | fstream::trunc;
         vChgFile.open (("vChg.dat"), mode);
         double valChg;
         for (int x = 0; x < mesh(0); x++)  {
            for (int y = 0; y < mesh(1); y++)  {
               for (int z = 0; z < mesh(2); z++)  {
                  ssize_t idx = mesh.getMeshIdx(x,y,z, SxMesh3D::Positive);
                  valChg = chargePotential(idx) * HA2EV;
                  vChgFile << SxString(valChg) << endl;
               }
            }
         }
         vChgFile.close ();
         cout << "done. " << endl;
      }
   }
   // --- show bandstructure
   if (hamiltonian->containsGroup("bandstructure"))  {
      int iStep = 0;
      SxString bsOutFile;
      SxVector3<int> rCoord;
      bandstructure = hamiltonian->getGroup("bandstructure");
      cout << "+------ Band structure calculation requested. -------" << endl;
      if (bandstructure->contains("onlyBS")) onlyBS = true;
      bsOutFile = bandstructure->get("outFile")->toString();
      rCoord = SxVector3<int>(bandstructure->get("rCoord")->toIntList());
      for (i = 0; i < 3; i++)
         if (rCoord(i) >= mesh(i))  {
            cout << "| Real space coordinate is not inside mesh." << endl;
            SX_EXIT;
         }
      cout << "| Calculate band structure at real space grid point " << rCoord
           << endl;
      cout << "| Material composition at this point: " << endl << "| ";
      double composition;
      int meshIdx = (int)mesh.getMeshIdx(rCoord, SxMesh3D::Positive);
      for (int iMat = 0; iMat < nMat; iMat++)  {
         composition = materials(iMat)(meshIdx);
         cout << matNames(iMat) << ": " << composition << "\t";
      }
      cout << endl;

      SxVector3<double> step, k;
      cout << "| Set up path for band structure" << endl;
      SxStack<SxVector3<double> > kStack;
      for (path  = bandstructure->getGroup("path");
               path != NULL;
               path  = path->nextSibling ("path"))
         {
            bsStart   = SxVector3<double>(path->get("from")->toList());
            bsEnd     = SxVector3<double>(path->get("to")->toList());
            int stepsBs   = path->get("steps")->toInt();
            step = 1. / stepsBs * (bsEnd - bsStart);
            // --- set up full path for bandstructure
            for (iStep = 0; iStep < stepsBs; iStep++)  {
               k = bsStart + iStep * step;
               kStack << k;
            }
         }
      k = bsStart + iStep * step;
      kStack << k;
      bsPts = kStack;
      cout << "| Calculate band structure... " << endl;
      showBandstructure(bsOutFile, rCoord);
      if (onlyBS) {
         printTiming ();
         cout << "Only band structure calculation requested.\n"
              << "Exiting now." << endl; 
         exit(0); 
      }
   }

   // --- potentialLandscape
   if (hamiltonian->containsGroup("potentialLandscape"))  {
      SxString plOutFile;
      plandscape = hamiltonian->getGroup("potentialLandscape");
      cout << "+------ Potential landscape calculation requested. -------" << endl;
      if (plandscape->contains("onlyPL")) onlyPL = true;
      plOutFile = plandscape->get("outFile")->toString();
      plLayer(0) = -1;
      plLayer(1) = -1;
      plLayer(2) = -1;
      if (plandscape->contains("xCoord"))
         plLayer(0) = plandscape->get("xCoord")->toInt ();
      if (plandscape->contains("yCoord"))
         plLayer(1) = plandscape->get("yCoord")->toInt ();
      if (plandscape->contains("zCoord"))
         plLayer(2) = plandscape->get("zCoord")->toInt ();

      SxVector3<double> kvals;
      kvals(0) = 0;
      kvals(1) = 0;
      kvals(2) = 0;
      if (plandscape->contains("kval")) 
         kvals = SxVector3<double>(plandscape->get("kval")->toList());

      cout << "Computing potential landscape at k = (" << kvals(0) << ", "
              << kvals(1) << ", " << kvals(2) << ").\n";
      potentialLandscape(plOutFile,kvals);
      if (onlyPL) {
         printTiming ();
         cout << "Only potential landscape calculation requested.\n"
              << "Exiting now." << endl; 
         exit(0); 
      }
   }
   // --- info regarding parameter fit
   int nFitParam = int(fitFileNames.getSize());
   if (nFitParam > 0) {
      matParamOriginal = matParam;  // save the original
      cout <<  "+-----------------------------------------------------------------------------" << endl;
      cout <<  "| Parameter fit requested: " << endl;
      cout <<  "+-----------------------------------------------------------------------------" << endl;
      int iFitParam = 0;
      for (iFitParam = 0; iFitParam < nFitParam; iFitParam++)  {
         nMat = fitParams(iFitParam);
         SxString material = matNames(iFitParam);
         SxString fileName = fitFileNames(iFitParam);
         int fitSteps = nFitSteps(iFitParam); 
         int fitPoints = nFitPoints(iFitParam); 
         // --- load file, read number of lines and columns, extract k's and nBands
         SxFileParser fp(fileName);
         SxStack<double> data;
         // --- first line: determine number of columns
         int inputBSBands = -3; // columns = 3x k-vector + bands
         while (!fp.endOfLine ())  {
            data << fp.getDouble ();
            inputBSBands++;
         }
         int nKValues = 1;
         // --- now read lines expecting the same number of columns
         //     (without checking, though)
         while (!feof (fp.fp))  {
            for (int iBand = 0; iBand < inputBSBands + 3; iBand++)
               data << fp.getDouble ();
            fp.skipWhite ();
            nKValues++;
         }
         fp.close ();
         SxVector<double> inputBandstructure = data;
         // current ordering is each k-point in one matrix column
         inputBandstructure.reshape (3 + inputBSBands, nKValues);
         // ... so transpose to have k-points in rows
         inputBandstructure = inputBandstructure.transpose ();
         data.removeAll (); // throw away temp. data

         double kWgt = 1;
         if (kPriority.getSize() > 0)
         {
            kWgt = kPriorityWeights(iFitParam);
         }

         cout << "| Fit required for " << material << "." << endl;
         cout << "| Input file: " << fileName << " has " << nKValues
                 << " k-values with " << inputBSBands << " bands. " << endl;
         cout << "| Number of iterations: " << fitSteps
                 << ", number of points in Sobol sequence: " << fitPoints
                 << "." << endl;
         cout << "| Band priorities: " << bandPriority << endl;
         if (kPriority.getSize() > 0)  {
            cout << "| k-priorities: " << kPriority << endl;
            cout << "| k-priority weight: " << kWgt << endl;


         }
         fitBandstructure(inputBandstructure, nKValues, inputBSBands);
         cout << "+-----------------------------------------------------------------------------" << endl;
      }
      writeUpdatedParameters ();
      printTiming (!onlyFit); // final timing if onlyFit, otherwise continue
      cout << "| Parameter fit finished." << endl;
      if (onlyFit) {
         cout << "| Exiting now." << endl;
         exit(0); 
      }

   }

}

void SxKdotP::fitBandstructure (SxVector<double> inBands, int nKs, int nBands)
{
   SX_CLOCK (Timer::FitBandsTotal);
   cout << "| Fitting band structure for " << nMat << "'th material." << endl;
   cout << "| Write input band structure to inputBands.dat." << endl;

   SxVector<double> params(nParam);
   // strain
   SxVector<double> strain(6);
   strain.set (0.);
   // potentials
   double vPol = 0., vExtern = 0., vCharge = 0.;

   // --- write input bandstructure
   int ik, ib;
   SxString outfile;
   outfile="inputBS_" + matNames(nMat) + ".dat";
   SX_OUTPUT_FILE (outfile); 
   try  {
      SxFileIO outio;
      outio.open (outfile, "a");
      for (ik = 0; ik < nKs; ik++) {
      
         SxString outLine = "" + SxString(inBands(ik,0)) + "  " + 
            SxString(inBands(ik,1)) + "  " + SxString(inBands(ik,2));
         for (ib = 0; ib < nBands; ib++)
            outLine = outLine + "   "
               + SxString(inBands(ik,ib+3));
         outLine = outLine + "\n";
         outio.write (outLine);
      }
      outio.close ();
   } catch (const SxException &e) {
      e.print ();
      SX_EXIT;
   }

         bandPriority = bandPrioritiesList(nMat);
         if (bandPriority.getSize() != nBands) {
            cout <<
            "| No band priorities given or number not matching input bands."
            << endl << "| All priorities set to 1." << endl;
            bandPriority.resize(nBands);
            for (ib = 0; ib < nBands; ib++)  {
               bandPriority(ib) = 1;
            }

         }
         kPriorities.resize(nKs);
         for (ik = 0; ik < nKs; ik++)  {  // default value 1 for all k
            kPriorities(ik) = 0;
         }

         kPriority = kPrioritiesList(nMat);
         double kWgt = 1;
         if (kPriority.getSize() > 0)
         {
            cout << "| setting up k-wise scaling..." << endl;
            kWgt = kPriorityWeights(nMat);
            int iKPts;
            double xk;
            for (iKPts = 0; iKPts < kPriority.getSize(); iKPts++)  {
               for (ik = 0; ik < nKs; ik++)  {
                  xk = abs(kPriority(iKPts)-ik);
                  if (kWgt > 1)  {kPriorities(ik) = kPriorities(ik) + exp(-xk*xk/kWgt);}
                  else {kPriorities(ik) = 1;}
               }
            }
         }

   // --- compute bandstructure using initial parameters
   outfile="initBS_" + matNames(nMat) + ".dat";
   SxVector3<double> kVector; // list of all points of bandstructure
   SxVector<SxComplex16> hMat;
   SxVector<SxComplex16> eigs;
   double initialDiff = 0;
   double insideGap;
   double insideGapPrice = gapPriceList(nMat);
   int discrBZ = BZSamplingList(nMat);
   double maxKx = 0., maxKy = 0., maxKz = 0.;
   for (ik = 0; ik < nKs; ik++) {  // determine max k values in band structure
      if (abs(inBands(ik,0)) > maxKx) {maxKx = abs(inBands(ik,0));}
      if (abs(inBands(ik,1)) > maxKy) {maxKy = abs(inBands(ik,1));}
      if (abs(inBands(ik,2)) > maxKz) {maxKz = abs(inBands(ik,2));}
   }
   double VBO = eVBList(nMat);
   double CBO = eCBList(nMat);
   cout << "+-----------------------------------------------------------------------------" << endl;
   // --- sanity checks for brute-force ellipticity tests
   if (VBO > CBO) {
      cout << "| >>> ERROR: negative band gap in fitting process." << endl
           << "| >>> Pease correct order of band gap values (lowest first)."
           << endl << "| >>> EXITING." << endl;
      SX_EXIT;
   }
   if ((VBO < CBO && discrBZ == 0)) {
      cout << "| >>> WARNING: number of discretisation steps in BZ = 0." << endl
           << "| >>> Energies inside band gap will be ignored during fit." << endl;
   }
   if (VBO < CBO && insideGapPrice == 0) {
      cout << "| >>> WARNING: price for energies in gap = 0." << endl
           << "| >>> Energies inside band gap will be ignored during fit." << endl;
   }

   if (VBO == 0 && CBO == 0) {
      cout << "| States in band gap not penalised." << endl;
   }
   else  {
      cout << "| Forbidden energy range: " << VBO << " to " << CBO << endl;
      cout << "| Sampling points in BZ to be evaluated: " << discrBZ << endl;
      cout << "| Price for energies inside gap: " << insideGapPrice << endl;
   }
   double percentageInsideGap = 0;
   // --- brute-force ellipticity check

   SxStack<SxComplex16> workspace;

   SX_OUTPUT_FILE (outfile); 
   try  {
      SxFileIO outio;
      outio.open (outfile, "a");
      int iPar;
      for (iPar = 0; iPar < nParam; iPar++)  {
         params(iPar) = matParam(nMat,iPar) = matParamOriginal(nMat,iPar);
      }
      for (ik = 0; ik < nKs; ik++) {
         kVector(0) = inBands(ik,0);
         kVector(1) = inBands(ik,1);
         kVector(2) = inBands(ik,2);
         hMat = hMatrixBS(kVector, strain, vPol, vExtern, vCharge, params,
                          workspace);

         SxString outLine = "" + SxString(kVector(0)) + "  " + 
            SxString(kVector(1)) + "  " + SxString(kVector(2));
         eigs = symEigenvalues (hMat);
         for (ib = 0; ib < nBands; ib++)  {
            outLine = outLine + "   "
               + SxString(eigs(ib).re * HA2EV);
            initialDiff = initialDiff + bandPriority(ib)
                    * kPriorities(ik) 
                    * abs(eigs(ib).re * HA2EV - inBands(ik,ib+3));

         }

         // --- count energies inside gap: brute-force ellipticity

         outLine = outLine + "\n";
         outio.write (outLine);
      }
      outio.close ();
      if (discrBZ > 0 && insideGapPrice > 0)   {
         insideGap = 0;
         for (int ix = 0; ix <= discrBZ; ix++)  {
            kVector(0) = 0 + ix * maxKx / double (discrBZ);
            for (int iy = 0; iy <= discrBZ; iy++)  {
               kVector(1) = 0 + iy * maxKy / double (discrBZ);
               for (int iz = 0; iz <= discrBZ; iz++)  {
                  kVector(2) = 0 + iz * maxKz / double (discrBZ);
                  
                  hMat = hMatrixBS(kVector, strain, vPol, vExtern, vCharge,
                                   params, workspace);

                  eigs = symEigenvalues (hMat);
                  for (ib = 0; ib < nBands; ib++)  {
                     if (   eigs(ib).re * HA2EV > VBO
                         && eigs(ib).re * HA2EV < CBO) insideGap++;
                  }
               }
            }
         }
         int nBZ = (discrBZ+1) * (discrBZ+1) * (discrBZ+1);
         percentageInsideGap = insideGap / double(nBZ * nBands);
         cout << "| Eigenvalues inside gap, initial: " << percentageInsideGap*100 << " %" << endl;
         initialDiff = initialDiff * (1 + percentageInsideGap * insideGapPrice);
      }


   } catch (const SxException &e) {
      e.print ();
      SX_EXIT;
   }
   cout << "+-----------------------------------------------------------------------------" << endl;

   // --- determine nParamFit of fit parameter space
   int iParam;
   int paramNr;
   SxList<int> paramsToFit;
   SxList<double> rangeForFit;
   for (iParam = 0; iParam < nParam; iParam++)  {
      if (fitParam(nMat,iParam) != 0) {
         paramsToFit << iParam;
         rangeForFit << fitParam(nMat,iParam);
      }
   }
   int nParamFit = int(paramsToFit.getSize());
   cout << "| Initial diff: " << initialDiff << endl;
   cout << "| Dimension of parameter space to fit: " << nParamFit << endl;
   cout << "| Parameters to fit: " << endl;
   for (iParam = 0; iParam < nParamFit; iParam++)  {
      paramNr = paramsToFit(iParam);
      cout << paramNames(paramNr) << " = " << matParam(nMat, paramNr) << " +- "
              << rangeForFit(iParam) << endl;
   }
   
   // --- set up point sequence
   int fitPoints = nFitPoints(nMat); 
   int fitSteps  = nFitSteps(nMat); 
   SxVector<double> PointSequence(fitPoints,nParamFit);

   int iPts, iDim;

   // --- Sobol sequence 
   long long int seed = 1;
   SxVector<double> r(nParamFit);
   for ( iPts = 0; iPts < fitPoints; iPts++ )
   {
      i8_sobol ( nParamFit, &seed, r.elements );
      PointSequence.rowRef (iPts) = r;
   }
   SX_VALIDATE_VECTOR(PointSequence);
   
   time_t start = time(0);

   int iFitStep;
   int minPos = -1; // position of minimum in point set
   double minDiff = initialDiff; // best so far: initial values
   double percGapAtMin = -1.;
   SxVector<double> searchSpace(fitPoints,nParamFit);
   double span = 1.0; // parameter to reduce search span once frame is fix

   SX_START_TIMER (Timer::FitBands);
   for (iFitStep = 0; iFitStep < fitSteps; iFitStep++)  { // loop over iterations
      cout << "step: " << iFitStep << ": ";
      // --- map point sequence to parameter space
      for (iPts = 0; iPts < fitPoints; iPts++)  {
         for (iDim = 0; iDim < nParamFit; iDim++)  {
            paramNr = paramsToFit(iDim);
            searchSpace(iPts,iDim) = matParam(nMat,paramNr)
                    - span * rangeForFit(iDim)
                    + (2.0 * span) * rangeForFit(iDim) * PointSequence(iPts,iDim);
            if ((nonzeros(nMat,paramNr) == 1) && (searchSpace(iPts,iDim) == 0)) {
               searchSpace(iPts,iDim)=1.e-4;
            }
            SX_CHECK_NUM (searchSpace(iPts,iDim));
         }
      }

   
      // --- determine best suited set
      minPos = -1;
      for (iPts = 0; iPts < fitPoints; iPts++)  {  // loop over all points in sequence
         double diff = 0.;
         for (iDim = 0; iDim < nParamFit; iDim++)  {  // update all parameters that require fit
            paramNr = paramsToFit(iDim);
            params(paramNr) = searchSpace(iPts, iDim);
         }  

         for (ik = 0; ik < nKs; ik++) {  // compare band structures
            kVector(0) = inBands(ik,0);
            kVector(1) = inBands(ik,1);
            kVector(2) = inBands(ik,2);
            hMat = hMatrixBS(kVector, strain, vPol, vExtern, vCharge, params,
                             workspace);


            SX_START_TIMER (Timer::FitBandsDiag);
            eigs = symEigenvalues (hMat);
            SX_STOP_TIMER (Timer::FitBandsDiag);

            // --- compute diff for specific parameter set, k, iB
            for (ib = 0; ib < nBands; ib++)  {
               diff = diff + bandPriority(ib)
                       * kPriorities(ik)  
                       * abs(eigs(ib).re * HA2EV - inBands(ik,ib+3));
            }
         }

         if (discrBZ > 0 && insideGapPrice > 0 && diff <= minDiff) {
            // --- count energies inside gap: brute-force ellipticity
            insideGap = 0;
            for (int ix = 0; ix <= discrBZ; ix++)  {
               kVector(0) = 0 + ix * maxKx / double (discrBZ);
               for (int iy = 0; iy <= discrBZ; iy++)  {
                  kVector(1) = 0 + iy * maxKy / double (discrBZ);
                  for (int iz = 0; iz <= discrBZ; iz++)  {
                     kVector(2) = 0 + iz * maxKz / double (discrBZ);

                     hMat = hMatrixBS(kVector, strain, vPol, vExtern, vCharge,
                                      params, workspace);
                     SX_START_TIMER (Timer::FitBandsDiag);
                     eigs = symEigenvalues (hMat);
                     SX_STOP_TIMER (Timer::FitBandsDiag);

                     for (ib = 0; ib < nBands; ib++)  {
                        if (   eigs(ib).re * HA2EV > VBO
                            && eigs(ib).re * HA2EV < CBO) insideGap++;
                     }
                  }
               }
            }
            int nBZ = (discrBZ+1) * (discrBZ+1) * (discrBZ+1);
            percentageInsideGap = insideGap / double(nBZ * nBands);
            diff = diff * (1 + percentageInsideGap * insideGapPrice);
         }

         if (diff <= minDiff)  {
            // if better than min-so-far: new min and mark point set
            minDiff = diff;
            minPos = iPts;
            percGapAtMin = percentageInsideGap;
         }
      }

      // --- move and/or reduce search interval size
      if (minPos > -1)  {
         for (iDim = 0; iDim < nParamFit; iDim++)  {  // update all parameters that require fit
            paramNr = paramsToFit(iDim);
            matParam(nMat, paramNr) = searchSpace(minPos, iDim);
         } 
      }
      if (minPos  <= 0)  {span = span * 0.5;} // if search frame does not move any more, reduce search interval around
      cout << "diff = " << minDiff << ", minPos: " << minPos << ", span: " << span << endl;
   }
   SX_STOP_TIMER (Timer::FitBands);
   time_t finish = time(0);
   time_t totaltime = finish-start;
   cout << "+-----------------------------------------------------------------------------" << endl;
   cout << "| best fit: minDiff = " << minDiff << " for " << minPos
        << "'th set in sequence. " << endl
        << "|           percentage in gap: " << percGapAtMin * 100 << endl
        << "|           total time: " << totaltime <<  " s." << endl;
   cout << "+-----------------------------------------------------------------------------" << endl;
   cout << "| Updated parameters: " << endl;
   for (iParam = 0; iParam < nParamFit; iParam++)  {
      paramNr = paramsToFit(iParam);
      cout << paramNames(paramNr) << " = " << matParam(nMat, paramNr) << endl;
   }
   cout <<  "+-----------------------------------------------------------------------------" << endl;

   cout << "| Write output files..." << endl;
   // --- compute bandstructure using fit parameters
   outfile="fitBS_" + matNames(nMat) + ".dat";
   SX_OUTPUT_FILE (outfile); 
   try  {
      SxFileIO outio;
      outio.open (outfile, "a");
      for (ik = 0; ik < nKs; ik++) {
         kVector(0) = inBands(ik,0);
         kVector(1) = inBands(ik,1);
         kVector(2) = inBands(ik,2);
         hMat = hMatrixBS(kVector, strain, vPol, vExtern, vCharge, params,
                          workspace);
        
         SxString outLine = "" + SxString(kVector(0)) + "  " + 
            SxString(kVector(1)) + "  " + SxString(kVector(2));
         eigs = symEigenvalues (hMat);

         for (ib = 0; ib < nBands; ib++)  {
            outLine = outLine + "   "
               + SxString(eigs(ib).re * HA2EV);
         }
         outLine = outLine + "\n";
         outio.write (outLine);
      }
      outio.close ();
   } catch (const SxException &e) {
      e.print ();
      SX_EXIT;
   }
   cout << "| Done. " << endl;
   // safe updated parameter files in Sphinx format

}

void SxKdotP::printEnergies () const
{
   cout << "printEnergies()" << endl;
   std::streamsize oldPrec = cout.precision ();
   cout.precision(3);
   const SxPWSet &waves = getWavesRef();
   int i, iComp, nStates;
   const SxRBasis &R = *rBasisPtr;
   nStates = waves.getNStates(); 
   cout << "composition of the wave functions of each state:" << endl;
   SxArray<SxVector<double> > rho;
   SxMatrix3<double> cell = R.cell;
   SxMesh3D mesh = R.getMesh();
   rho.resize(nStates);

   if (moreIO)  
   for (i=0; i < nStates; i++)  {
      rho(i) = zero;
      cout << " " << i << " :\t";

      PsiG psiI = waves(i,0,0); // TODO: PsiRef ?
      SxVector<SxComplex16> psiIR = (R | psiI);
      SxArray<PsiR> rComp(nComp);

      for (iComp = 0; iComp < nComp; iComp++)  {
           
         const SxVecRef<PrecCoeffG,SubMatrix> &comp = psiI.compRef(iComp);
         const SxVecRef<PrecCoeffG,SubMatrix> &compR = psiIR.compRef(iComp);
         rComp(iComp) = compR;
         
         cout << dot(comp, comp).re << "\t";
         // --- calculate charge density
         rho(i) += compR.absSqr();
      }

      SxBinIO ioRho ("rho-" + SxString(i) + ".sxb",
            SxBinIO::BINARY_WRITE_ONLY);
      ioRho.writeMesh (rho(i), cell, mesh);
      ioRho.setMode (SxBinIO::WRITE_DATA);
      ioRho.writeMesh (rho(i), cell, mesh);
      ioRho.close();

      // --- open output files
      ofstream file, wfile, file1d, outParFile, sum1d, vChgFile, vChg1d;
      fstream::openmode mode = fstream::out | fstream::trunc;
      file.open (("rho_" + SxString(i) + ".dat").ascii (), mode);
      vChgFile.open (("vChg.dat"), mode);
      vChg1d.open (("vChg1d.dat"), mode);
      wfile.open(("psi_" + SxString(i) + ".dat").ascii (), mode);
      file1d.open (("rho1d_" + SxString(i) + ".dat").ascii (), mode);
      sum1d.open (("sum1d_" + SxString(i) + ".dat").ascii (), mode);
      if (outPar != "") outParFile.open ((outPar + ".dat").ascii (), mode);

      int x,y,z,idx;
      double val, val2, val0, valChg;
      val0=0;z=0;
      SxVector3<int> pos;
      SxString waveStr;
      cout << "size of chargePotential: " << chargePotential.getSize() << endl;
      if ((chargePotential.getSize() > 1)) { cout << "write charge potential" << endl;}
      else {cout << "no charge potential" << endl;}
      for ( x = 0; x < mesh(0); x++)  {
         pos(0) = x;
         for ( y = 0; y < mesh(1); y++)  {
            pos(1) = y;
            for ( z = 0; z < mesh(2); z++)  {
               pos(2) = z;
               idx = (int)mesh.getMeshIdx(pos, SxMesh3D::Positive);
               val = rho(i)(idx);
               if (z == 0) val0 = val;
               file << SxString(val) << endl;
               // --- DEBUGGING ONLY, REMOVE LATER
               if ((chargePotential.getSize() > 1)) {
                  valChg = chargePotential(idx) * HA2EV;
                  vChgFile << SxString(valChg) << endl;
                  if ((y == mesh(0) - x) && (z == 1))  vChg1d << SxString(valChg) << endl;
               }
               // --- 'til 'ere
               waveStr = "";
               for (iComp = 0; iComp < nComp; iComp++)
                  waveStr = waveStr
                          + SxString((rComp(iComp)(idx)).re)
                          + "   "
                          + SxString((rComp(iComp)(idx)).im)
                          + "   ";
               wfile << waveStr << endl;
//               if ((x == mesh(0)/2) && (y == mesh(1)/2))  {
                  file1d << SxString(z) << "\t" << SxString(val) << endl;
                  val2 = outMesh(idx);
                  if (outPar != "") {
                     outParFile << SxString(val2) << endl;
                  }

//               }
            }
         }
      }
      file1d << SxString(z) << "\t" << SxString(val0) << endl;
      for ( z = 0; z < mesh(2); z++)  {
         pos(2) = z;
         val = 0.;
         for ( x = 0; x < mesh(0); x++)  {
            pos(0) = x;
            for ( y = 0; y < mesh(1); y++)  {
               pos(1) = y;
               idx = (int)mesh.getMeshIdx(pos, SxMesh3D::Positive);
               val = val + rho(i)(idx);
            }
         }
         sum1d << SxString(z) << "\t" << SxString(val) << endl;
      }
      cout << endl;
      // --- close files
      file.close ();
      vChgFile.close ();
      vChg1d.close ();
      wfile.close ();
      file1d.close ();
      sum1d.close ();
      if (outPar != "") outParFile.close ();
   }
   cout << endl;
   cout.precision (oldPrec);
}

SxRho &SxKdotP::getRho ()
{
   return *this;
}

void SxKdotP::computeRho (const SxPsiSet &wavesIn, const SxFermi &fermi)
{  
   SX_CHECK (dynamic_cast<const SxPW *> (&wavesIn));
   const SxPW &wavesRef = *dynamic_cast<const SxPW *> (&wavesIn);
   SxRho::computeRho (fermi.focc, wavesRef);
}

PrecEnergy SxKdotP::getEnergy (const SxPsiSet &psiSet,
                                  const SxFermi &)
{
   SX_CHECK (dynamic_cast<const SxPWSet *>(&psiSet));   
   wavesPtr = static_cast<const SxPWSet *> (&psiSet);
   
   return 0.; 
}


void SxKdotP::set (const SxPsiSet &psiSet, const SxFermi &)
{
   SX_CHECK (dynamic_cast<const SxPWSet *> (&psiSet));
   wavesPtr = dynamic_cast<const SxPWSet *>(&psiSet);
}

void SxKdotP::readRho (const SxBinIO &io)
{
   SxRho::readRho (io);
}

void SxKdotP::normalizeRho ()
{
   SxRho::normalizeRho ();
}

void SxKdotP::writeRho (const SxString &filename) const
{
   SxRho::writeRho (filename);
}

SxVector<double> 
SxKdotP::preconditioner (const SxVecRef<PrecCoeffG> &psi,
                                 Preconditioner ) const
{
   SxVector<double> x, x2, x3, x4, n, K;
   x.resize(psi.getSize());
   double kin = dot(D.real (),psi.absSqr());
   x = D / kin;
   x2 = x.sqr();
   x3 = x.cub();
   x4 = x2.sqr();
   n  = 27. + 18.*x + 12.*x2 + 8.*x3;
   K  = n / (n + 16.*x4);

   return K; 
}
