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

#include <SxTimeCorr.h>
#include <SxFFT1d.h>
#include <SxFFT3d.h>
#include <SxParser.h>


//----------------------------------------------------------------------------
//    Constructors and Destructors
//----------------------------------------------------------------------------


SxTimeCorr::SxTimeCorr ()
   : needsUpdateRC (false), 
     needsUpdateRN (false),
     needsUpdateNF (false),
     needsUpdateDM (false),
     needsUpdateCM (false),
     needsUpdateModes (false),
     needsUpdateFreqs (false)
{
   trajectoryL.resize (0);
   tL.resize (0);
   nSteps = 0;
   //TODO: clean up emty Constructor in SxFFT
   fft = SxFFT1d (SxFFT1d::Both, 1, 1);
   // scheme is hardcoded here, should be user controlled if necessary
   //scheme = SxString("realSpace");
   scheme = SxString("fft");
}

SxTimeCorr::~SxTimeCorr ()
{
   // empty
}

void SxTimeCorr::setType (SxString inType)
{
  type = inType;
}

void SxTimeCorr::setObservationRange (int start, int end) 
{
   if ((start != startIndex)||(end != endIndex)) {
      needsUpdateRC = needsUpdateRN = needsUpdateNF 
         = needsUpdateDM = needsUpdateCM = needsUpdateFreqs 
         = needsUpdateModes = true;
   }
   startIndex = start;
   endIndex = end;
   length = endIndex - startIndex + 1;
}

void SxTimeCorr::push (const SxAtomicStructure &X, double time_)
{
   SxAtomicStructure toAppend;
   toAppend.copy (X);
   trajectoryL.append (toAppend);
   tL.append (time_);
   nSteps += 1;
   needsUpdateRC = needsUpdateRN = needsUpdateNF 
   = needsUpdateDM = needsUpdateCM = true;
   nDoF = X.getNAtoms () * 3;
}

void SxTimeCorr::loadTrajectory 
(const SxString &filename, const SxString &quantity) 
{
   //--- reads in whole file before pushing to trajectory (should be refined)
   SxElemDB elemDB;
   
   ssize_t i, j;
   ssize_t nAtoms = 0;
   
   ifstream inStream;
   string line;
   int fileLength;
   //char *buffer;
   double t;
   SxAtomicStructure X;
   SxVector3<double> coord;

   SxList<SxString> xStringTok, coords;
   SxString step, timeString, xString;
   SxList<SxString>::Iterator itSteps;
   SxList<SxString>::Iterator itX;
   
   setType (quantity);
   inStream.open (filename.ascii ());
   inStream.seekg (0,ios::end);
   fileLength = int(inStream.tellg ());
   
   SxArray<char> buffer(fileLength);
   inStream.seekg (0, ios::beg);
   inStream.read (buffer.elements, fileLength);
   
   SxString whole = SxString(buffer.elements);
   
   //--- tokenizes into steps
   whole += SxString("\n@");
   SxList <SxString> stepList = whole.tokenize('@');

   itSteps = stepList.begin ();
   
   for (i = 0; i < stepList.getSize(); i++) {
      step = (*itSteps);
      timeString = (step.right ("time#")).left ("\n#");
      t = timeString.toDouble ();
      step = step + SxString ("#\n");
      
      xString = step.right(SxString(quantity + SxString("#"))).left("\n#");
      
      xStringTok = xString.tokenize ('\n');
      itX = xStringTok.begin ();
      
      if (i == 0) {
         nAtoms = xStringTok.getSize ();
         masses.resize (nAtoms);
      }
      
      //X.resize (nAtoms);
      X = SxAtomicStructure ();
      X.startCreation ();
      
      for (j = 0; j < nAtoms; j++) {
         coords = (*itX).tokenize(' ');
         
         if (i == 0) {
            masses (j) = elemDB.getAtomicWeight(coords(0));
         }
         
         coord(0) = coords(1).toDouble();
         coord(1) = coords(2).toDouble();
         coord(2) = coords(3).toDouble();
         X.addAtom (coord);
         itX++;
      }
      X.endCreation ();
      itSteps++;
      push (X, t);
   }
   updateTrajectory ();
}


   SxArray<double> SxTimeCorr::getPairCorrelationFunctionAtom 
(int iAtom, int jAtom) 
{
   
   int i, k;
   
   updateTrajectory ();
   
   SxArray<double> returnValue (length);
   SxArray<double> toAdd (length);
   
   for (i = 0; i < length; i++) returnValue (i) = 0.;
   for (k = 0; k < 3; k++) {
      toAdd = getPairCorrelationFunctionDoF
         (3*iAtom + k, 3*jAtom + k);
      
      for (i = 0; i < length; i++) returnValue (i) += toAdd(i);
   }
   return returnValue;
}

   SxArray<double> SxTimeCorr::getPairCorrelationFunctionDoF 
                  (int iAtom, int jAtom) 
{
   
   int tp, t, index, i;
   
   updateTrajectory ();
   
   SxArray<double> returnValue (length);
   if (scheme == SxString("realSpace")) {
      for (t = startIndex; t < endIndex; t++) {
         returnValue(t - startIndex) = 0.;
      }
      for (tp = 0; tp < length; tp++) {
         returnValue(tp) = 0.;
         for (t = startIndex; t < endIndex; t++) {
            index = t + tp;
            if (index > endIndex) index -= length;
            returnValue(tp)
               += (trajectory(t).coordRef()(iAtom)
                     * trajectory(index).coordRef()(jAtom));
         }
         returnValue(tp) /= (double)length;
      }
      return returnValue;
   } 
   
   if (scheme == SxString("fft")) {
      renewFFT (length, 1.);
      SxVector<SxComplex16> fftRI (length);
      SxVector<SxComplex16> fftRJ (length);
      SxVector<SxComplex16> fftGI (length);
      SxVector<SxComplex16> fftGJ (length);
      SxVector<SxComplex16> fftRErg (length);
      SxVector<SxComplex16> fftGErg (length);
      
      fftRI.set (0.);  fftRJ.set (0.);  fftGI.set (0.); 
      fftGJ.set (0.);  fftRErg.set (0.);  fftGErg.set (0.); 
      
        for (t = startIndex; t <= endIndex; t++) {
           fftRI(t - startIndex).re = trajectory(t).coordRef()(iAtom);
           fftRJ(t - startIndex).re = trajectory(t).coordRef()(jAtom);
           fftRJ(t - startIndex).im = fftRI(t - startIndex).im = 0.;
        }
         fft.fftReverse (length, fftRI.elements, fftGI.elements);
         fft.fftReverse (length, fftRJ.elements, fftGJ.elements);
      
         fftGErg = fftGJ*fftGI.conj ();
      
      
      fft.fftForward (length, fftGErg.elements, fftRErg.elements);

      for (i = 0; i < length; i++) {
         returnValue (i) = fftRErg(i).re;
      }
      return returnValue;

   }
   cout << "SxTimeCorr: Unknown convolution scheme: " << scheme << endl;
   SX_EXIT;
   return returnValue;
} 

double SxTimeCorr::getPairCovarianceDoF 
                  (int iDoF, int jDoF) 
{

   int  t;
   updateTrajectory ();
   double returnValue = 0.;
   for (t = 0; t < length; t++) {
      
      returnValue
         += (trajectory(t + startIndex).coordRef()(iDoF)
               * trajectory(t + startIndex).coordRef()(jDoF));
   }
   returnValue /= (double)length;
   return returnValue;
} 


SxArray<double> SxTimeCorr::getAutoCorrelationFunction (int iAtom) 
{
   return getPairCorrelationFunctionAtom(iAtom, iAtom);
}

SxArray<double> SxTimeCorr::getAveragedAutoCorrelationFunction () 
{
   int iAtom, j;
   int nAtoms = trajectory(0).getNAtoms ();
   SxArray<double> summand;
   SxArray<double> returnValue(length);
   
   for (j = 0; j < length; j++) {
      returnValue(j) = 0.;
   }
   for (iAtom = 0; iAtom < nAtoms; iAtom++) {
      summand = getAutoCorrelationFunction (iAtom);
      for (j = 0; j < length; j++) {
         returnValue(j) += summand(j)/(double)nAtoms;
      }
   }
   return returnValue;
}

SxArray<double> SxTimeCorr::getGeneralizedFrequencySpectrum () 
{
   int i;
   SxArray<double> returnValue(length);
   SxArray<double> autoCorr(length);
   SxVector<SxComplex16> fftIn (length);
   SxVector<SxComplex16> fftOut (length);
   
   updateTrajectory ();
   autoCorr = getAveragedAutoCorrelationFunction ();
   renewFFT (length, 1.);
   for (i = 0; i < length; i++) {
      fftIn(i).re = autoCorr(i);
      fftIn(i).im = 0.;
   }
   fft.fftReverse (length, fftIn.elements, fftOut.elements);
   for (i = 0; i < length; i++) {
      returnValue(i) = fftOut(i).re;
   }
   return returnValue;
}

double SxTimeCorr::getDFreq () 
{
   double dfreq = 1./((double)length*dt*AU2S)/(CVEL*100.);
   return dfreq;
}

int SxTimeCorr::getIndex (double time_) 
{
   int index = 0;
   if ((time_ < startTime) || (time_ > endTime)) {
      cout << "Error in SxTimeCorr: getting outside recording "
           << "time window" << endl;
      SX_QUIT;
   } else {
      index = (int) ((time_ - startTime)/dt);
   }
   return index;
}

void SxTimeCorr::updateTrajectory () 
{
   SxList<SxAtomicStructure>::Iterator itTraj;
   SxList<double>::Iterator itT;
   int i;

   if (needsUpdateRC) {
      trajectory.resize (trajectoryL.getSize ());
      time.resize (tL.getSize ());
      itTraj = trajectoryL.begin ();
      itT = tL.begin ();
      for (i = 0; i < nSteps; i++) {
         trajectory(i).copy (*itTraj);
         time(i) = *itT;
         itTraj++;
         itT++;
      }
      needsUpdateRC = false;
   }
   startTime = time(0);
   endTime = time(nSteps - 1);
   dt = time(1) - time(0);

}
      
void SxTimeCorr::renewFFT (int size, double omega)
{
  if (size != fft.meshSize) {
     cout << "Renewing fft: " << size << endl;
    fft.setMesh (size, omega, 1);
    
     cout << " fft: " << fft.meshSize << endl;
  }
}

SxVector<double> SxTimeCorr::getCrossCorrelationMatrix ()
{
   int nDoF_ = trajectory(0).getNAtoms () * 3;
   int i, j, k;
   
   SxArray<double> cf (length);
   SxVector<double> returnValue (nDoF_, nDoF_);
   returnValue.set (0.);
   
   for (i = 0; i < nDoF_; i++) {
      for (j = 0; j < nDoF_; j++) {
         cf = getPairCorrelationFunctionDoF (i, j);
         for (k = startIndex; k <= endIndex; k++) {
           returnValue(i, j) += cf(k)/((double)length*dt);
         } 
      }
   }
   return returnValue;
}
   
const SxVector<double>& SxTimeCorr::getCovarianceMatrix ()
{
   int i, j;
   
   if (needsUpdateCM) {
      cout << "Evaluating covariance Matrix ...";
      fflush(stdout);
      covarianceMatrix = SxVector<double> (nDoF, nDoF);
      for (i = 0; i < nDoF; i++) {
         for (j = 0; j < nDoF; j++) {
            covarianceMatrix(i, j) = getPairCovarianceDoF (i, j);
         }
      }
      needsUpdateCM = false;
      cout << " ready !" << endl;
   }
  
   return covarianceMatrix;
}

const SxVector<double>& SxTimeCorr::getNormalModes ()
{
   int i, j;

   const SxVector<double> &covarianceMatrix_ = getCovarianceMatrix ();
   SxVector<double> massWeighted(nDoF, nDoF);
   for (i = 0; i < nDoF; i++) {
      for (j = 0; j < nDoF; j++) {
         massWeighted(i, j) 
            = covarianceMatrix_(i, j)*sqrt(masses(i/3)*masses(j/3));
      }
   }
   
   if (needsUpdateModes) {
      cout << "Evaluating normal Modes ..."; fflush(stdout);
      SxSymEigensystem<double> eig(massWeighted);
      normalModes = std::move(eig.vecs);
      cout << " ready !" << endl;
      needsUpdateModes = false;
   }
   return normalModes;
}

double SxTimeCorr::checkHermiticity (const SxVector<SxComplex16> &m)
{
   ssize_t i, j, size;
   size = (m.getNRows ());
   double returnValue = 0.;

   for (i = 0; i < size; i++) {
      for (j = 0; j <= i; j++) {
         returnValue += (m(i, j) - m(j, i).conj()).abs ();
      }
   }
   return returnValue;
}

void SxTimeCorr::sortByEigenvectors (const SxSymEigensystem<double> &eig1,
                                           SxSymEigensystem<double> *eig2)
{
   int i, j;
   SxVector<double> newVecs (nDoF, nDoF);
   SxVector<double> newVals (nDoF);

   newVecs.set (0.);
   newVals.set (0.);

   double maxSimilarity = 1.e15;
   double similarity = 0.;
   int maxSimilarityIndex = 0;

   //(*eig1).vecs.orthogonalize ();
   //(*eig2).vecs.orthogonalize ();

   for (i = 0; i < nDoF; i++) {
      const SxVecRef<double> &aVec = (*eig2).vecs.colRef (i);
      maxSimilarity = 0.; 
      for (j = 0; j < nDoF; j++) {
         cout << j << endl;
         const SxVecRef<double> &bVec = eig1.vecs.colRef (j);
         similarity = fabs(dot(aVec, bVec));
         if (similarity > maxSimilarity  && fabs(newVals(j)) < 1e-15)  {
            maxSimilarity = similarity;
            maxSimilarityIndex = j;
            cout << similarity << endl;
         }
      }
      cout << "-----------------" << i << endl;
      if (fabs(newVals(maxSimilarityIndex)) > 1e-15) {
         cout << "Error: sorting Eigensystems: more then one vektor of system 2"
           << "are very similar to one vektor of system 1" << endl;
         cout << "vector system 1: " << eig1.vecs.colRef(maxSimilarityIndex)
              << endl;
         cout << "vector A system 2 (value=" << newVals(maxSimilarityIndex)
              << "): " << newVecs.colRef(maxSimilarityIndex) << endl;
         cout << "overlap="
              << fabs(dot(aVec,newVecs.colRef(maxSimilarityIndex))) << endl;
         cout << "vector B system 2 (value=" << eig2->vals(i)
              << "): " << aVec << endl;
         cout << "overlap=" << maxSimilarity << endl;
         SX_EXIT;
      } else {
         newVals(maxSimilarityIndex) = (*eig2).vals(i);
         newVecs.colRef(maxSimilarityIndex) = (*eig2).vecs.colRef(i);
      }
   }

   (*eig2).vecs = std::move(newVecs);
   (*eig2).vals = std::move(newVals);
}

SxVector<double> SxTimeCorr::getHessian (SxTimeCorr &tcVel, SxTimeCorr &tcForce)
{
   SxVector<double> covarianceVel, covarianceForce, dynamicalStruc;
   int i, j;

   covarianceVel = tcVel.getCovarianceMatrix ();
   for (i = 0; i < nDoF; i++) {
      for (j = 0; j < nDoF; j++) {
         covarianceVel(i, j) = covarianceVel(i, j)
                            *sqrt(tcVel.masses(i/3)*tcVel.masses(j/3));
      }
   }
   SxSymEigensystem<double> eigVel(covarianceVel);

   covarianceForce = tcForce.getCovarianceMatrix ();
   for (i = 0; i < nDoF; i++) {
      for (j = 0; j < nDoF; j++) {
         covarianceForce(i, j) = covarianceForce(i, j)
                            /sqrt(tcForce.masses(i/3)*tcForce.masses(j/3));
      }
   }

   SxSymEigensystem<double> eigAcc (covarianceForce);
   tcForce.sortByEigenvectors(eigVel, &eigAcc);

   SxVector<double> dynamical;
   SxVector<double> hessian(nDoF, nDoF);
   SxVector<double> diag(nDoF, nDoF);
   diag.set(0.);
   for (i = 0; i < nDoF; i++) {
      diag (i, i) = eigAcc.vals(i)/eigVel.vals(i);
      if (fabs(eigAcc.vals(i)) < 1e-8) diag(i, i) = 0.;
   }
   dynamical = eigVel.vecs ^ diag ^ eigVel.vecs.transpose ();

   for (i = 0; i < nDoF; i++) {
      for (j = 0; j < nDoF; j++) {
         hessian(i, j) = sqrt(tcForce.masses(i/3)*tcForce.masses(j/3))
                       * dynamical(i, j);
      }
   }

   return hessian;

}

void SxTimeCorr::setMasses (const SxVector<double> &massesIn) {
int i;
masses.resize (nDoF/3);
if (massesIn.getSize () == nDoF) {
      for (i = 0; i < nDoF/3; i++) {
         masses(i) = massesIn(i*3);
      }
      return;
   }
   if (massesIn.getSize () == nDoF/3) {
      for (i = 0; i < nDoF/3; i++) {
         masses(i) = massesIn(i);
      }
      return;
   }
   cout << "SxTimeCorr: Input Mass Vector has wrong size: "
      << massesIn.getSize () 
      << " !" << endl;
   SX_QUIT;
}


//--- non-tested routines (project under development)
//    L. Ismer (03/17/05)
/*
void SxTimeCorr::updateTrajectoryNF () 
{
   if (needsUpdateNF) {
      int i, j;
      SxVector<SxComplex16> fftIn (length);
      SxVector<SxComplex16> fftOut (length);
      renewFFT (length, 1.);
      trajectoryNF.resize(0);
      trajectoryNF.resize(length);
      for (j = 0; j < length; j++) {
         trajectoryNF(j) = SxVector<SxComplex16> (nDoF);
      }
      updateTrajectoryN ();
      cout << "Evaluating fourier analysis of trajectory "
           << "in normal coordinates ...";
      for (i = 0; i < nDoF; i++) {
         for (j = 0; j < length; j++) {
            fftIn(j).re = trajectoryN(j).coordRef()(i);
            fftIn(j).im = 0.;
         }
         fft.fftReverse (length, fftIn.elements, fftOut.elements);
         for (j = 0; j < length; j++) {
            trajectoryNF(j)(i) =  fftOut(j);
         }
      }
      cout << " ready !" << endl;
      needsUpdateNF = false;
   }
}


void SxTimeCorr::updateTrajectoryN () 
{
   SxVector<double> modes = getNormalModes ();
   SxVector<double> projection (nDoF);

   updateTrajectory ();
   cout << "Evaluating trajectory in normal coordinates ...";
   trajectoryN.resize (0);
   trajectoryN.resize (length);
   for (int i = startIndex; i <= endIndex; i++) {
      projection = modes^(trajectory(i).coordRef ());
      trajectoryN(i - startIndex).set(projection);
   }
   cout << " ready !" << endl;
} 


SxVector<double> SxTimeCorr::getDynamicalMatrix ()
{
   SxVector<double> modes = getNormalModes ();
   SxVector<double> id (nDoF, nDoF);
   SxVector<SxComplex16> freqs = getNormalFrequencies ();

   id.set (0.);
   for (int i = 0; i < freqs.getSize (); i++) id(i) = freqs(i)*freqs(i);

   cout << "Evaluating dynamical Matrix ...";
   dynamicalMatrix = (modes.transpose ()^id^modes);
   cout << " ready !" << endl;
   return dynamicalMatrix;
}

SxArray<double> SxTimeCorr::getNFProjection(int i)  
{
   int j;
   updateTrajectoryNF ();
   SxArray<double> returnValue(length);
   for (j = 0; j < length; j++) {
      //returnValue(j) = trajectoryNF(j)(i);
      SxComplex16 a = trajectoryNF(j)(i);
      returnValue(j) = a;
   }
   return returnValue;
}
      
SxVector<SxComplex16> SxTimeCorr::getNormalFrequencies ()
{
   updateTrajectoryNF ();
   normalFrequencies.resize (nDoF);
   double dFreq = getDFreq ();
   SxComplex16 norm;
   SxComplex16 average;

   SxVector<SxComplex16> fftIn (length);
   SxVector<SxComplex16> fftOut (length);
   renewFFT (length, 1.);
   cout << "Evaluating normal Frequencies ...";
   for (int i = 0; i < nDoF; i++) {
      average = 0.;
      norm = 0.;
      for (int j = 0; j < length; j++) {
         norm +=  trajectoryNF(j)(i);
         average += trajectoryNF(j)(i)*(double)j;
      }
      average *= dFreq/norm;
      normalFrequencies(i) = average;
   }
   cout << " ready !" << endl;
   return normalFrequencies;
}
            
*/ 
