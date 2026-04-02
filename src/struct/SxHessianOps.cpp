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

#include <SxHessianOps.h>
#include <SxTextIO.h>

//----------------------------------------------------------------------------
//    Constructors and Destructors
//----------------------------------------------------------------------------

SxHessianOps::SxHessianOps ()
{
   //empty
} 

SxHessianOps::SxHessianOps 
(  const SxHessianOps &in)
{
   nDoF = in.getNDoF ();
   hessian = in.hessian;
   dynamical = in.dynamical;
   massVec = in.massVec;
   eigDyn   = in.eigDyn;
   eigHessian = in.eigHessian;
}

SxHessianOps::SxHessianOps 
(  const SxVecRef<SxComplex16> &hessianIn, const SxVecRef<double> &massVecIn)
{
   int i, j;
   nDoF = (int)massVecIn.getSize ();
   setMasses (massVecIn);
   setHessian (hessianIn);
   dynamical = SxVector<SxComplex16>(nDoF, nDoF);
   for (i = 0; i < nDoF; i++) {
      for (j = 0; j < nDoF; j++)  
         dynamical(i, j) = hessian(i, j)/sqrt (massVec(i)*massVec(j));
   }
   update ();
}

SxHessianOps::~SxHessianOps ()
{
   // empty
}

void SxHessianOps::set (const SxVecRef<SxComplex16> &inH, const SxVecRef<double> &inM)
{
   int i, j;
   nDoF = (int)inM.getSize ();
   setMasses (inM);
   setHessian (inH);
   dynamical = SxVector<SxComplex16>(nDoF, nDoF);
   for (i = 0; i < nDoF; i++) {
      for (j = 0; j < nDoF; j++)  
         dynamical(i, j) = hessian(i, j)/sqrt (massVec(i)*massVec(j));
   }
   update ();
}

void SxHessianOps::setFromRefinement 
(const SxHessianOps &basisHOps, const SxVecRef<double> &responsesIn)
{
   nDoF = basisHOps.getNDoF ();
   const SxVector<double> &T = basisHOps.eigHessian.vecs;
  
   setMasses (basisHOps.massVec);
   
   //--- projection on displacement vectors
   SxVector<double> responses = T.transpose () ^ responsesIn;
 
   //--- transformation to cartesian basis
   setHessian (T ^ (responses ^ T.transpose ()));
   dynamical = SxVector<SxComplex16> (nDoF, nDoF);
   for (int i = 0; i < nDoF; i++)
      for (int j = 0; j < nDoF; j++)
         dynamical(i, j) = hessian(i, j)/sqrt (massVec(i)*massVec(j));
   update ();
}


SxHessianOps &SxHessianOps::operator= (const SxHessianOps &in)
{
   if ( &in == this )  return *this;

   nDoF = in.getNDoF ();
   hessian = SxVector<SxComplex16> (nDoF, nDoF); hessian.copy (in.hessian);
   dynamical = SxVector<SxComplex16> (nDoF, nDoF); dynamical.copy (in.dynamical);
   massVec = SxVector<double> (nDoF); massVec.copy (in.massVec);
   
   //--- using the constructors as shown in the below two lines
   //   doesn't work properly ... therefore commented out 
   //   however, updating by re-calculating the eigensystems 
   //   is a bit inefficient and should be replaced soon 
   
  // eigDyn   = SxVector<double>::Eigensystem (in.eigDyn);
  // eigHessian 
  //    = SxVector<double>::Eigensystem (in.eigHessian);
   
   
   needsUpdate = true;
   update ();
   return *this;
}

//----------------------------------------------------------------------------
//    Member Functions
//----------------------------------------------------------------------------

SxVector<SxComplex16> SxHessianOps::getHessian () {
   return hessian;
}

SxVector<SxComplex16> SxHessianOps::getDynamical () {
   return dynamical;
}

int SxHessianOps::getNDoF () const
{
   return nDoF;
}

SxVector<double> SxHessianOps::getEigenVelocities ()
{
   int i, j;
   update ();
   
   SxVector<double> invSqrtM = 1./sqrt(massVec);
   SxVector<double> returnValue (nDoF, nDoF);
   for (i = 0; i < nDoF; i++)
      for (j = 0; j < nDoF; j++)
         returnValue(j,i) = eigDyn.vecs(j,i) * invSqrtM(j);

   return returnValue;
}

SxVector<SxComplex16> SxHessianOps::getEigenFrequencies ()
{
   int i;
   double eV;
   SxVector<SxComplex16> returnValue (nDoF);
   
   update ();
   returnValue.set (0.);
   
   for (i = 0; i < nDoF; i++) {
      eV = eigDyn.vals(i);
      if (eV >= 0.)
         returnValue(i).re = sqrt(eV)*5123.75;
      else 
         returnValue(i).im = sqrt(-eV)*5123.75;
   } 
   return returnValue;
}

SxVector<double> SxHessianOps::getReducedMasses ()
{
   SxVector<double> returnValue (nDoF);
   update ();
   for (int i = 0; i < nDoF; i++) {
      double norm = 0.;
      for (int j = 0; j < nDoF; j++)  {
         double c = eigDyn.vecs(j,i);
         norm += c*c/massVec(j);
      }
      returnValue(i) = 1./norm;
   }
   return returnValue;
}

void SxHessianOps::setMasses (const SxVecRef<double> &massVecIn)
{
   needsUpdate = true;
   massVec = SxVector<double> (nDoF);
   massVec.copy (massVecIn);
}

void SxHessianOps::setHessian (const SxVecRef<SxComplex16> &hessianIn)
{
   needsUpdate = true;
   hessian = SxVector<SxComplex16> (nDoF, nDoF);
   hessian.copy (hessianIn);
}


void SxHessianOps::update () 
{
   SxVector<double> a (nDoF, nDoF);
   SxVector<double> b (nDoF, nDoF);
   if (needsUpdate) {
      eigDyn = SxSymEigensystem<SxComplex16>(dynamical);
      eigDyn.vecs = getOrthonormalizedMatrix(eigDyn.vecs);
      eigHessian = SxSymEigensystem<SxComplex16>(hessian);
      eigHessian.vecs = getOrthonormalizedMatrix(eigHessian.vecs);
      
      needsUpdate = false;
   }
}
      
void SxHessianOps::setCurvaturesDynBasis 
                     (SxList<SxList<double> > &toReplace) 
{
   int i, j, index;
   SxVector<double> diagonal (nDoF, nDoF);
   SxVector<double> reducedMasses = getReducedMasses ();
   cout << reducedMasses << endl;
   SxVector<double> column(nDoF); SxVector<double> orth(nDoF);
   diagonal.set (0.);
   for (i = 0; i < nDoF; i++)  diagonal(i, i) = eigDyn.vals (i);
   needsUpdate = true;

   //--- to avoid inconsitencies a symmetrization of the 
   //    new curvatures according to the degenracy of 
   //    the old curvatures is necessary 
   bool changed = true;
   double k1, k2;
   while (changed) { 
      changed = false;
      for (i = 0; i < (toReplace.getSize () - 1); i++) {
         index = (int)(toReplace(i)(0) - 1);
         k1 = toReplace(i)(1);
         k2 = toReplace(i + 1)(1);
         if (   (fabs (diagonal (index, index) 
                    - diagonal(index + 1, index + 1)) < 1e-7) 
             && (fabs (k1 - k2) > 1e-7) ) {
            changed = true;  
            k1 = (k1 + k2)/2.;
            toReplace(i)(1) = k1;
            toReplace(i + 1)(1) = k1;
         }
      }
   }
   

   for (i = 0; i < toReplace.getSize (); i++) {
      index = (int)(toReplace(i)(0) - 1);
      if ((index < 0) || (index >= nDoF)) {
         cout << "Ill-defined index for replacing curvature: "
              << index << endl;
         cout << "Should be in the DoF-range : 1 ..." << nDoF << endl;
         SX_QUIT;
      }
      diagonal (index, index) = toReplace(i)(1);
   } 
     
   const SxVector<double> &BW = eigDyn.vecs;
   dynamical = BW ^ (diagonal ^ BW.inverse ());
   
   for (j = 0; j < nDoF; j++) {
      for (i = 0; i < nDoF; i++) {
         hessian(i, j) = sqrt(massVec(i)*massVec(j))*dynamical(i, j);
      }
   }

}

void SxHessianOps::setFrequencies 
                     (SxList<SxList<SxComplex16> > &toReplace) 
{
   int i, index;
   SxVector<double> reducedMasses = getReducedMasses ();
   SxList<SxList<double> > curvatures;
   SxList<double> listEntry;
   double k;
   SxComplex16 freq;

   for (i = 0; i < toReplace.getSize (); i++) {
      //rescale to curvatures consumerable by setCurvaturesDynBasis
      freq = toReplace(i)(1);
      index = (int)(toReplace(i)(0).re - 1);
      k = (freq.abs () * freq.abs ())
          /5123.75/5123.75;
      //k = k/reducedMasses(index)/reducedMasses(index);
      
      listEntry.resize (0);
      listEntry.append (index + 1);
      listEntry.append (k);
      curvatures.append (listEntry);
   } 
   setCurvaturesDynBasis (curvatures);
}

void SxHessianOps::write (const SxString &fileName)
{
   SxTextIO io;
   io.open (fileName, SxTextIO::WriteSx);
   io.writeMat ("matrix", hessian);

   SxVector<SxComplex16> vals = getEigenFrequencies ();
   io.writeVec ("frequencies (unit 1/cm)", vals);
   io.close();

   SxString fileNameRaw = fileName + SxString (".out");
   io.open (fileNameRaw, SxTextIO::WriteOnly);
   io.writeMat (hessian);
   io.close ();
}

SxVector<double> SxHessianOps::getDisplacementVector (int DoF)
{
   SxVector<double> returnValue = eigHessian.vecs.colRef (DoF);
   
   //--- displacement length is normed according to curvature of PES 
   double curv = fabs(eigHessian.vals(DoF));
   
   //--- threshold value (if you change this, please also change the according 
   //    value in getHessianFromRefinement)
   
   if (curv < 1e-9) curv = 1e-9;

   returnValue = returnValue/sqrt(curv);
   return returnValue;
}

SxVector<double> SxHessianOps::getSymmetrizedMatrix (const SxVecRef<double> &h)
{
   int i, j;
   int size = (int)h.getNRows ();
   double avg;
   SxVector<double> returnValue (size, size);
   for (i = 0; i < size; i++) {
      for (j = 0; j <= i; j++) {
         avg = (h(i, j) + h(j, i))/2.;
         returnValue(i, j) = returnValue(j, i) = avg;
      }
   }
   return returnValue;
}

SxVector<SxComplex16>
SxHessianOps::getOrthonormalizedMatrix (const SxVecRef<SxComplex16> &h)
{
   int i, j, k;
   bool found = false;
   int size = (int)h.getNRows ();
   SxVector<SxComplex16> returnValue (size, size);

   for (i = 0; i < size; i++) {
      SxVecRef<SxComplex16> col = returnValue.colRef(i);
      col <<= h.colRef(i);
      //--- for portability purposes (not all platforms return the samle signum)
      k = 0;
      found = false;
      while (!found) {
         if (fabs(col(k).re) > 1e-15) {
            if (col(k).re <= 0.)
               col *= -1.;
            found = true;
         } else k++;
      }
         
		
      // Gram-Schmidt orthogonalization
      for (j = 0; j < i; j++) {
         const SxVecRef<SxComplex16> &orth = returnValue.colRef(j);
         col.plus_assign_ax(-dot(orth,col), orth);
      }
      col.normalize ();
   }
   return returnValue;
}
         
void SxHessianOps::printMolden (const SxString &fileName, 
                               const SxAtomicStructure &tau,
                               int REPLICZ)
                               
{
   FILE *fp = NULL;
   int i,  zReplic, j;
   SxAtomicStructure tauLoc; tauLoc.copy (tau);
   double zExt = tau.cell(2)(2);
    
   SxString line, numberOne, numberTwo, numberThree, numberZero, toSave,
            innerLine;
   SxVector<double> help(nDoF);
   SxList<SxString> chemName;
   SxVector<SxComplex16> eigenFrequencies = getEigenFrequencies ();
   SxVector<double> eigenVelocities = getEigenVelocities ();
   
   //--- parsing in chemical elements
   SxElemDB elemDB;
   
   for (i = 0; i < nDoF/3; i++) 
      chemName.append (elemDB.getChemSymbol (massVec(i*3)));
   
   line = "[Molden Format]\n";
   
   line += "[FREQ]\n";
   for (i = 0; i < nDoF; i++) {
      if (eigenFrequencies(i).re > eigenFrequencies(i).im)
         line += SxString ((float) eigenFrequencies(i).re);
      else {
         line += SxString ((float) eigenFrequencies(i).im);
         line += "*I";
      }
      line += "\n";
   }
   
   line += "[FR-COORD]\n";

   // --- the input-variable zReplic contains how often the structure should
   //     be replicated in the z-Direction, this is a very special option
   //     for being able to vizualize infinite helical molecules 

   for (zReplic = 1; zReplic <= REPLICZ; zReplic++) {  
      for (i = 0; i < nDoF/3; i++) {
            // --- has to be replaced by something like element(is) 
            //     (chemical formula)
            numberZero = chemName(i);
            numberOne = SxString(tauLoc.ref(i)(0));
            numberTwo = SxString(tauLoc.ref(i)(1));
            numberThree = SxString(tauLoc.ref(i)(2) 
                  + (double)(zReplic - 1.)*zExt) ;
            // NgString numberThree (tau(is)(ia)(2)); 
            
            toSave = 
               numberZero + SxString(" ") + numberOne + SxString(" ")
               + numberTwo + SxString(" ")
               + numberThree + SxString("\n");
            line += toSave;
         }
   } 
   
   
   line += "[FR-NORM-COORD]\n";
   
   for (i = 0; i < nDoF; i++) {
      line += SxString("dummy ") + SxString(i) + SxString("\n");
      help <<= eigenVelocities.colRef(i);
     
      innerLine = SxString();
      for (zReplic = 1; zReplic <= REPLICZ; zReplic++) {  
         for (j = 0; j < nDoF; j++) {
            innerLine += SxString(help(j)); 
            innerLine += SxString(" ");
            j++;
            innerLine += SxString(help(j)) + SxString(" ");
            j++;
            innerLine += SxString(help(j)) + SxString("\n");
         }
      }
      line += innerLine;
   }
   
   if ( !(fp = fopen (fileName.ascii (), "w")) ) {
      sxprintf ("Can't open file %s", fileName.ascii ()); 
      SX_EXIT; 
   }     
   fprintf (fp, "%s\n", line.ascii ());
   
   fclose (fp);
}
 
