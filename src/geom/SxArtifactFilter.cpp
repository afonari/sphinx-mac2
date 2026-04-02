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
//   AUTHOR: Lars Ismer, ismer@fhi-berlin.mpg.de

#include <SxArtifactFilter.h>
#include <SxEigensystem.h>

SxArtifactFilter::SxArtifactFilter ()
{
     // empty
}


void SxArtifactFilter::set (const SxAtomicStructure &tauIn,
                                    const SxString &freezeModeIn,
                                    bool massWeightingIn)
{
   needsUpdateDM = true;
   setStructure (tauIn);
   setFreezeMode (freezeModeIn);
   setMassWeighting (massWeightingIn);
   DMatrix = getDMatrix ();
   // empty
}

SxArtifactFilter::~SxArtifactFilter ()
{
   // empty
}

SxVecRef<PrecTauR>
SxArtifactFilter::operator*(const SxVecRef<PrecTauR> &in) const
{
  if (freezeMode != SxString ("notr")) {
     SxVector<PrecTauR> DDT = DMatrix ^ DMatrix.transpose ();
     if ((in.getNRows () > 1) && (in.getNCols () > 1)) {
        return DDT ^ in ^ DDT;
     } else {
        return DDT ^ in;
     }
  } else {
     //  cout << "Artifact Filter has no effect !" << endl;
     return const_cast<SxVecRef<PrecTauR>&>(in);
  }
}

void SxArtifactFilter::setFreezeMode (const SxString &freezeModeIn)
{
   freezeMode = freezeModeIn;
   //needsUpdateDM = true;
}

void SxArtifactFilter::setStructure (const SxAtomicStructure &structureIn)
{
   tau.copy (structureIn);
   //needsUpdateDM = true;
}

void SxArtifactFilter::setMassWeighting (bool isActive)
{
   massWeighting = isActive;
   //needsUpdateDM = true;
}

void SxArtifactFilter::optimizeDMatrix (const SxVecRef<double> &/*m*/)
{
   SX_EXIT;
   /*
   cout << "Optimizing D-Matrix ..." << endl;

   int freeRots = 0;
   if (freezeMode == SxString("z")) freeRots = 1;
   if (freezeMode == SxString("xyz")) freeRots = 3;
   if (freezeMode == SxString("trans"))  freeRots = 0;
   */

}



SxVector<double> SxArtifactFilter::getDMatrix ()
{
      int i, j, k;
      int freeRots = 0;
      int size = tau.getNAtoms ()*3;
      SxVector<double> tauVec = tau.coordRef ();
      SxVector3<double> COM;
      SxVector<double> inertiaTensor (3, 3);
      SxVector<double> D (size, size);

      const SxVecRef<double> &masses = getMasses ();

      if (freezeMode == SxString("z")) freeRots = 1;
      if (freezeMode == SxString("xyz")) freeRots = 3;
      if (freezeMode == SxString("trans"))  freeRots = 0;

      D.set (0.);

      if (freezeMode != SxString ("notr")) {

         // --- columns of the D-Matrix according to translational degrees
         //     of freedom
         for (j = 0; j < 3; j++) {
            for (int i3 = 0; i3 < size; i3+=3)
               D(i3 + j, j) = sqrt (masses(i3));
         }

         if (freezeMode == SxString ("xyz") || freezeMode == SxString ("z")) {

            // --- structure needs to be shifted to center of mass
            //     to calculate inertia tensor
            COM = translateToCOM(&tauVec);
            inertiaTensor = getInertiaTensor(tauVec);
            SxSymEigensystem<double> IEig (inertiaTensor);

            // --- if freezing only around the z-axis is used to generate
            //     rotational movements (see below) for xyz the eigenvectors
            //     of the inertia tensor are used

            if (freezeMode == SxString ("z")) {
               IEig.vecs.set (0.);
               IEig.vecs(0,2) = 1.;
            }

            // --- the columns of D according to rotational movements are
            //     generated
            for (k = 0; k < freeRots; k++) {
               const SxVector3<double> &IV
                  = SxVector3<double>::toVec3Ref(&IEig.vecs(0,k));
               for (i = 0; i < size/3; i++) {
                  const SxVector3<double> &R
                     = SxVector3<double>::toVec3Ref(&tauVec(i*3));

                  // set D(3*i .. 3*i+2, 3+k) = (R x IV) * sqrt(M_i)
                  SxVector3<double>::toVec3Ref(&D(i*3, 3+k))
                     = R.x (IV) * sqrt(masses(3*i));

               }
            }
         }

         // --- Normalisation of the translational and rotational columns
         //     of D
         for (i = 0; i < (3 + freeRots); i++)
            D.colRef(i).normalize ();

         //--- gram-schmidt orthonormalisation yields the other
         //    columns of D, which are then purely internal movements

         for (i = 1; i < (3 + freeRots); i++) {
            SxVecRef<double> column = D.colRef(i);
            for (j = 0; j < i; j++) {
               const SxVecRef<double> &orth = D.colRef(j);
               column.plus_assign_ax(-dot (column, orth), orth);
            }
            column.normalize ();
         }

         int shift = 1; // CF ??? -> should be 0, I think.
         for (i = 0; i < (size - (3 + freeRots)); i++) {
            SxVecRef<double> column = D.colRef(i + 3 + freeRots);
            // set D to (i+shift)-th unit vector
            column.set (0.);
            column(i + shift) = 1.;
            // --- orthogonalize to previous degrees of freedom
            for (j = 0; j < (3 + freeRots + i); j++) {
               const SxVecRef<double> &orth = D.colRef(j);
               column.plus_assign_ax (-dot(column, orth), orth);
               if (column.norm () < 0.1) {
                  SX_EXIT; // CF 2020-05-29 -> this would crash
                  i--;
                  shift++;
                  break;
               }
               column.normalize ();
            }
         }

         //--- the D matrix further used is build up by the
         //    the columns which belong to purely internal
         return D.blockRef (0,            3 + freeRots,   // offset
                            size, size - (3 + freeRots)); // sizes
   } else {
      SxVector<double> D_ (size, size - (3 + freeRots));
      D_.set(0.);
      return D_;
   }
}

const SxVecRef<double> SxArtifactFilter::getMasses () const
{
   if (massWeighting) return const_cast<SxVector<double>&>(massVec);

   SxVector<double> masses(tau.getNAtoms () * 3);
   masses.set (1.);
   return masses;
}

void SxArtifactFilter::setMasses (const SxSpeciesData &speciesData)
{
   massVec.resize (tau.getNAtoms () * 3);
   for (int is = 0, counter = 0; is < tau.getNSpecies (); is++)
      for (int ia = 0; ia < tau.getNAtoms (is); ia++)
         for (int i = 0; i < 3; i++)
            massVec(counter++) = speciesData.ionicMass(is);
}

SxVector3<double> SxArtifactFilter::translateToCOM (SxVecRef<double> *tauVec)
{
   int size = (int)(*tauVec).getSize ();
   SxVector3<double> COM(0., 0., 0.);
   const SxVecRef<double> &masses = getMasses ();
   double M = 0.;

   for (int i = 0; i < size;) {
      M += masses(i);
      COM(0) += (*tauVec)(i)*masses(i); i++;
      COM(1) += (*tauVec)(i)*masses(i); i++;
      COM(2) += (*tauVec)(i)*masses(i); i++;
   }
   COM /= M;
   /// subtract center of mass from tauVec
   for (int i = 0; i < size; i+= 3)
      SxVector3<double>::toVec3Ref (&(*tauVec)(i)) -= COM;
   return COM;
}

SxVector<double>
SxArtifactFilter::getInertiaTensor (const SxVecRef<double> &tauVec)
{
   int size = tau.getNAtoms () * 3;
   SxVector<double> T(3, 3);
   const SxVecRef<double> &masses = getMasses ();
   T.set (0.);
   for (int i3 = 0; i3 < size; i3+=3) {
      double M = masses(i3);
      T(0, 0) += M * (  tauVec(i3 + 1)*tauVec(i3 + 1)
                      + tauVec(i3 + 2)*tauVec(i3 + 2));
      T(1, 1) += M * (  tauVec(i3 + 0)*tauVec(i3 + 0)
                      + tauVec(i3 + 2)*tauVec(i3 + 2));
      T(2, 2) += M * (  tauVec(i3 + 0)*tauVec(i3 + 0)
                      + tauVec(i3 + 1)*tauVec(i3 + 1));

      T(0, 1) -= M * (tauVec(i3 + 0)*tauVec(i3 + 1));
      T(0, 2) -= M * (tauVec(i3 + 0)*tauVec(i3 + 2));
      T(1, 2) -= M * (tauVec(i3 + 1)*tauVec(i3 + 2));
   }
   T(1, 0) = T(0, 1);
   T(2, 0) = T(0, 2);
   T(2, 1) = T(1, 2);
   return T;
}

void SxArtifactFilter::validate (const SxVecRef<int> &/*equivalentIdx*/)
{
   // empty
}

