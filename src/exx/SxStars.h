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

#ifndef _SX_STARS_H_
#define _SX_STARS_H_

#include <SxExx.h>
#include <SxConfig.h>
#include <SxGkBasis.h>


/** \brief Calculates the star(s) of a k-point set.

    \b SxStars = S/PHI/nX Stars of k-points

    \author Matthias Wahnovic, wahn@fhi-berlin.mpg.de */
class SX_EXPORT_EXX SxStars
{
   public:

      SxStars ();
      SxStars (const SxCell                      &cell_,
               const SxArray<SxVector3<PrecG> > &kVec_,
               const SxVecRef<PrecWeights>      &weightsOrig_,
                     bool                         useInvSymmetry_=true);
      ~SxStars ();

      void compute ();

      void computeStar (int ik);

      /// Get a star (in relative coordinates)
      SxArray<SxVector3<PrecG> > getStar(int ik) const;
      PrecWeights getWeight (int ik) const;
      SxArray<SxArray<int> > getRotations () const;
      int getRotation (int ik, int iRot) const;

      int  getSizeOfMaxStar () const;

   protected:

      SxCell                                  cell;
      SxArray<SxVector3<PrecG> >             kVec;
      bool                                    useInvSymmetry;

   public:

      int                                     nkOrig;

      /// The stars (in relative coordinates)
      SxArray<SxArray<SxVector3<PrecG> > >   stars;
      SxArray<PrecWeights>                    weights;
      SxArray<SxArray<int> >                  rotations;

   protected:
      SxVector<PrecWeights>                weightsOrig;

};

#endif /* _SX_STARS_H_ */
