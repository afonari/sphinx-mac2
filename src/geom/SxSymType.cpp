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

#include <SxSymType.h>

SxSymType::SxSymType (Classification c_, 
                      int axisCount_,
                      const SxVector3<double> &opCoord_)
   : classification(c_), 
     axisCount(axisCount_), 
     opCoord(opCoord_),
     identifier("")
{
   updateIdentifier ();
}


SxSymType::SxSymType (const SxSymType &in)
   : classification (in.classification),
     axisCount      (in.axisCount),
     opCoord        (in.opCoord),
     identifier     ()
{
   updateIdentifier ();
}


void SxSymType::updateIdentifier ()
{
   const SxVector3<double> &v = opCoord;
   switch (classification)  {
      case Identity        :  identifier = "E";
                              break;
      case Inversion       :  identifier = "-1";
                              break;
      case Mirror          :  identifier = SxString() + "m [" 
                              + v(0) + "," + v(1) + "," + v(2) + "]";
                              break;
      case Rotation        :  identifier = SxString() + "C" 
                              + axisCount + " ["
                              + v(0) + "," + v(1) + "," + v(2) + "]";
                              break;
      case MirrorRotation  :  identifier = SxString() + "S" 
                              + axisCount + " ["
                              + v(0) + "," + v(1) + "," + v(2) + "]";
                              break;
      default              :  identifier = "unknown";
   }
}


bool SxSymType::operator== (Classification in) const
{
   return (classification == in);
}


bool SxSymType::operator== (const SxSymType &in) const
{
   return (   classification == in.classification
           && axisCount      == in.axisCount
           && (opCoord - in.opCoord).absSqr().sum() < 1e-10 );
}

int SxSymType::getRank () const
{
   switch (classification)  {
      case Identity:
         return 1;
      case Inversion:
         return 2;
      case Mirror:
         return 3;
      case Rotation:
         return 4 + 2*axisCount;
      case MirrorRotation:
         return 5 + 2*axisCount;
      default:
         return 0;
   }
}


