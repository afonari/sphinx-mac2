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

#include <SxOverlap.h>

/// Normalize states x
void SxOverlapBase::normalize (SxVecRef<SxComplex16> *xPtr) const
{
   SX_CHECK (xPtr);
   ssize_t n = xPtr->getNCols ();
   if (n <= 1)  {
      double nrm2 = normSqr (*xPtr);
      SX_CHECK (nrm2 > 0., nrm2);
      *xPtr *= 1. / sqrt(nrm2);
   } else {
      for (ssize_t i = 0; i < n; ++i)  {
         SxVecRef<SxComplex16> Xi = xPtr->colRef(i);
         double nrm2 = normSqr (Xi);
         SX_CHECK (nrm2 > 0., nrm2);
         Xi *= 1. / sqrt(nrm2);
      }
   }
}
                         
