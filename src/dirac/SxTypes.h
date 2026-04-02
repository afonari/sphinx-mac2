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

#ifndef _SX_TYPES_H_
#define _SX_TYPES_H_

#include <SxMatrix3.h>
#include <SxVector.h>
#include <SxArray.h>
#include <SxPrecision.h>

// --- wave functions
typedef SxVector<PrecCoeffG>                  PsiG;
typedef SxVector<PrecCoeffR>                  PsiR;
typedef SxVecRef<PrecCoeffR>                  PsiRef;
// --- potentials, charge density
typedef SxVector<PrecRhoG>                    SxMeshG;
typedef SxVector<PrecRhoR>                    SxMeshR;
typedef SxArray<SxMeshR>                         VxcR;    // :iSpin
typedef SxArray<SxMeshG>                         RhoG;    // :iSpin
typedef SxArray<SxMeshR>                         RhoR;    // :iSpin


#endif /* _SX_TYPES_H_ */
