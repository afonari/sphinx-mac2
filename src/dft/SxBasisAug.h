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

#ifndef _SX_BASIS_AUG_H_
#define _SX_BASIS_AUG_H_

#include <SxDFT.h>
#include <SxBasis.h>

class SxPAWBasis;

/** \brief Basis class for augmented basis sets


    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxBasisAug : public SxExtraBasis<SxComplex16>
{
   public:
      /** \brief placeholder for any <R|X> projector */
      VIRTUAL_PROJECT_FROM (SxPAWBasis, SxComplex16);

};

#endif /* _SX_BASIS_AUG_H_ */
