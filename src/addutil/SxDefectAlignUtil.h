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

#ifndef _SX_DEFECTALIGN_UTIL_H_
#define _SX_DEFECTALIGN_UTIL_H_

#include <SxAddUtil.h>
#include <SxAtomicStructure.h>
#include <SxTypes.h>
#include <SxFFT1d.h>
#include <SxCLI.h>

// --- Utility functions for sxdefectalign tools

namespace SxDefectAlignUtil {

enum FileType {
   None,
   sxb,
   VASP_LOCPOT,
   Socorro,
   QuantumEspresso
};


SX_EXPORT_ADD_UTIL SxMeshR getPot
   (SxCell &cell, SxMesh3D &mesh, SxString &fileName,
    FileType fileType, bool isHartree, SxAtomicStructure *structure);

class SX_EXPORT_ADD_UTIL SxRemesh1D {
   public:
      SxFFT1d oldFFT, newFFT;
      int mul;
      /// Constructor
      SxRemesh1D (ssize_t oldN, ssize_t newN, int mulIn = 1)
         : oldFFT(SxFFT::Reverse, (int)oldN),
           newFFT(SxFFT::Forward, (int)newN),
           mul(mulIn)
      {
         SX_CHECK (mul > 0);
      }
      /// Do the actual thing
      SxMeshR remesh (const SxMeshR &pot);
};

/** \brief Fourier interpolate periodic 1D function into a new mesh
    @param pot the 1D data
    @param N   the new mesh
    @param mul repeat pot multiple times during interpolation
    @return data on the new mesh

    Use this for few-time uses. For others, set up SxRemesh1D once
    to keep the FFT plan and avoid SxFFT plan cache explosion.
*/
inline SxMeshR remesh1D (const SxMeshR &pot, ssize_t N, int mul = 1)
{
   return SxRemesh1D ((int)pot.getSize (), N, mul)
          .remesh (pot);
}

SX_EXPORT_ADD_UTIL SxVector<double> readLine
                          (const SxCell &potCell,
                           const SxMesh3D &potMesh,
                           const SxMeshR &potData,
                           int idir,
                           ssize_t n,
                           const SxCell &cell,
                           const SxString &file);

SX_EXPORT_ADD_UTIL SxVector<double> average (const SxVecRef<double> &x, double w);

SX_EXPORT_ADD_UTIL int getMappedAtom (const Coord &pos,
                                 SxAtomicStructure &structure,
                                 const int iSpecies);

SX_EXPORT_ADD_UTIL FileType getFileType (SxCLI &cli);

}

#endif /* _SX_DEFECTALIGN_UTIL_H_ */
