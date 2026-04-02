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

#ifndef _SX_WANNIER_HAM_H_
#define _SX_WANNIER_HAM_H_


#include <SxConfig.h>
#include <SxDirTensor3.h>
#include <SxPW.h>
#include <SxGkBasis.h>
#include <SxFermi.h>
#include <SxCell.h>
#include <SxExt.h>


/** \brief Computes the Hamilton operator in the Wannier function basis set.

    \b SxWannierHam = S/PHI/nX Hamiltonian within the Wannier basis

    \author Matthias Wahn, wahn@fhi-berlin.mpg.de */
class SX_EXPORT_EXT SxWannierHam 
{
   public:
      /** empty constructor */
      SxWannierHam ();

      /** destructor */
      ~SxWannierHam ();

      /** dump the matrix \f$ \langle 0\,n | H | \mathbf{R}\,m \rangle \f$
          belonging to the lattice vector \f$ \mathbf{R} \f$ */
      static void printMatrix (const SxVector<SxComplex16> &M);
      /** write eigenenergies obtained from Wannier function representation
          to output file */
      static void printEps (const SxDirTensor3<double> &bs,
                            const SxString &file);
      /** write eigenenergies obtained from Wannier function representation
          on screen */
      static void printEps (const SxDirTensor3<double> &bs);
};

#endif /* _SX_WANNIER_HAM_H_ */
