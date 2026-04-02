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

#ifndef _SX_HESSIAN_H_
#define _SX_HESSIAN_H_

#include <SxStruct.h>
#include <SxBinIO.h>
#include <SxAtomicStructure.h>

/** \brief Hessian matrix

    \b SxClass = S/PHI/nX Hessian matrix

    \author C. Freysoldt */
class SX_EXPORT_STRUCT SxHessian
{
   protected:

      /// The Hessian matrix elements (iAtom:iNeighbor)
      SxArray<SxArray<SxMatrix3<double> > > hessian;

      /// The sparsity
      SxArray<SxVector<int> > neighbors;

      /// The cell mapping (optional)
      SxArray<SxArray<int> > offset;

   public:
      /// Empty constructor
      SxHessian () { }
      /// Constructor from full matrix
      SxHessian (const SxVecRef<double> &full);

      /// Construct full matrix
      SxVector<double> getFull () const;

      /// Write full matrix to binary file
      static void writeFull (SxBinIO &io,
                             const SxVecRef<double> &full);

      static void write (const SxString &fileName,
                         const SxVecRef<double> &hessianFull,
                         const SxAtomicStructure &str);

      /// Read full matrix from binary file
      static SxVector<double> readFull (const SxBinIO &io);

};

#endif /* _SX_HESSIAN_H_ */
