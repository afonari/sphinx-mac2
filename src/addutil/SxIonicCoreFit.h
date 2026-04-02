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

#ifndef _SX_IONIC_CORE_FIT_H_
#define _SX_IONIC_CORE_FIT_H_

#include <SxAddUtil.h>
#include <SxAtomicStructure.h>
#include <SxGrid.h>

/** \brief Fit ionic cores to ions on a mesh

    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_ADD_UTIL SxIonicCoreFit
{
   public:
      /// delta r
      double dr;
      /// maximum radius
      double rMax;
      /// order of splines
      int splineOrder;
      /// number of radial points
      int nPoints;

   protected:
      /// 3D mesh on which data is given
      SxMesh3D mesh;
      /// real size of element of the mesh
      SxCell meshCell;

      /// The atomic structure
      SxAtomicStructure structure;
      /// The corresponding grid for neighbor searches
      SxGrid grid;

   public:
      /// List of all elements
      SxArray<SxString> chemNames;

      /// Mapping of species for current structure to chemName
      SxArray<int> chemId;

      /// Set structure & mesh
      void set (const SxAtomicStructure &structureIn, const SxMesh3D &meshIn);

      /// Set shape parameters
      void set (int nPts, double rMaxIn, int order);

      /// Get parameter derivative wrt data
      SxVector<double> getParamDeriv ();

      /// Evaluate shape function given by Cn at point x
      double getVal (double x, const SxVecRef<double> &Cn,
                     SxVector<double> &work) const;
      inline double getVal (double x, const SxVecRef<double> &Cn) const
      {
         SxVector<double> work(splineOrder + 1);
         return getVal (x, Cn, work);
      }

      /// Evaluate overlap of spherical functions on mesh
      SxVector<double> getFunc () const;

      /// Evaluate overlap of spherical functions on mesh
      SxVector<double> getFunc (const SxAtomicStructure &structureIn,
                                const SxMesh3D &meshIn)
      {
         set (structureIn, meshIn);
         return getFunc ();
      }

      void addFitData (const SxVector<double> &data,
                       const SxMesh3D &meshIn,
                       const SxAtomicStructure &structureIn,
                       double weight = 1.);
      void computeFit (bool zeroAvg = false);
   protected:
      /// Fitting variables
      SxVector<double> fitA, fitb, param;
   public:
      /// Get parameters after computeFit
      const SxVector<double> &getParam () const
      {
         return param;
      }

      /// Get parameters after computeFit
      const SxVecRef<double> getParam (ssize_t iSpecies) const
      {
         return const_cast<SxVector<double>&>(param).colRef (iSpecies);
      }


};

#endif /* _SX_IONIC_CORE_FIT_H_ */
