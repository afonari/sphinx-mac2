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

#ifndef _SX_SPACESX_SEPARATOR_H_
#define _SX_SPACESX_SEPARATOR_H_

#include <SxBinIO.h>
#include <SxMesh3D.h>
#include <SxNeighbors.h>
#include <SxRBasis.h>
#include <SxRho.h>

class SxSpaceSeparator
{
   public:

      /// Constructor
      SxSpaceSeparator () { /*empty*/ };

      SxSpaceSeparator (const SxAtomicStructure &structureIn);

      ~SxSpaceSeparator () { /*empty*/ };

      SxArray<SxVector<double> > voronoi (const SxRBasis &R);

      SxArray<SxVector<double> > bader (const SxVecRef<double> &rhoIn);

      SxArray<SxVector<double> > baderTrinkle (
            const SxVecRef<double> &rhoIn);

      ssize_t getNextPerGrad (ssize_t index, SxVector3<double> &corr);

      SxArray<int> getNextAtom (ssize_t index);

      double getShareFactor(ssize_t index, 
            const SxVecRef<int> &atomicSpace, int iAtom);

      int findMaximum (ssize_t index, SxVecRef<int> &currentMap,
                       SxVector3<double> &corr, bool setTrajectory);
      int getCriticalType (ssize_t index);

      bool isBorderPoint (ssize_t index, SxVecRef<int> &currentMap);

      SxVector<double> getFilter (const SxVecRef<int> &atomicSpace, int iAtom);

      SxArray<SxVector<double> > computeGradRho ();
      SxArray<SxVector<double> > computeGradGradRho ();
      SxVector<double>& getGradRho(int dir);
      SxVector<double>& getGradGradRho(int dir);

      ssize_t getNMaxima() {return maximaIndices.getSize();};
      SxArray<ssize_t> getMaximaIndices() {
         if (getNMaxima () > 0) return SxArray<ssize_t>(maximaIndices);
         else return SxArray<ssize_t> (0);
      }
      ssize_t getNMinima() {return minimaIndices.getSize ();};
      SxArray<ssize_t> getMinimaIndices() {
         if (getNMinima () > 0) return SxArray<ssize_t>(minimaIndices);
         else return SxArray<ssize_t> (0);

      };
      ssize_t getNRing() {return ringIndices.getSize ();};
      SxArray<ssize_t> getRingIndices() {
         if (getNRing () > 0) return SxArray<ssize_t>(ringIndices);
         else return SxArray<ssize_t> (0);
      };
      ssize_t getNBond() {return bondIndices.getSize ();};
      SxArray<ssize_t> getBondIndices() {
         if (getNBond () > 0) return SxArray<ssize_t>(bondIndices);
         else return SxArray<ssize_t> (0);
      };

      void plotDensity(const SxString &filename);
      void plotGradDensity(const SxString &filename);

      SxVector<double> PAWDist;

   protected:

      enum pointType {Maximum, Interior, Boundary};

      SxVector<double> rho;

      SxArray<SxVector<double> > gradRhoComponents;

      SxArray<SxVector<double> > gradGradRhoComponents;

      SxAtomicStructure structure;

      SxMesh3D mesh;

      SxStack<ssize_t> maximaIndices;
      SxStack<ssize_t> ringIndices;
      SxStack<ssize_t> bondIndices;
      SxStack<ssize_t> minimaIndices;
      
      SxArray<ssize_t> findNeighbors(ssize_t index, int shell);

      double ramp (double u);

      SxVector<double> getProbabilityFlux (ssize_t index,
            const SxVecRef<double> &rhoIn);

      void printWeights (const SxArray<SxVector<double> > &weights, ssize_t idx);

};

namespace Timer {
   enum SpaceSeparatorTimer {
      Voronoi,
      Bader,
      BaderTrinkle
   };
}

SX_REGISTER_TIMERS(Timer::SpaceSeparatorTimer)
{
   using namespace Timer;
   regTimer (Voronoi,     "Voronoi");
   regTimer (Bader,       "Bader bond critical");
   regTimer (BaderTrinkle,"Bader-Trinkle");
}

#endif /* _SX_SPACESX_SEPARATOR_H_ */
