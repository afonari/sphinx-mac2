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

#ifndef _SX_VDW_H_
#define _SX_VDW_H_

#include <SxDFT.h>
#include <SxAtomicStructure.h>
#include <SxSpeciesData.h>
#include <SxSymbolTable.h>

/** This class is for the empirical van der waals correction
    \author Lars Ismer
    \author Michael Ashton, ashton@mpie.de
 */
class SX_EXPORT_DFT SxVDW
{
   public:

      enum { Vdw_D2, Vdw_D3, Vdw_TS, Vdw_D3_TS, Vdw_TSI, Vdw_D3_TSI } method;
      enum { Vdw_GouldBucko, Vdw_Tang, Vdw_Default } combinationRule;
      enum { Vdw_Fermi, Vdw_BJ } damping;
      SxAtomicStructure coord;
      SxArray<int> speciesNum;
      SxArray<SxArray<int> > neighbors;
      SxAtomicStructure forces;
      SxArray<SxArray<double> > dist;
      SxArray<SxArray<SxVector3<double> > > unitVectors;
      SxArray<double> polarizability;
      SxArray<double> C6D2;
      SxArray<double> C6TS;
      SxArray<double> vdwRadiusD2;
      SxArray<double> vdwRadiusTS;
      SxArray<double> covalentRadius;
      SxArray<double> effectiveVolume;
      SxArray<double> d3Coordination;

      int nPairs;  // number of atom-atom pairs contributing to totalEnergy
      double totalEnergy;

      SxVDW ();
      SxVDW (const SxAtomicStructure &, const SxSymbolTable *);
     ~SxVDW ();

     /// Update the atomic structure (no compute)
     void update (const SxAtomicStructure &newcoord);
     /// Update the effective volumes (no compute)
     void update (const SxArray<double> &effectiveV);

   protected:
     void updateNeighbors ();
     void updateD3Coordination ();
   public:

     /// Whether effective volumes are needed
     bool needEffVolumes () const
     {
        return (method == SxVDW::Vdw_TS || method == SxVDW::Vdw_D3_TS);
     };

     void compute ();

     const SxAtomicStructure& getForces () const
     {
        return forces;
     }

   protected:
     /** \brief Compute damping function
       @param R     actual distance
       @param Rij   sum of van-der-Waals radii
       @param deriv here, the first derivative with respect to R is stored 
       @return the damping function

       @note: currently, only a Fermi-like function is implemented
       */
     double getDampingFunction (double R, double Rij, double *deriv);

     double getC6ij (int, int);
     double getRij (int, int);
     /// Get cutoff for van-der-Waals interactions
     double getVdwCutoff () const
     {
        return 94.5;  // 50 Angstroms -> Bohr
     }

   public:

     double getTotalEnergy () const
     {
        return totalEnergy;
     }

};


#endif /* _SX_VDW_H_ */
