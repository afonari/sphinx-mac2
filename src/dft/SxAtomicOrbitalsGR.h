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

#ifndef _SX_ATOMIC_ORBITALSGR_H_
#define _SX_ATOMIC_ORBITALSGR_H_

#include <SxDFT.h>
#include <SxVector.h>
#include <SxArray.h>
#include <SxBinIO.h>
#include <SxAtomicStructure.h>
#include <SxRadialBasis.h>

class SX_EXPORT_DFT SxAtomicOrbitalsGR
{
   public:

      SxAtomicOrbitalsGR ();
      SxAtomicOrbitalsGR (const SxArray<SxArray<SxVector<double> > >  &in,
                         const SxConstPtr<SxRadialBasis> &radialBasisPtrIn,
                         bool splineRepIn);
      
      /// \brief Copy Constructor
      SxAtomicOrbitalsGR (const SxAtomicOrbitalsGR &);

      /// \brief Destructor
      virtual ~SxAtomicOrbitalsGR ();

      void set (double val);

      /// \brief Assignment operator
      void operator= (const SxAtomicOrbitalsGR &);

      void operator+= (const SxAtomicOrbitalsGR &);
      void operator-= (const SxAtomicOrbitalsGR &);

      SxAtomicOrbitalsGR operator+ (const SxAtomicOrbitalsGR &) const;
      SxAtomicOrbitalsGR operator- (const SxAtomicOrbitalsGR &) const;

      SxAtomicOrbitalsGR operator* (double skalar) const;
      SxAtomicOrbitalsGR operator/ (double skalar) const;

      void sumMPI() const;

      /// \brief Get number of Species
      int getNSpecies () const;

      int getNOrbTypes () const;
      /** \brief Gets the maximum angular component for a certain species.

          This function returns the index of the maximum angular component
          of a specified species.

          \b Note: The returned value is not the number of components. Thus,
          if used in for() loops the '<=' operator has to be used, e.g.
          \code
             // --- Loop over all angular components
             for (int l=0; l <= orbitals.getLMax(); l++)  {
                ...
             }
          \endcode
       */
      int getNOrbTypes (int iSpecies) const;

      SxArray<SxQuantumNumbers> getReducedOrbitalMap () const;
      
      SxArray<SxQuantumNumbers> getOrbitalMap (const SxAtomicStructure &structure) const;

      //// \brief get muSet(is)(iot)
      SxVector<double> &operator() (int is, int iot);
      const SxVector<double> &operator() (int is, int iot) const;

      SxVector<double> operator() (int is, int ia, int iot, int l, int m) const;
      
      void print (const SxString &fileIn) const;

      double getNormSqr (int is, int iot) const;
      double getNormSqrSum () const;

      double dot (const SxAtomicOrbitalsGR &orbitalsIn, int is, int iot, int js, int jot) const;
      double dot (const SxAtomicOrbitalsGR &orbitalsIn) const;

      double sum (const SxAtomicOrbitalsGR &orbitalsIn, int is, int iot, int js, int jot) const;
      double sum (const SxAtomicOrbitalsGR &orbitalsIn) const;

      void normalize ();

      void setBasis (const SxConstPtr<SxRadialBasis> &radialBasisPtrIn);

      SxConstPtr<SxRadialBasis> getRadGBasisPtr () const { return radialBasisPtr; };

      void toSpline () const;

      void toVec () const;

      SxVector<double> toVec (int is, int iot) const;

      bool isSpline () const {return splineRep;};

      void createFuncLMap ();

      SxRadialBasis::TPsi &getFuncL (int is, int l, int ifl);

      const SxRadialBasis::TPsi &getFuncL (int is, int l, int ifl) const;

      int getFuncPerL (int is, int l) const;

      SxArray<SxArray<SxVector<double> > > orthogonalize ();

      SxArray<SxArray<SxVector<double> > > getOverlap () const;

      SxVector<double> getOverlap (int is, int l) const;

      SxArray<SxArray<SxArray<int> > > funcLMap;

      int getLMax () const;

      int getLMax (const int iSpecies) const;

      void orthogonalizeOn(SxAtomicOrbitalsGR &basis);

      void rotate(const SxArray<SxArray<SxVector<double> > > &rotMat);

      int getNOrbitals (const SxAtomicStructure &structure) const;

      int getIOT (int is, int n, int l) const;

      int getOrbitalIdx (int is, int ia, int iot, int l, int m, SxArray<SxQuantumNumbers> map) const;

   protected:

      mutable bool splineRep;

      SxConstPtr<SxRadialBasis>  radialBasisPtr;

      /** \brief The sampling points of each orbitals.

        <b>Strorage order</b>: 
        -# iSpecies
        -# l
        -# r */
      mutable SxArray<SxArray<SxVector<double> > >  muSet;  // :is,:iot,:r

};

inline SX_EXPORT_DFT SxAtomicOrbitalsGR operator* (double skalar, const SxAtomicOrbitalsGR &in)
{
   return in * skalar;
}

#endif /* _SX_ATOMIC_ORBITALSGR_H_ */
