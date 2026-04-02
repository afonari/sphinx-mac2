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

#ifndef _SX_RAD_BASIS_H_
#define _SX_RAD_BASIS_H_

#include <SxPrecision.h>
#include <SxQuantumNumbers.h>
#include <SxBasis.h>
#include <SxGBasis.h>
#include <SxVector.h>
#include <SxArray.h>
#include <SxTimer.h>
#include <SxDiracLib.h>

/** This class describes a radial mesh basis. It will be used e.g. for
    atomic wavefunctions and in PAW part. In contrast to most of the
    other basis modules it contains a intermediate function to convert
    a wavefunction given on a radial grid to a plane wave grid.

    \b SxRadBasis = S/PHI/nX Radial Basis

    \ingroup group_dft
    \author  Sixten Boeck, Christoph Freysoldt
    */
class SX_EXPORT_DIRAC SxRadBasis : public SxBasis
{
   public:

      typedef double              CoeffType;
      typedef SxVector<CoeffType> TPsi;

      SxRadBasis ();
      SxRadBasis (const SxArray<SxVector<double> > &radFunc,
                  const SxArray<Real8>               &logDr);
      SxRadBasis (double rMin, double rMax, int nPoints);
      SxRadBasis (const SxVector<double> &rMin, const SxVector<double> &rMax,
            const SxVector<int> &nPoints);
      SxRadBasis (const SxString &file);
      SxRadBasis (const SxBinIO &io);
      void set   (const SxArray<SxVector<double> > &radFunc,
                  const SxArray<Real8>               &logDr);
      void set   (double rMin, double rMax, int nPoints);
      void set (const SxVector<double> &rMin, const SxVector<double> &rMax,
            const SxVector<int> &nPoints);

      virtual ~SxRadBasis ();
      void read (const SxString &file);
      void read (const SxBinIO &io);
      SxRadBasis readMesh (const SxString &file);
      SxRadBasis readMesh (const SxBinIO &io);

      /** \brief Add radial mesh
        \param radFuncIn radial functions for new species for different l
        \param logDrIn   \f$ \ln \frac{r(i+1)}{r(i)}\f$ for different l
        \return new species id

        The different radial meshes in Dirac projections are given by the
        auxData.is information of the Dirac vector. The return
        value of this function gives you the 'is' to access the radial mesh
        that is added by this function.

        Note that adding meshes in this way is less efficient than
        giving all in the constructor, and that all internal cashes
        are cleaned. Avoid calling this routine very often.

        \todo Single-shot routine for toPWBasis ()
        */
      int addMesh (const SxVecRef<double> &radFuncIn,
                   double             logDrIn);

      /// Dummy function, since number of elements may vary
      virtual ssize_t getNElements () const {
         SX_EXIT;
         return -1;
      }

   protected:
      void registerMemoryObservers ();

   public:
      int getNSpecies () { return int(radFunc.getSize()); }
      /// Return basis type
      virtual SxString getType () const { return "|r>"; }

      REGISTER_PROJECTOR (SxRadBasis,  changeRadBasis);
      REGISTER_PROJECTOR (SxGBasis,    toPWBasis);
      REGISTER_PROJECTOR (SxRadialBasis, toRadialBasis);

      SxVector<CoeffType> changeRadBasis (const SxRadBasis *basis,
                                          const SxVecRef<CoeffType> &vec) const;
      /**
        \brief Project to plane-wave basis

        \Note is the iAtom field in the input vector set, i.e., other
                 than -1, the translation operator is applied,
                 otherwise not...
       */
      SxVector<SxComplex16> toPWBasis (const SxGBasis *,
                                       const SxVecRef<CoeffType> &) const;
      SxVector<double> toRadialBasis (const SxRadialBasis *,
                                      const SxVecRef<CoeffType> &) const;

      /// Integrate over radial space with r^2 dr
      virtual Real8 tr (const SxVecRef<double> &) const;

      Real8 cosTheta (int ig, int jg, const SxGBasis &g);
      Real8 pl (int l, Real8 x);

      double integrate (const SxVecRef<double> &) const;

      /// Radial mesh, iSpecies:l:ir
      SxArray<SxVector<double> >  radFunc;  // :is,:r
      SxArray<Real8>                logDr;    // :is
      double getRMax (int is) const { return radFunc(is)(radFunc(is).getSize() - 1);};
      double getRMax () const;
      double getRMin (int is) const { return radFunc(is)(0);};
      ssize_t getNPoints (int is) const { return radFunc(is).getSize();};

      int getNSpecies () const {
         SX_CHECK (logDr.getSize () == radFunc.getSize (),
                   logDr.getSize (), radFunc.getSize ());
         return int(radFunc.getSize ());
      }

      /// \name Cashing
   protected:
      //@{
      /** The most expensive part in a <G+k|RY> projection are the integrals
        \f[ \int_0^\infty r^2 dr j_{l}(|G+k|\cdot r) \cdot R(r) \f]
        because they have to be calculated for each G vector (well, at least
        for each possible length |G+k|). Because this doesn't depend on m, we
        should cache these integrals for the current l and vector and all
        known G+k bases.
        \brief Cashed spherical Bessel integrals ik:ig
        */
      mutable SxArray<SxVector<CoeffType> > cashedRl;
      /// The vector the integrals of which are cashed
      mutable SxVector<double> cashedVec;
      /// The known G bases, :ik:
      // no SxArray, because it's just for finding ik, where lists are
      // perfect
      mutable SxList<const SxGBasis *> gBases;

      //@}

   public:

      enum { S=0, IPY, IPZ, IPX, DXY, DYZ, DZ2, DXZ, DX2_Y2 };

      /** \brief Get real, unnormalized Ylm for all l<=lmax and all m
           for a shifted G basis
           \todo Move to SxGBasis (as member)
        */
      static SxVector<double> realYlm (int lmax, const SxGBasis &G,
                                         const Coord &dG);

      /**
        \brief Calculates the projection of the radial part onto plane waves

        This routine should not be used from the outside. Set up a proper
        radial basis containing all the radial meshes you need and use the
        Dirac projector notation.

        \Note The normalization factor of \f$\frac{4\pi}{\sqrt\Omega}\f$
              is NOT included here, because it is combined with the Ylm
              normalization factors later on.
        \todo Shouldn't be called from outside, so make it protected
        */
      SxVector<double>
      toPWBasis (const SxVecRef<double> &rad,
                 const SxVecRef<double> &psi,
                 const SxGBasis &pwBasis,
                 int l, Real8 logDr) const;

      /**
        \brief cutoffvalue for (Rad|G)(G|Psi) Transformation
        */
      mutable double cutoff;

   public:
      /** \brief Specific radial basis for projections

          SxRadBasis is a container class for the radial meshes of all
          species. Sometimes you'd like to pick a particular mesh.
          This is what special basis is doing for you - it combines the
          multi-species container with meta-data (species, l, m).

          Don't treat this as a real basis.
        */
      class SpecialBasis
      {
         friend class SxRadBasis;
         /// The full radial basis
         const SxRadBasis &radBasis;
         public:
            /// Specific data for radial basis
            int is, l, m;
         protected:
            /// Constructor
            inline
            SpecialBasis (const SxRadBasis &fullBasis, int isIn,
                          int lIn = -1, int mIn = 0)
            : radBasis(fullBasis), is(isIn), l(lIn), m(mIn)
            {
               SX_CHECK (is >= 0 && is < radBasis.radFunc.getSize (),
                         is, radBasis.radFunc.getSize ());
               SX_CHECK (l < 0 || abs(m) <= l, l, m);
            }
         public:
            /// Get number of elements
            int getSize () const  {
               return (int)radBasis.radFunc(is).getSize ();
            }

            /// SxVector constructor
            inline
            operator SxVector<double> () const
            {
               SxVector<double> res(getSize ());
               res.setBasis (&radBasis);
               res.auxData.is = is;
               res.auxData.l  = char(l);
               res.auxData.m  = char(m);
               return res;
            }

      };
      /// Return species-specific radial basis
      inline SpecialBasis operator() (int is) const {
         SX_CHECK(is < radFunc.getSize (), is, radFunc.getSize ());
         return SpecialBasis (*this, is);
      }
};


namespace Timer {
   enum RadBasisTimer {
      JsbInt,
      JsbCalc,
      Phase,
      RadTotal,
      rad2radG,
      rad2radR
   };
}

SX_REGISTER_TIMERS (Timer::RadBasisTimer)
{
   using namespace Timer;
   regTimer (JsbInt,   "jsb Integrals");
   regTimer (JsbCalc,  "jsb Setup");
   regTimer (RadTotal, "radBasis");
   regTimer (rad2radG, "rad to radG Basis");
   regTimer (rad2radR, "rad to radR Basis");
}

/** \brief Stand alone radial interpolation function
     @param psi function to be interpolated
     @param r0     start of original logarithmic grid
     @param logDr  logDr of original logarithmic grid
     @param newRad new grid
     @return interpolated function

     Interpolate psi from a logarithmic grid (defined by r0 and logDr)
     to new grid given in newRad.
     \author C. Freysoldt
  */
SxVector<double> SX_EXPORT_DIRAC
interpolateRad(const SxVecRef<double> &psi,
               double r0, double logDr,
               const SxVecRef<double> &newRad);
#endif /* _SX_RAD_BASIS_H_ */
