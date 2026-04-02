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

#ifndef _SX_QUAMOL_H_
#define _SX_QUAMOL_H_

#include <SxExt.h>
#include <SxVector.h>
#include <SxAtomicStructure.h>
#include <SxGkBasis.h>
#include <SxRadBasis.h>
#include <SxRadialBasis.h>
#include <SxAtomicOrbitals.h>
#include <SxQuantumNumbers.h>
#include <SxConstants.h>
#include <SxTimer.h>
#include <SxPW.h>

typedef SxArray<SxVector<SxComplex16> > SxOrbitals; //(ik,iOrbital:ig)

class SX_EXPORT_EXT SxQuamolRedSubSpace   {
   public:
      /// \brief Standart Constructor
      SxQuamolRedSubSpace ();
      /// \brief Standart Destructor
      ~SxQuamolRedSubSpace ();
      /// \brief set function for jsbBasis
      void set ( SxString fileName,
                 SxConstPtr<SxRadBasis> radBasisPtrIn,
                 SxPW &waveIn,
                 SxGkBasis *GkPtrIn,
                 SxConstPtr<SxRadialBasis> radGBasisPtrIn,
                 SxPtr<SxOverlapBase> SPtrIn,
                 SxAtomicStructure structureIn,
                 double rCut);
      /// \brief set Function for Basis reduction
      void set ( SxString fileName,
                 SxAtomicOrbitals initGuess,
                 SxConstPtr<SxRadBasis> radBasisPtrIn,
                 SxPW &waveIn,
                 SxGkBasis *GkPtrIn,
                 SxConstPtr<SxRadialBasis> radGBasisPtrIn,
                 SxPtr<SxOverlapBase> SPtrIn,
                 SxAtomicStructure structureIn);
      /// \brief QUasi AtoMic OrbitaLS (QUAMOLS) in radBasis
      SxAtomicOrbitals basis,improvedBasis;

      SxOrbitals basisGk;

      /// \brief Waves
      SxPW waves;
      /// \brief Plane-Wave Quamol Projection \f$\langle\Psi|b\rangle \f$
      SxArray<SxVector<SxComplex16> > beta,
      /** \brief Expansion coefficients matrix 
         \f$ \beta_{m,\mu} = \sum\limits_n{[S^b_{m,n}]^{-1}
                       \langle b_n|c_\mu \rangle} \f$ 
      */             
                          Sb;
      SxArray<SxArray<SxVector<SxComplex16> > > P;

      bool restartCG;

                          

      SxAtomicStructure structure;

      /// \brief RadBasisPtr
      SxConstPtr<SxRadBasis> radBasisPtr;
      /// \brief pointer to radial G basis
      SxConstPtr<SxRadialBasis> radGBasisPtr;
      /// \brief GkBasisPtr
      SxGkBasis *GkBasisPtr;

      /// \brief Overlap Operator
      SxPtr<SxOverlapBase> SPtr;

      void info (SxString fileName, SxConstPtr<SxRadBasis> radBasisPtrIn);

      void checkInitialGuess ();

      SxArray<SxArray<SxVector<SxComplex16> > > getProjections ();

      SxArray<SxVector<SxComplex16> > getCoefficients ();

      SxArray<SxVector<SxComplex16> > getOverlapC ();

      SxArray<SxVector<SxComplex16> > getOverlapB ();
      
      /** \brief Calculate the squarenorm: 
          \f$\sum\limits_{states,k}{w_k (P\beta[S^c]^{-1}\beta^{\dagger}P^{\dagger})_{states,k;states,k} 
      */   
      
      double calcNormB ();
      double calcNormC ();


      /** \brief calculate Gradient for SxOrbitals orbitals and transform back into SxRadials:
        \f$\sum\limits_{states,k}   {
         (P\beta[S^c]^{-1})_{states,k;\pi})P^{\dagger}_{p;staes,k}
         -(P\beta[S^c]^{-1})_{staes,k;\pi}
         (S^b\beta[S^c]^{-1}\beta P^{\dagger})_{p;states,k}  
         }\f$
      */   
      SxArray<SxVector<SxComplex16> > calcGradient ();

      void keepSymmetry (SxVector<SxComplex16> *MatIN);

      /// \brief Update orbitals with (Gk|radials)
      SxOrbitals radialsToOrbitals (const SxAtomicOrbitals &radialsIN);

      //void betaToRadial ();

      int getOrbital(SxArray<SxQuantumNumbers> &map, int is, int n, int l) const;

      SxVector<SxComplex16> reduce (SxArray<SxVector<SxComplex16> > &matIn) const;

      SxArray<SxVector<SxComplex16> > expand (SxVector<SxComplex16> &matIn) const;

      void redBetaToRadial ();

      /*
      /// \brief Update radials with (rad|orbitals)
      SxAtomicOrbitals orbitalsToRadials (
            const SxOrbitals &orbitalsIN,
            const SxAtomicOrbitals &radialsIN);
            */

      /// \brief \f$\langle\mu|\nu\rangle\f$
      double dotproductRad (
            const SxVecRef<double> &vec1,
            const SxVecRef<double> &vec2);

      double error;

      int maxSteps;

      void compute ();

      SxArray<SxVector<SxComplex16> > addBeta (
            const SxArray<SxVector<SxComplex16> > &beta1,
            const SxArray<SxVector<SxComplex16> > &beta2);

      SxArray<SxVector<SxComplex16> > skalarMultBeta (
            const double skalar,
            const SxArray<SxVector<SxComplex16> > &betaIN);
      
      SxComplex<double> lineMin (const SxArray<SxVector<SxComplex16> > &dir, double sw);

      SxComplex<double> parabelFit (double,double,double,double,double,double);

      SxAtomicOrbitals getJSB (double rCut, int lMax, int nZeros);

      SxArray<SxVector<double> > getJSBZeros (int lMax, int nZeros);

      double findRootJSB (int l, double guess);

      SxVector<double> jsb (int l, const SxVecRef<double> &z) const;

};

namespace Timer {
   enum QuamolTimer {
      normCalc,
      gradCalc,
      rad2Orb,
      setup,
      overlapB,
      projections
   };
}

SX_REGISTER_TIMERS (Timer::QuamolTimer)
{
   using namespace Timer;
   regTimer (setup,      "Basis setup");
   regTimer (rad2Orb,    "Radial to Orbitals");
   regTimer (overlapB,   "<b_i|b_j>");
   regTimer (projections,"<Psi|b_i>");
   regTimer (normCalc,   "Norm Calculation");
   regTimer (gradCalc,   "Gradient Calculation");
}

#endif /* _SX_QUAMOL_H_ */
