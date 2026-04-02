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

#ifndef _SX_R_BASIS_H_
#define _SX_R_BASIS_H_

#include <SxString.h>
#include <SxBasis.h>
#include <SxCell.h>
#include <SxSymGroup.h>
#include <SxFFT3d.h>
#include <SxBinIO.h>
#include <SxTypes.h>
#include <SxDiracLib.h>


/** 
  \author Sixten Boeck
  */
class SxGBasis;

/** \brief Realspace basis

  \b SxRBasis = S/PHI/nX Realspace Basis

  \ingroup group_dft
  \author  Sixten Boeck
  */
class SX_EXPORT_DIRAC SxRBasis : public SxBasis
{
   public:

      typedef PrecCoeffR             CoeffType;
      typedef SxVector<PrecCoeffR>   TPsi;
      typedef SxVector<PrecRhoR>     TRho;

      /// Return basis type
      virtual SxString getType () const { return "|R>"; }

      SxRBasis ();
      SxRBasis (int nx, int ny, int nz, const SxCell &);
      SxRBasis (const SxVector3<int> &, const SxCell &);

      /// Initialize basis
      void set (const SxVector3<int> &, const SxCell &);
      /// Initialize basis
      void set (int nx, int ny, int nz, const SxCell &);

      virtual ~SxRBasis ();

      template<class T>
      SxVecRef<T> symmetrize (const SxVecRef<T> &) const;
      void writeMesh3d (const SxString &file, const TRho &mesh3d) const;
      RhoR readMesh3d  (const SxString &file) const;

      int getMeshSize () const { return fft3d.meshSize; }
      SxVector3<int> getMesh () const { return fft3d.mesh; }

      virtual void print () const;

      /** \brief \f$ \Delta \Omega\f$ used for integration

          Sometimes realspace entities (like the charge density SxRho) have
          to be integrated. Since the realspace is represented on a regular
          grid (FFT grid) such integrations can numerically easily evaluated.
          If X(R) is an entity living on the realspace grid its integral 
          is
          \f[
             \int X(R) \Omega \Rightarrow \sum_i X(R_i) \Delta \Omega
          \f] and
          \f[ 
             \Delta \Omega = \frac{\Omega}{n_\mathrm{FFT}}
          \f]
          The number of realspace (FFT) mesh points are \f$ n_\mathrm{FFT}.\f$
          
       */
      Real8             dOmega;

      /// Return number of mesh points
      virtual ssize_t getNElements () const { return fft3d.meshSize; }
//   protected:

      /// \brief The FFT object used to transform objects between |G> and |R> 
      mutable SxFFT3d   fft3d;

      /// The unit cell
      SxCell            cell;

   protected:
      /// The dual G-basis
      mutable const SxGBasis* gBasisPtr;

   public:

      /** Register the dual G basis */
      inline void registerGBasis (const SxGBasis &gBasis) const;

      /** Get the dual G basis */
      const SxGBasis& getGBasis () const
      {
         SX_CHECK(gBasisPtr);
         return *gBasisPtr;
      }

      /// Check if dual G basis is available
      bool hasGBasis () const
      {
         return gBasisPtr;
      }

      /** Deregister a basis */
      virtual void deregister (const SxBasis *basis) const;

      REGISTER_PROJECTOR (SxRBasis, identity);
      REGISTER_PROJECTOR (SxGBasis, toGSpace);

      virtual Real8 tr (const SxVecRef<double> &) const;
      virtual Real8 tr (const SxVecRef<SxComplex16> &) const;

      SxVecRef<CoeffType> identity (const SxRBasis *,
                                    const SxVecRef<CoeffType> &) const;

      SxVector<PrecCoeffG> toGSpace (const SxGBasis *,
                                     const SxVecRef<CoeffType> &) const;
};

namespace Timer {
   enum RGBasisTimer {
      GBasis,
      Sym,
      GtoR, GtoR_FFT,
      RtoG, RtoG_FFT
   };
}

SX_REGISTER_TIMERS (Timer::RGBasisTimer)
{
   using namespace Timer;
   regTimer (GBasis,     "|G> Setup");
   regTimer (Sym,        "Symmetrization");
   regTimer (GtoR,       "G->R routine");
   regTimer (GtoR_FFT,   "G->R FFT call");
   regTimer (RtoG,       "R->G routine");
   regTimer (RtoG_FFT,   "R->G FFT call");
}

//------------------------------------------------------------------------------
// Symetrization of a 3d realspace mesh
//------------------------------------------------------------------------------
template<class T>
SxVecRef<T> SxRBasis::symmetrize (const SxVecRef<T> &meshIn) const
{
   SX_CHECK (meshIn.getBasisPtr() == this);
   SX_CLOCK (Timer::Sym);

   // --- do we have to symmetrize mesh???
   if (!cell.symGroupPtr || cell.symGroupPtr->getNSymmorphic () <= 1)  {
      return const_cast<SxVecRef<T>&>(meshIn);
   }

   const SxSymGroup &S = *cell.symGroupPtr;
   int nOp = S.getNSymmorphic ();

   SxCell meshCell (cell(0) / fft3d.mesh(0),
                    cell(1) / fft3d.mesh(1),
                    cell(2) / fft3d.mesh(2));

   SxArray<SxMatrix3<int> >  symOpRel (nOp);
   for (int iOp=0; iOp < nOp; iOp++)
      symOpRel(iOp) = meshCell.carToRel (S.getSymmorphic(iOp));

   SxVector3<int> mesh = fft3d.mesh;
   SxVector<T> meshOut (meshIn.getSize());
   meshOut.setBasis (meshIn.getBasisPtr());


   double nOpInv = 1. / double(nOp);
#ifdef USE_OPENMP
#pragma omp parallel
#endif
   {
      SxArray<ssize_t> idxRot(nOp);
      const int nj = mesh(0), nk = mesh(1), nl = mesh(2);
#ifdef USE_OPENMP
#pragma omp for collapse(2)
#endif
      for (int j=0; j < nj; j++)  {
         for (int k=0; k < nk; k++)  {
            for (int l=0; l < nl; l++)  {
               SxVector3<int> i(j,k,l);
               ssize_t iOut = fft3d.mesh.getMeshIdx(i, SxMesh3D::Positive);

               // --- collect all indices
               bool done = false;
               for (int iOp = 0; iOp < nOp; iOp++)  {
                  SxVector3<int> rot = symOpRel(iOp) ^ i;
                  ssize_t idx = fft3d.mesh.getMeshIdx(rot, SxMesh3D::Unknown);
                  // only lowest index does the computation
                  if (idx < iOut) { done = true; break; }
                  idxRot(iOp) = idx;
               }
               // only lowest index does the computation
               if (!done)  {
                  // --- compute sum of all equivalent elements ...
                  T res = 0.;
                  for (int iOp = 0; iOp < nOp; iOp++)
                     res += meshIn(idxRot(iOp));
                  res *= nOpInv;
                  // --- distribute result
                  for (int iOp = 0; iOp < nOp; iOp++)
                     meshOut(idxRot(iOp)) = res;
               }
            }
         }
      }
   }

   return meshOut;
}

#include <SxGBasis.h>

void SxRBasis::registerGBasis (const SxGBasis &gBasis) const
{
   registerBasis (gBasis);
   gBasisPtr = &gBasis;
}
#endif /* _SX_R_BASIS_H_ */
