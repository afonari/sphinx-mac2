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

#include <SxPAWBasis.h>
#include <SxRBasis.h>
#include <SxAOBasis.h>

void SxPAWBasis::print () const
{
   gBasis->print ();
   pBasis->print ();
}

SxVector<PrecCoeffG>
SxPAWBasis::toPWBasis (const SxGBasis *gIn,
                       const SxVecRef<SxComplex16> &psi) const
{
   SX_CHECK (gBasis);
   SX_CHECK (gIn == gBasis.getPtr ());
   int nPsi = (int)psi.getNCols ();
   if (nPsi == 0) nPsi = 1;
   SX_CHECK (psi.getSize () > 0);
   PsiG res (gBasis->getNElements (), nPsi);
   res.auxData = psi.auxData;
   res.setBasis (gIn);
   // --- copy G-vector part into new vector
   res <<= psi.getRef<SubMatrix> (0, gBasis->ng, psi.getNCols ());
   return res;
}

SxVecRef<SxComplex16, SubMatrix>
SxPAWBasis::toBasis (const SxGBasis *gIn,
                     const SxVecRef<SxComplex16> &psi) const
{
   SX_CHECK (gBasis);
   SX_CHECK (gIn == gBasis.getPtr ());
   SX_CHECK (psi.getNRows () == getNElements (),
             psi.getNRows (), getNElements ());
   SxVecRef<SxComplex16, SubMatrix> res
      = psi.getRef<SubMatrix> (0, gBasis->ng, psi.getNCols ());
   res.auxData = psi.auxData;
   res.setBasis (gIn);
   return res;
}

SxVecRef<SxComplex16, SubMatrix>
SxPAWBasis::toBasis (const SxPartialWaveBasis *pIn,
                     const SxVecRef<SxComplex16> &psi) const
{
   SX_CHECK (pBasis);
   SX_CHECK (pIn == pBasis.getPtr ());
   SX_CHECK (psi.getNRows () == getNElements (),
             psi.getNRows (), getNElements ());
   SxVecRef<SxComplex16, SubMatrix> res
      = psi.getRef<SubMatrix> (gBasis->ng, pBasis->getNElements (),
                               psi.getNCols ());
   res.auxData = psi.auxData;
   res.setBasis (pIn);
   return res;
}

SxVector<PrecCoeffG>
SxPAWBasis::toRBasis (const SxRBasis *rIn,
                      const SxVecRef<SxComplex16> &psi) const
{
   SX_CHECK (psi.getSize () > 0);
   ssize_t N = getNElements (), ng = gBasis->ng;
   ssize_t nPsi = psi.getNCols ();
   if (nPsi == 1)  {
      PsiG res = gBasis->toRealSpace (rIn, psi.getRef(ng));
      res.auxData = psi.auxData;
      res.setBasis (rIn);
      return res;
   }
   // --- rare case: multi-state projection
   PsiG res (rIn->getNElements (), nPsi);
   res.auxData = psi.auxData;
   res.setBasis (rIn);
   for (ssize_t i = 0; i < nPsi; ++i)  {
      res.colRef (i) <<= gBasis->toRealSpace (rIn, psi.getRef(ng, N * i));
   }
   return res;
}


SxComplex16 SxPAWBasis::scalarProduct (const SxVecRef<SxComplex16> &x,
                                       const SxVecRef<SxComplex16> &y) const
{
   SX_CHECK (x.getSize () == getNElements (),
             x.getSize (), getNElements ());
   ssize_t ng = gBasis->ng;
   if (&x == &y) return x.getRef (ng).normSqr ();
   return dot (x.getRef (ng), y.getRef (ng));
}

SxVector<SxComplex16> SxPAWBasis::toAO (const SxAOBasis *basis,
                                        const SxVecRef<SxComplex16> &in) const
{
   int ik = in.auxData.ik;
#ifndef NDEBUG
   int nk = (int)basis->refOrbitals.getSize ();
   SX_CHECK (ik >= 0 && ik < nk, ik, nk);
#endif
   if (basis->SPtr) return basis->fromPWBasis(basis->SPtr->apply(in), ik);
   return basis->fromPWBasis(in, ik);
}

PsiRef SxPAWBasis::projectFrom (const SxGBasis*,
                                const PsiRef &psi) const
{
   SX_CHECK (psi.getBasisPtr () == gBasis.getPtr ());
   SX_CHECK (psi.getSize () > 0);
   SX_CHECK (pBasis);
   SX_CHECK (pBasis->projectors);
   SxVector<PrecCoeffG> proj = pBasis->projectors->fromPWBasis (psi);

   ssize_t nPsi = psi.getNCols ();
   ssize_t N = getNElements ();
   ssize_t ng = psi.getNRows ();
   ssize_t nProj = pBasis->getNElements ();
   SxVector<SxComplex16> res (N, nPsi);
   res.auxData = psi.auxData;
   res.setBasis (this);
   // --- assemble data: G-vector coefficients, partial wave coefficients
   for (ssize_t i = 0; i < nPsi; ++i)  {
      res.getRef (ng,    i * N     ) = psi.colRef (i);
      res.getRef (nProj, i * N + ng) = proj.colRef (i);
   }
   return res;
}


SxVecRef<SxComplex16> SxPAWBasis::projectTo (const SxBasis *basis,
                                       const SxVecRef<SxComplex16> &in) const
{
   // --- identity ?
   if (dynamic_cast<const SxPAWBasis*> (basis))  {
      SX_CHECK (basis == this);
      // cast away constness
      return const_cast<SxVecRef<SxComplex16>&>(in);
   }
   // --- G basis?
   if (const SxGBasis* gk = dynamic_cast<const SxGBasis*> (basis))
      return toPWBasis (gk, in);

   if (const SxAOBasis* ao = dynamic_cast<const SxAOBasis*> (basis))
      return toAO (ao, in);

   SX_EXIT; return SxVecRef<SxComplex16>();
}
                       

