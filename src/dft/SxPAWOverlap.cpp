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

#include <SxPAWOverlap.h>
#include <SxPAWBasis.h>
#include <SxPAWHamiltonian.h>
#include <SxEigensystem.h>
#include <SxDiagMat.h>

SxComplex16 SxPAWOverlap::dot (const PsiRef &x, const PsiRef &y) const
{
   SX_CLOCK (Timer::PAWdot);
   SX_CHECK (pBasis);
   SX_START_TIMER (Timer::PAWdotP);
   SxVecRef<SxComplex16> px = *pBasis | x,
                         py = *pBasis | y;
   SX_STOP_TIMER (Timer::PAWdotP);
   return (x | y) + dotPartial (px, py);
}

SxComplex16 SxPAWOverlap::dotPartial (const SxVecRef<SxComplex16> &px,
                                      const SxVecRef<SxComplex16> &py) const
{
   const SxPAWPot &pawPot = pBasis->getPot ();
   SxComplex16 res = 0.;
   // --- partial wave correction (ref. 1, Eq. (30))
   int nProjGlobal = (int)px.getSize ();
   for (int ipg = 0; ipg < nProjGlobal; ++ipg)  {
      const SxAOBasis::OrbitalIndex &id = pBasis->projectors->orbitalMap(ipg);
      const SxArray<int>   &lPhi = pawPot.lPhi(id.is);
      const SxArray<int> &offset = pawPot.offset(id.is);
      int ipt = pBasis->projectors->refOrbMap(id.is)(id.io).n;
      const SxVector<double> &dover = pawPot.deltaS(id.is);
      int nProjLocal = int(lPhi.getSize ());
      // --- loop over projectors with same (l,m) on this atom
      for (int jpt = 0; jpt < nProjLocal; ++jpt)  {
         if (lPhi(jpt) == lPhi(ipt))  {
            // same l
            int offJ = offset(jpt) - offset(ipt);
            res += dover(ipt, jpt) * px(ipg).conj () * py(ipg + offJ);
         }
      }
   }
   return res;
}

double SxPAWOverlap::normSqr (const PsiRef &x) const
{
   SX_CHECK (pBasis);
   SX_CLOCK (Timer::PAWdot);
   SX_START_TIMER (Timer::PAWdotP);
   SxVecRef<SxComplex16> px = *pBasis | x;
   SX_STOP_TIMER (Timer::PAWdotP);
   return (x | x) + dotPartial (px, px);
}

/// Normalize states x
void SxPAWOverlap::normalize (SxVecRef<SxComplex16> *xPtr) const
{
   SX_CHECK (xPtr);
   int n = (int)xPtr->getNCols ();
   if (n == 0) n = 1;
   SX_CLOCK (Timer::PAWdot);
   SX_START_TIMER (Timer::PAWdotP);
   SxVecRef<SxComplex16> px = *pBasis | *xPtr;
   SX_STOP_TIMER (Timer::PAWdotP);

   for (int i = 0; i < n; ++i)  {
      SxVecRef<SxComplex16> Xi = xPtr->colRef(i),
                            pi = px.colRef(i);
      double nrm2 = (Xi | Xi) + dotPartial(pi, pi);
      SX_CHECK (nrm2 > 0., nrm2);
      Xi *= 1. / sqrt(nrm2);
   }
}
                         
SxVector<SxComplex16> 
SxPAWOverlap::getMatrix (const SxVecRef<SxComplex16> &x,
                         const SxVecRef<SxComplex16> &y) const
{
   SX_CHECK (x.getBasisPtr () == y.getBasisPtr ());
   SX_CHECK (pBasis);
   // --- get G basis
   const SxGBasis *gPtr = dynamic_cast<const SxGBasis *> (x.getBasisPtr ());
   if (!gPtr) gPtr = x.getBasis<SxPAWBasis> ().gBasis.getPtr ();
   SX_CHECK (gPtr);

   if (&x == &y)  {
      // --- self-overlap
      PsiRef px = *pBasis | x;
      return x.overlap (x, gPtr->ng) + getDeltaS (px, px);
   }

   if (x.getNRows () == gPtr->ng)  {
      // --- G basis
      if (x.getNCols () == y.getNCols () && x.getNCols () < 16)  {
         // --- put x and y into one matrix and apply S
         int ng = gPtr->ng, nCol = (int)x.getNCols ();
         SxIdx xPart(0, nCol * ng - 1),
               yPart(nCol * ng, 2*nCol*ng - 1);
         SxVector<PrecCoeffG> xy(ng , 2 * nCol);
         xy.setBasis (*gPtr);
         xy(xPart) <<= x;
         xy(yPart) <<= y;
         SxVector<PrecCoeffG> pxy = std::move(*pBasis | xy);
         return x.overlap (y) 
                + getDeltaS (pxy(xPart).reshape (ng, nCol),
                             pxy(yPart).reshape (ng, nCol));
      } else if (x.getNCols () > y.getNCols ())  {
         // apply S to y
         return x.overlap (apply(y));
      } else {
         // apply S to x
         return apply (x).overlap (y);
      }
   }
   // PAW basis
   return x.overlap (y, gPtr->ng) + getDeltaS (*pBasis | x, *pBasis | y);
}

SxVector<PrecCoeffG> 
SxPAWOverlap::getDeltaS (const PsiRef &px, const PsiRef &py) const
{
   const SxPAWPot &pawPot = pBasis->getPot ();
   SxVector<SxComplex16> Sni(px.getNCols (), px.getNRows ());
   Sni.set (0.);

   int offset = 0;
   int nOrb = int(pBasis->getNElements ());
   for (int is = 0; is < pawPot.getNSpecies (); ++is)  {
      int nProjLocal = pawPot.getNProj (is);

      // get 1-center corrections
      const SxVector<double> &dS = pawPot.deltaS(is);

      // --- loop over atoms
      for (int ia = 0; 
           offset < nOrb && pBasis->projectors->orbitalMap(offset).is == is;
           ++ia, offset+=nProjLocal)
      {

         // --- loops over projectors
         for (int ipl = 0; ipl < nProjLocal; ++ipl)  {
            const SxAOBasis::AoIndex &idI 
               = pBasis->projectors->refOrbMap(is)(ipl);
            for (int jpl = 0; jpl < nProjLocal; ++jpl)  {
               const SxAOBasis::AoIndex &idJ 
                  = pBasis->projectors->refOrbMap(is)(jpl);
               if (idI.l != idJ.l || idI.m != idJ.m) continue;

               // --- Ref. 1, Eq. 30, |pi>(...)<pj term
               // compute sum_j S(R)ij <pj|psi'>
               // TODO:
               // pxAdj = px.adjoint ()
               // Sni.colRef (offset + ipl)
               // .plus_assign_ax (dS(..), pxAdj.colRef(offset + jpl)
               for (int iState = 0; iState < px.getNCols (); ++iState)  {
                  Sni(iState, offset + ipl) += px(offset+jpl,iState).conj () 
                                               * dS(idI.n, idJ.n);
               }
            }

         }
      }
   }
   // compute sum_j [sum_i <psi|pi> S(R)ij] <pj|psi'>]
   return Sni ^ py;
}

void SxPAWOverlap::orthonormalize (SxVecRef<PrecCoeffG> *psiPtr, 
                                   SxOrthoMethod how) const
{
   SX_CLOCK (Timer::ortho);
   SX_CHECK (psiPtr);
   PsiRef &psi = *psiPtr;

   // --- get G basis
   const SxGBasis *gPtr = dynamic_cast<const SxGBasis *> (psi.getBasisPtr ());
   if (!gPtr) gPtr = psi.getBasis<SxPAWBasis> ().gBasis.getPtr ();
   SX_CHECK (gPtr);

   // --- compute current overlap matrix S
   SxVector<PrecCoeffG> S = getMatrix (psi, psi);

   if (how == GramSchmidt)  {
      // --- Gram/Schmidt orthogonalization
      SxVector<PrecCoeffG> U = S.choleskyDecomposition ().adjoint ();
      psi.rotate (U.inverse ());
   }  else  {
      // --- diagonalization, according to LOEWDIN orthogonalization scheme
      SxVector<PrecCoeffG> U;  // S being the overlap matrix,
                               // U is S^-1/2, and
      SxSymEigensystem<PrecCoeffG> eig (S);
      SxDiagMat<double> Ieps(1. / sqrt(eig.vals));
      U = eig.vecs ^ Ieps ^ eig.vecs.adjoint();

      psi.rotate (U);

   }
}

void SxPAWOverlap::setOrthogonal (PsiRef *X, const PsiRef &orthoStates) const
{
   SX_CHECK (X);
   if (orthoStates.getSize () == 0) return;

   SX_CHECK (X->getNRows () == orthoStates.getNRows (),
             X->getNRows (), orthoStates.getNRows ());
   SX_CHECK (X->getBasisPtr () == orthoStates.getBasisPtr ());

   SX_CLOCK (Timer::ortho);

   PsiRef &psi = *X;
   int nPsi = (int)psi.getNCols ();
   if (nPsi == 0) nPsi = 1;
   int nOrtho = (int)orthoStates.getNCols ();

   // for low n, the direct (convential) BLAS1 approach is faster
   // because there is a certain overhead, and BLAS2 routines may
   // not be fastest for BLAS1-like tasks
   // the cross-over point should be tested
   // this should be the zdotc / gemv crossover for nPsi == 1
   // this should be the gemv / gemm crossover for nPsi > 1
   bool blas1 = (nOrtho <= 150 && nPsi == 1);
   const SxPAWBasis *pawBasis 
      = dynamic_cast<const SxPAWBasis*> (psi.getBasisPtr ());
   if (pawBasis)  {
      if (blas1) {
         for (int j=0; j < nOrtho; j++)  {
            // psi -= |j><j|S|psi>
            PrecCoeffG scp = dot(orthoStates.colRef(j), psi);
            psi.plus_assign_ax (-scp, orthoStates.colRef(j));
         }
      } else {
         const SxGBasis &gk = *pawBasis->gBasis;
         const SxPartialWaveBasis &p = *pawBasis->pBasis;
         SxVector<PrecCoeffG> scps;
         //scps = (gk | orthoStates).adjoint () ^ (gk | psi);
         scps = orthoStates.overlap (psi, gk.ng);
         scps += getDeltaS (p | orthoStates, p | psi);
         psi -= (orthoStates ^ scps);
      }
   } else {
      // --- |G+k> basis only: do minimal amount of projections on the fly
      if (nPsi < nOrtho)  {
         // apply S operator to psi
         PsiG Spsi = apply (psi);
         if (blas1)  {
            for (int j=0; j < nOrtho; j++)  {
               const PsiRef &psiJ = orthoStates.colRef(j);
               // psi -= |j><j|S|psi>
               psi.plus_assign_ax (-::dot(psiJ, Spsi), psiJ);
            }
         } else {
            SxVector<PrecCoeffG> scps;
            scps = Spsi.overlap (orthoStates).adjoint ();
            psi -= (orthoStates ^ scps);
         }
      } else {
         // apply S operator to orthoStates
         PsiG Sortho = apply (orthoStates);
         if (blas1) {
            for (int j=0; j < nOrtho; j++)  {
               // psi -= |j><j|S|psi>
               PrecCoeffG scp = ::dot(Sortho.colRef(j), psi);
               psi.plus_assign_ax (-scp, orthoStates.colRef(j));
            }
         } else {
             SxVector<PrecCoeffG> scps;
             scps = Sortho.overlap (psi);
             psi -= (orthoStates ^ scps);
         }
      }
   }
}

SxVecRef<SxComplex16> SxPAWOverlap::apply (const PsiRef &psi) const
{
   const SxPAWPot &pawPot = pBasis->getPot ();
   int ik = psi.auxData.ik;
   SxVecRef<SxComplex16> p = *pBasis | psi;
   
   SxVector<SxComplex16> Sin(p.getNRows (), p.getNCols ());
   Sin.set (0.);
   int nOrb = int(pBasis->getNElements ());
   int offset = 0;
   for (int is = 0; is < pawPot.getNSpecies (); ++is)  {
      int nProjLocal = pawPot.getNProj (is);

      // get 1-center corrections
      const SxVector<double> &dS = pawPot.deltaS(is);

      // --- loop over atoms
      for (int ia = 0;
           offset < nOrb && pBasis->projectors->orbitalMap(offset).is == is;
           ++ia, offset+=nProjLocal)
      {

         // --- loops over projectors
         for (int ipl = 0; ipl < nProjLocal; ++ipl)  {
            const SxAOBasis::AoIndex &idI 
               = pBasis->projectors->refOrbMap(is)(ipl);
            for (int jpl = 0; jpl < nProjLocal; ++jpl)  {
               const SxAOBasis::AoIndex &idJ 
                  = pBasis->projectors->refOrbMap(is)(jpl);
               if (idI.l != idJ.l || idI.m != idJ.m) continue;

               // --- Ref. 1, Eq. 30, |pi>(...)<pj term
               // compute sum_j S(R)ij <pj|psi>
               for (int iState = 0; iState < p.getNCols (); ++iState)  {
                  Sin(offset+ipl, iState) += dS(idI.n, idJ.n) 
                                           * p(offset+jpl,iState);
               }
            }
         }
      }
   }
   Sin.setBasis (pBasis.getPtr ());
   Sin.auxData.ik = ik;
   const SxGBasis *gBasis
      = dynamic_cast<const SxGBasis *>(psi.getBasisPtr ());
   if (!gBasis)
      gBasis = psi.getBasis<SxPAWBasis> ().gBasis.getPtr ();
   SX_CHECK (gBasis);
   PsiG res = *gBasis | Sin;
   res += (*gBasis | psi);
   res.auxData = psi.auxData;
   res.setBasis (gBasis);
   return res;
}



