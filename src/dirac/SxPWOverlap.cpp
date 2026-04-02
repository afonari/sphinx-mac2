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

#include <SxPWOverlap.h>
#include <SxTypes.h>
#include <SxTimer.h>
#include <SxEigensystem.h>
#include <SxDiagMat.h>

double SxPWOverlap::normSqr (const SxVecRef<SxComplex16> &psi) const
{
   return psi.normSqr ();
}

SxComplex16 SxPWOverlap::dot (const SxVecRef<SxComplex16> &x,
                              const SxVecRef<SxComplex16> &y) const
{
   return ::dot (x, y);
}

/// Apply overlap operator
SxVecRef<SxComplex16> 
SxPWOverlap::apply (const SxVecRef<SxComplex16> &psi) const
{
   // nothing to be done
   return const_cast<SxVecRef<SxComplex16>&>(psi);
}

/// Set states x orthogonal to states y
void SxPWOverlap::setOrthogonal (SxVecRef<PrecCoeffG> *X, const PsiRef &orthoStates) const
{
   // by C. Freysoldt
   // This routine is heavily commented on technical details
   // This is because I'm implementing it now, but the XPress library
   // will need changes here, but maybe not by me
   SX_CLOCK (Timer::ortho);

   SX_CHECK (X);
   if (orthoStates.getSize () == 0) return;

   SX_CHECK (X->getNRows () == orthoStates.getNRows (),
             X->getNRows (), orthoStates.getNRows ());

   PsiRef &psi = *X;
   int nPsi = (int)psi.getNCols ();
   if (nPsi == 0) nPsi = 1;
   int nOrtho = (int)orthoStates.getNCols ();

   int j;
   PrecCoeffG scp;
   if (nOrtho > 150 || nPsi > 1)  {
      // for low n, the direct (convential) BLAS1 approach is faster
      // because there is a certain overhead, and BLAS2 routines may
      // not be fastet for BLAS1-like tasks
      // the cross-over point should be tested
      // this should be the zdotc / gemv(transpose='c') crossover for nPsi == 1
      // this should be the gemv / gemm crossover for nPsi > 1

      // --- get <j|psi> for all 0 <= j < nOrtho
      SxVector<PrecCoeffG> scps = orthoStates.overlap (psi);
      
      // --- subtract |j><j|psi>
      // BLAS2 (nPsi == 1) or BLAS3 version
      psi -= (orthoStates ^ scps);
      
   } else {
      // --- conventional BLAS1 approach
      for (j=0; j < nOrtho; j++)  {
         const PsiRef &psiJ = orthoStates.colRef(j);
//       psi -= psiJ * scp;
         psi.plus_assign_ax (-::dot(psiJ, psi), psiJ); // psi -= |j><j|psi>
      }
   }
}

/// Orthonormalize states x
void SxPWOverlap::orthonormalize (SxVecRef<PrecCoeffG> *psiPtr, 
                                  SxOrthoMethod how) const
{
   SX_CHECK (psiPtr);
   PsiRef &psi = *psiPtr;

   if (how == GramSchmidt)  {
      // --- Gram/Schmidt orthogonalization
      SX_CLOCK (Timer::ortho);
      SxVector<PrecCoeffG> S = psi.overlap (psi);
      SxVector<PrecCoeffG> U = S.choleskyDecomposition ().adjoint ();
      psi.rotate (U.inverse ());
   }  else  {
   // --- diagonalization, according to LOEWDIN orthogonalization scheme
      SX_CLOCK (Timer::ortho);
      SxVector<PrecCoeffG> S, U; // S being the overlap matrix,
                                 // U is U^-1/2, and
      {
         SX_CLOCK (Timer::loewdinS);
         S = psi.overlap (psi);
      }

      SxSymEigensystem<PrecCoeffG> eig(S);
      SxDiagMat<double> Ieps(1. / sqrt(eig.vals));
      U = eig.vecs ^ Ieps ^ eig.vecs.adjoint();

      psi.rotate (U);

   }
}
