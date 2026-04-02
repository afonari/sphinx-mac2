#include <SxDFTConfig.h>
#include <SxAOBasis.h>
#include <SxProjMatrix.h>
#include <SxPWOverlap.h>
#include <SxRadialBasis.h>
#include <SxRadBasis.h>
#ifdef USE_SXGEMMM
#include <SxGemmm.h>
#endif

class SxAOBasisProj
   : public SxProjMatrix<PrecCoeffG>::SaveProjections
{
   public:
      const SxArray<SxAOBasis::OrbitalIndex> &orbitalMap;
      const SxArray<PsiG> &refOrbitals;
      mutable SxVecRef<PrecCoeffG> phase;
      SxAOBasisProj(const SxAOBasis &ao, int ik)
         : SxProjMatrix<PrecCoeffG>::SaveProjections (
               ao.getNOrb (),
               ao.blockSize,
               (int)ao.refOrbitals(ik)(ao.orbitalMap(0).is).getNRows ()),
           orbitalMap(ao.orbitalMap),
           refOrbitals(ao.refOrbitals(ik))
      { /* empty */ }

      virtual ~SxAOBasisProj () {}

      virtual void getProjector(int iOrb, SxVecRef<PrecCoeffG> *target) const
      {
         SX_CHECK (iOrb >= 0 && iOrb < orbitalMap.getSize (),
                   iOrb, orbitalMap.getSize ());
         const SxAOBasis::OrbitalIndex &idx = orbitalMap(iOrb);
         const SxVecRef<PrecCoeffG> &refOrb
            = refOrbitals(idx.is).colRef (idx.io);

         if (phase.getSize () <= 0
             || phase.auxData.is != idx.is
             || phase.auxData.ia != idx.ia)
         {
            phase.unref ();
            phase = refOrb.getBasis<SxGBasis> ().getPhaseFactors (idx.is, idx.ia);
         }
         SX_CHECK (phase.auxData.is == idx.is,
                   phase.auxData.is, idx.is);
         SX_CHECK (phase.auxData.ia == idx.ia,
                   phase.auxData.ia, idx.ia);

#ifdef USE_OPENMP
#pragma omp parallel for if (nElements > sxChunkSize)
#endif
         for (int ig = 0; ig < nElements; ++ig)
            (*target)(ig) = phase(ig) * refOrb(ig);
      }

      HAS_TARGET_GETPROJECTOR;
      NO_GETFACTOR;
};

SxAOBasis::SxAOBasis ()
{
   cacheRefOrb = Unknown;
   init ();
}

namespace Timer {
   // timer ids
   enum AOTimers{ AoRefOrbitalSetup, AoOverlap, AoSInversion, AoProjection,
                  AoGradient, AoGradProject, AoTotalTime} ;
}

SX_REGISTER_TIMERS(Timer::AOTimers)
{
   using namespace Timer;
   regTimer (AoRefOrbitalSetup,"Ref. orbital setup");
   regTimer (AoOverlap,"Overlap setup");
   regTimer (AoSInversion,"Overlap inversion");
   regTimer (AoProjection,"AO projection");
   regTimer (AoGradient,"AO gradient");
   regTimer (AoGradProject,"AO projection d/dR");
   regTimer (AoTotalTime, "AOBasis total");
}

void SxAOBasis::init ()  {
   cacheOverlap = cacheInverse = Unknown;
   refOrbCachedK = overlapCachedK = invOverlapCachedK = -1;
   blockSize = 64;
}

SxAOBasis::~SxAOBasis ()
{
   deregisterAll ();
}

SxAOBasis::SxAOBasis (const SxGkBasis &gk,
                      const SxRadBasis &rad,
                      const SxArray<SxArray<SxVector<double> > > &psiRad,
                      const SxConstPtr<SxOverlapBase> SPtrIN)
{
   init ();
   set (gk, rad, psiRad, SPtrIN);
}

SxAOBasis::SxAOBasis (const SxVecRef<SxComplex16> &YlmGl, int iSpecies)
{
   init ();
   const SxAtomicStructure &structure = YlmGl.getBasis<SxGBasis> ().getTau ();
   int nRef = (int)YlmGl.getNCols ();
   int nOrb = nRef * structure.getNAtoms (iSpecies);

   orbitalMap.resize (nOrb);
   for (int ia = 0; ia < structure.getNAtoms (iSpecies); ++ia)  {
      for (int io = 0; io < nRef; ++io)  {
         orbitalMap(io + nRef * ia) = OrbitalIndex(iSpecies, ia, io);
      }
   }

   refOrbitals.resize (1);
   refOrbitals(0).resize (structure.getNSpecies ());
   refOrbitals(0)(iSpecies) = YlmGl;
   cacheRefOrb = CacheAll;

}

void SxAOBasis::set (const SxGkBasis &gk,
                     const SxRadBasis &rad,
                     const SxArray<SxArray<SxVector<double> > > &psiRad,
                     const SxConstPtr<SxOverlapBase> SPtrIN)
{
   SX_CLOCK (Timer::AoTotalTime);
   SPtr = SPtrIN;
   if (!SPtr) SPtr = SxPtr<SxPWOverlap>::create ();
   int nSpecies = int(psiRad.getSize ());
   SX_CHECK (gk.getNk () > 0, gk.getNk ());
   SX_CHECK (nSpecies > 0, nSpecies);
   cacheRefOrb = CacheAll;

   // --- setup orbital map io -> (n,l,m) for each species
   refOrbMap.resize (nSpecies);
   SxVector<int> nAtoms(nSpecies), nOrbital(nSpecies);

   int nOrb = 0;
   {
      for (int is = 0; is < nSpecies; ++is)  {
         int nOrbType = int(psiRad(is).getSize ());
         int nOrbSpecies = 0;
         for (int iOrbType = 0; iOrbType < nOrbType; ++iOrbType)
            nOrbSpecies += 2 * psiRad(is)(iOrbType).auxData.l + 1;
         refOrbMap(is).resize (nOrbSpecies);
         int io = 0;
         for (int iOrbType = 0; iOrbType < nOrbType; ++iOrbType) {
            int l = psiRad(is)(iOrbType).auxData.l;
            for (int m = -l; m <= l; ++m, ++io)  {
               refOrbMap(is)(io) = AoIndex(iOrbType, l, m);
            }
         }
         nAtoms(is) = gk.getTau ().getNAtoms(is);
         nOrbital(is) = io;
         nOrb += io * nAtoms(is);
      }
   }

   // -- setup orbital map iOrb -> (is, ia, io)
   orbitalMap.resize (nOrb);
   {
      int iOrb = 0, ia, is, io;
      for (is = 0; is < nSpecies; ++is)
         for (ia = 0; ia < nAtoms(is); ++ia)
            for (io = 0; io < nOrbital(is); ++io, ++iOrb)
               orbitalMap(iOrb) = OrbitalIndex (is, ia, io);
   }

   // --- setup reference orbitals
   computeRefOrbitals(gk, rad, psiRad);
   // clean overlap cache
   overlap.resize (0);
   invOverlap.resize (0);
   overlapCachedK = invOverlapCachedK = -1;
}

SxAOBasis::SxAOBasis (const SxGkBasis &gk,
                      const SxArray<SxVector<double> > &psi,
                      const SxArray<SxArray<int> > &lPsi)
{
   init ();
   SX_CLOCK (Timer::AoTotalTime);
   //SPtr = SPtrIN;
   if (!SPtr) SPtr = SxPtr<SxPWOverlap>::create ();
   int nSpecies = int(psi.getSize ());
   int nk = gk.nk;
   SX_CHECK (nk > 0, nk);
   SX_CHECK (nSpecies > 0, nSpecies);
   cacheRefOrb = CacheAll;
   // clean overlap cache
   overlap.resize (0);
   invOverlap.resize (0);
   overlapCachedK = invOverlapCachedK = -1;

   // --- setup orbital map io -> (n,l,m) for each species
   refOrbMap.resize (nSpecies);
   SxVector<int> nAtoms(nSpecies), nOrbital(nSpecies);

   int nOrb = 0;
   {
      for (int is = 0; is < nSpecies; ++is)  {
         int nOrbType = (int)psi(is).getNCols ();
         int nOrbSpecies = 0;
         for (int iOrbType = 0; iOrbType < nOrbType; ++iOrbType)
            nOrbSpecies += 2 * lPsi(is)(iOrbType) + 1;
         refOrbMap(is).resize (nOrbSpecies);
         int io = 0;
         for (int iOrbType = 0; iOrbType < nOrbType; ++iOrbType) {
            int l = lPsi(is)(iOrbType);
            for (int m = -l; m <= l; ++m, ++io)  {
               refOrbMap(is)(io) = AoIndex(iOrbType, l, m);
            }
         }
         nAtoms(is) = gk.getTau ().getNAtoms(is);
         nOrbital(is) = io;
         nOrb += io * nAtoms(is);
      }
   }

   // -- setup orbital map iOrb -> (is, ia, io)
   orbitalMap.resize (nOrb);
   {
      int iOrb = 0, ia, is, io;
      for (is = 0; is < nSpecies; ++is)
         for (ia = 0; ia < nAtoms(is); ++ia)
            for (io = 0; io < nOrbital(is); ++io, ++iOrb)
               orbitalMap(iOrb) = OrbitalIndex (is, ia, io);
   }

   // --- setup reference orbitals
   refOrbitals.resize (nk);
   double gMax = 0.;
   for (int ik = 0; ik < nk; ++ik)  {
      SX_MPI_LEVEL("waves-k");
      if (SxLoopMPI::myWork(ik))  {
         refOrbitals(ik).resize (nSpecies);
         for (int is = 0; is < nSpecies; ++is)  {
            refOrbitals(ik)(is).reformat (gk(ik).ng, nOrbital(is));
            refOrbitals(ik)(is).setBasis (gk(ik));
            refOrbitals(ik)(is).auxData.ik = ik;
         }
         gMax = max(gMax, sqrt(gk(ik).g2(gk(ik).ng - 1)));
      }
   }
   if (gMax == 0.) return; // no k-points

   gMax *= 1.1;
   double cutVol = FOUR_PI / 3. * gMax * gMax * gMax;
   double dg = 0.003;
   int ng1 = int(cutVol / gk.getTau ().cell.getReciprocalCell ().volume);
   int ng2 = int(gMax/dg) + 1;
   bool viaRadG =  (ng1 * nk > ng2);
   SxRadialBasis radG;
   if (viaRadG) {
      radG.set (0., gMax, ng2, SxRadialBasis::Linear);
      radG.realSpace = false;
   }

   for (int is = 0; is < nSpecies; ++is)  {

      const SxRadialBasis *myRadG
         = dynamic_cast<const SxRadialBasis*>(psi(is).getBasisPtr ());
      for (int n = 0, io = 0; n < lPsi(is).getSize (); ++n)  {
         int l = lPsi(is)(n);
         SxVecRef<double> proj
            = const_cast<SxVector<double>&>(psi(is)).colRef (n);
         proj.auxData.is = is;
         proj.auxData.ia = -1;
         proj.auxData.n  = char(n);
         proj.auxData.l  = char(l);
         //proj.setBasis (psi(is).getBasisPtr ());
         if (myRadG)  {
            const SxVecRef<double> proj_ = std::move (proj);
            proj.unref ();
            proj = myRadG->toSpline (proj_);
         } else if (viaRadG)  {
            const SxVecRef<double> proj_ = std::move (proj);
            proj.unref ();
            proj = radG.toSpline (radG | proj_);
         }
         for (int m = -l; m <= l; ++m, ++io)  {
            proj.auxData.m = char(m);
            for (int ik = 0; ik < nk; ++ik)  {
               if (refOrbitals(ik).getSize () > 0) {
                  refOrbitals(ik)(is).colRef (io) <<= gk(ik) | proj;
               }
            }
         }
      }
   }

}

void SxAOBasis::computeRefOrbitals
   (const SxGkBasis                              &gk,
    const SxRadBasis                             &rad,
    const SxArray<SxArray<SxVector<double> > > &psiRad)
{
   SX_CLOCK (Timer::AoRefOrbitalSetup);
   ssize_t nSpecies = psiRad.getSize ();
   ssize_t nk = gk.getNk ();

   double gMax = 0.;
   refOrbitals.resize (nk);
   for (int ik = 0; ik < nk; ++ik)  {
      SX_MPI_LEVEL("waves-k");
      if (SxLoopMPI::myWork(ik))  {
         refOrbitals(ik).resize(nSpecies);
         SX_LOOP(is)  {
            ssize_t nOrbital = refOrbMap(is).getSize ();
            if (nOrbital > 0)  {
               refOrbitals(ik)(is).reformat (gk(ik).ng, nOrbital);
               refOrbitals(ik)(is).setBasis (&gk(ik));
            }
         }
         gMax = max(gMax, sqrt(gk(ik).g2(gk(ik).ng - 1)));
      }
   }
   if (gMax == 0.) return; // no k-points

   double dg = 0.003;
   int nPoints = int(gMax/dg) + 1;
   SxRadialBasis gRad(0.0, gMax, nPoints, false, SxRadialBasis::Linear);

   SX_LOOP(is)  {
      SxVecRef<double> psiGRef;
      if (refOrbMap(is).getSize () == 0) continue;
      SX_LOOP(io) {
         AoIndex nlm = refOrbMap(is)(io);
         SxVecRef<double> psiRef
            = const_cast<SxVector<double>&>(psiRad(is)(nlm.n));
         psiRef.setBasis (&rad);
         psiRef.auxData.is = (int)is;
         psiRef.auxData.n = char(nlm.n);
         psiRef.auxData.l = char(nlm.l);
         psiRef.auxData.m = char(nlm.m);
         if (psiGRef.getSize () == 0
             || psiGRef.auxData.n != nlm.n
             || psiGRef.auxData.l != nlm.l)
         {
            psiGRef = gRad.toSpline (gRad | psiRef);
         }
         psiGRef.auxData.m = char(nlm.m);
         // store <G+k|mu> = sum(r) <G+k|r><r|mu>
         //(*refOrbPtr)(is).colRef(io) <<= ( g | psiRef );
         // store <G+k|mu> = sum(r,gRad) <G+k|gRad><gRad|r><r|mu>
         SX_LOOP(ik)
            if (refOrbitals(ik).getSize () > 0)
               refOrbitals(ik)(is).colRef(io) <<= gk(ik) | psiGRef;
      }
   }
}

SxVector<SxComplex16> SxAOBasis::calculateOverlap (int ik) const
{
#ifndef NDEBUG
   if (cacheRefOrb == CacheAll)  {
      int nk = int(refOrbitals.getSize ());
      SX_CHECK (nk > 0);
      SX_CHECK (ik >= 0 && ik < nk, ik, nk);
   } else if (cacheRefOrb == CacheCurrentK)  {
      SX_CHECK (refOrbCachedK == ik, refOrbCachedK, ik);
   }
#endif
   int nOrb = int(orbitalMap.getSize ());

   SX_CLOCK (Timer::AoTotalTime);
   SX_CLOCK (Timer::AoOverlap);
   SxVector<SxComplex16> ovlp;
   SxOverlap S(SPtr);
   const SxGBasis *gBasis = dynamic_cast<const SxGBasis*>
                           (refOrbitals(ik)(0).getBasisPtr ());
   SX_CHECK (gBasis);
   int ng = gBasis->ng;
   // --- high-mem
   /*
   //SxAOBasisProj aop(*this, ik);
   TPsi aoMu.reformat(ng, nOrb);
   for (iOrb = 0; iOrb < nOrb; ++iOrb)
      aoMu.colRef(iOrb) <<= S | getAOinG(ik, iOrb);
   ovlp = aop.getProjection(aoMu);
   */

   // ---medium-mem
   ovlp.reformat(nOrb, nOrb);
   int nblock = 64;
   TPsi muBlock;
   for (int iOrb = 0, ib = 0; iOrb < nOrb; /* inside */)  {
      if (ib == 0)  {
         if (iOrb + nblock > nOrb) nblock = nOrb - iOrb;
         muBlock.reformat (ng, nblock);
         muBlock.setBasis (gBasis);
         muBlock.auxData.ik = ik;
      }
      muBlock.colRef (ib) <<= getAOinG(ik, iOrb + ib);
      if (++ib == nblock)  {
         SxIdx idx(iOrb * nOrb, (iOrb + nblock) * nOrb - 1);
         ovlp(idx) <<= fromPWBasis (muBlock);
         iOrb += nblock;
         ib = 0;
      }
   }

   return ovlp;
}


SxVecRef<SxAOBasis::CoeffType>
SxAOBasis::identity  (const SxAOBasis *basis,
                      const SxVecRef<CoeffType> &psiAO) const
{
   SX_CHECK (basis == this);
   return const_cast<SxVecRef<CoeffType>&>(psiAO);
}

SxVector<PrecCoeffG>
SxAOBasis::toPWBasis (const SxGBasis *gPtr,
                      const SxVecRef<SxComplex16> &psiAO) const
{
   SX_CLOCK (Timer::AoTotalTime);
   SX_CLOCK (Timer::AoGradient);
   // no projection of empty vectors
   SX_CHECK (psiAO.getSize () > 0);

   // --- find k-point
   SX_CHECK (cacheRefOrb == CacheAll);
   int ik = psiAO.auxData.ik;
   int nk = int(refOrbitals.getSize ());
   if (ik >= 0 && ik < nk)  {
#ifndef NDEBUG
      int is = orbitalMap(0).is; // we might start at is > 0
      SX_CHECK (gPtr == refOrbitals(ik)(is).getBasisPtr ());
#endif
   } else if (cacheRefOrb == CacheAll)  {
      int is = orbitalMap(0).is; // we might start at is > 0
      SX_CHECK (is>=0, is);
      for (ik = 0; ik < nk; ++ik)  {
         if (refOrbitals(ik).getSize () == 0) continue; // different MPI task
         if (refOrbitals(ik)(is).getBasisPtr () == gPtr)
            break;
      }
      if (ik == nk)  {
         cout << "Projection from unregistered G-basis in SxAOBasis\n";
         SX_EXIT;
      }
   }


   int nPsi = (int)psiAO.getNCols ();
#ifndef USE_SXGEMMM
   /* This is the BLAS algorithm
      It collects the projectors and uses zgemm calls to do
      the data reduction.
      sxgemm3m is an alternative direct algorithm, see below.
    */
   const int SX_PWPROJ_BLOCK_MINPSI = 1;
   if (nPsi > SX_PWPROJ_BLOCK_MINPSI)  {
      // --- use projector blocking
      SxAOBasisProj aop(*this, ik);
      SxVector<SxComplex16> psiAOphase;
      if (extraPhase.getSize () > 0)  {
         const SxVector<SxComplex16> &phase = extraPhase(ik);
         SX_CHECK (phase.getSize () == getNElements (),
                   phase.getSize (), getNElements ());
         psiAOphase.reformat (getNElements (), nPsi);
         for (ssize_t iState = 0, i = 0; iState < nPsi; ++iState)
            for (ssize_t iOrb = 0; iOrb < getNElements (); ++iOrb, ++i)
               psiAOphase(i) = psiAO(i) * phase(iOrb);
      } else {
         psiAOphase = psiAO;
      }
      SxVector<PrecCoeffG> res = aop.gradient (psiAOphase);
      res.auxData = psiAO.auxData;
      res.setBasis (gPtr);
      res.auxData.ik = ik;
      return res;
   } else if (nPsi > 1)  {
      // use algorithm below for each psi
      SxVector<SxComplex16> res (gPtr->ng, nPsi);
      res.auxData = psiAO.auxData;
      res.setBasis (gPtr);
      res.auxData.ik = ik;

      for (int i = 0; i < nPsi; ++i)
         res.colRef(i) <<= toPWBasis (gPtr, psiAO.colRef(i));
      return res;
   }
#endif

   /* --- Use atom-blocking
      For a single psi, we follow this strategy:

      For
      \sum_{is,ia,ipl} T(is,ia) |p(is,ipl>  A(is, ipl, ia)
      we first compute
      |pA(is,ia)> = \sum_{ipl} |p(is,ipl)> A(is, ipl, ia)
      as a matrix operation.

      If the number of atoms becomes large (>64), we use blocking
      to limit the size of pA.

      --- Direct algorithm (sxgemmm algorithm)
      We also have a direct blocked algorithm that calculates
      p(ig,ipl) * T(ig,ia).conj () * psiAO(ipl, ia, iPsi)
      exploiting CPU vectorization and level caches.
      Its performance appears superior on modern CPUs.
   */
   int nOrb = getNOrb ();
   const SxArray<PsiG> &proj = refOrbitals(ik);

   SxVector<PrecCoeffG> res;
   res.reformat (gPtr->ng, nPsi);
   res.set (0.);
   // --- loop over species
   for (int iOrb = 0; iOrb < nOrb; /* empty */) {
      int is = orbitalMap(iOrb).is;

      // --- find number of projectors for this species
      int npl = (int)proj(is).getNCols ();

      // --- get number of atoms
      int nAtoms, iOrbStart = iOrb;
      for (nAtoms = 0; iOrb < nOrb; iOrb += npl, ++nAtoms)  {
         if (orbitalMap(iOrb).is != is) break;
         SX_CHECK (orbitalMap(iOrb).io == orbitalMap(iOrbStart).io);
      }
      iOrb = iOrbStart;

#ifndef USE_SXGEMMM
      // --- compute gradient (with blocking on atoms)
      if (npl >= nAtoms)  {
         SxVecRef<PrecCoeffG> pA;
         for (int ia = 0, iab = 0; ia < nAtoms; ++ia)  {
            if (iab == pA.getNCols ())  {
               iab = 0;
               // discard old block
               pA.unref ();
               // --- compute new pA block
               int nAtomBlock = min(nAtoms-ia, blockSize);
               SX_CHECK (nAtomBlock > 0, nAtomBlock);
               SxVecRef<SxComplex16> psiAOphase;
               if (extraPhase.getSize () > 0)  {
                  const SxVector<SxComplex16> &phase = extraPhase(ik);
                  psiAOphase = SxVector<SxComplex16> (npl, nAtomBlock);
                  for (ssize_t i = 0; i < psiAOphase.getSize (); ++i)
                     psiAOphase(i) = psiAO(iOrb + i) * phase(iOrb+i);
               } else {
                  SxIdx blockIdx(iOrb, iOrb + nAtomBlock * npl - 1);
                  psiAOphase = psiAO.getRef<Compact> (iOrb, npl, nAtomBlock);
               }
               pA = proj(is) ^ psiAOphase;
               iOrb += nAtomBlock * npl;
            }
            // apply translation
            //res += gk.getPhaseFactors (is, ia) * pA.colRef(iab++);
            int iAtom = orbitalMap(iOrbStart + npl * ia).ia;
            SxVecRef<PrecCoeffG> T     = gPtr->getPhaseFactors (is, iAtom),
                                 pAcol = pA.colRef(iab++);
#ifdef USE_OPENMP
#pragma omp parallel for if (res.getSize () > sxChunkSize)
#endif
            for (int ig = 0; ig < res.getSize (); ++ig)
               res(ig) += T(ig) * pAcol(ig);
         }
         SX_CHECK (iOrbStart + nAtoms * npl == iOrb,
                   iOrb - iOrbStart, nAtoms, npl);
      } else {
         SxVector<PrecCoeffG> phaseNl(gPtr->ng, npl);
         phaseNl.set (0.);
#else /* USE_SXGEMMM */
      {
         SxVector<double> projReal = proj(is).real ();
#endif
         SxVector<PrecCoeffG> T;
         for (int ia = 0, iab = 0; /* done inside */; ++ia)  {
            if (iab == T.getNCols ())  {
               if (iab > 0)  {
                  if (extraPhase.getSize () > 0)  {
                     const SxVector<SxComplex16> &phase = extraPhase(ik);
#ifdef USE_SXGEMMM
                     // --- multiply psiAO by orbital-dependent phase
                     ssize_t nOrbBlock = iab * npl;
                     SxVector<SxComplex16> psiAOphase (nOrbBlock, nPsi);
                     for (ssize_t iPsi = 0, iPh = 0; iPsi < nPsi; ++iPsi)  {
                        ssize_t iAO = iOrb + iPsi * nOrb;
                        for (ssize_t jOrb = 0; jOrb < nOrbBlock; ++jOrb)  {
                           psiAOphase(iPh++) = psiAO(iAO++) * phase(iOrb+jOrb);
                        }
                     }
                     // do sxgemm summation
                     sxpgemm3m(nPsi, iab, npl, gPtr->ng,
                               T.elements, projReal.elements,
                               psiAOphase.elements, npl, nOrbBlock,
                               res.elements);
#else
                     // --- multiply psiAO by orbital-dependent phase
                     //     and reorder to (iAtom, ipl)
                     SxVector<SxComplex16> psiAOphase (iab, npl);
                     for (ssize_t jab = 0; jab < iab; ++jab)
                        for (ssize_t ipl = 0; ipl < npl; ++ipl)
                           psiAOphase(jab, ipl) = psiAO(iOrb + ipl + jab * npl)
                                                * phase(iOrb + ipl + jab * npl);
                     phaseNl += T ^ psiAOphase;
#endif
                  } else {
#ifdef USE_SXGEMMM
                     sxpgemm3m(nPsi, iab, npl, gPtr->ng,
                               T.elements, projReal.elements,
                               psiAO.elements + iOrb, npl, nOrb,
                               res.elements);
#else
                     // sum over atoms
                     phaseNl += T ^ psiAO.getRef<Compact> (iOrb, npl, iab, 1)
                                                   .transpose ();
#endif
                  }
                  iOrb += iab * npl;
               }
               if (ia == nAtoms) break;
               iab = 0;
               int nAtomBlock = min(nAtoms-ia, blockSize);
               SX_CHECK (nAtomBlock > 0, nAtomBlock);
               T.reformat (gPtr->ng, nAtomBlock);
            }
            // collect phase factors
            int iAtom = orbitalMap(iOrbStart + npl * ia).ia;
            T.colRef (iab++) <<= gPtr->getPhaseFactors(is, iAtom);
         }

#ifndef USE_SXGEMMM
         // --- sum over local projectors
         int ng = gPtr->ng;
         /*
         for (int ipl = 0; ipl < npl; ++ipl)
            res += proj(is).colRef(ipl) * phaseNl.colRef(ipl);
         */
         for (int ipl = 0; ipl < npl; ++ipl)  {
            const SxVecRef<PrecCoeffG> &phiRef = proj(is).colRef(ipl);
#ifdef USE_OPENMP
#pragma omp parallel for if (ng > sxChunkSize)
#endif
            for (int ig = 0; ig < ng; ig++)
               res(ig) += phiRef(ig) * phaseNl(ig + ng * ipl);
         }
         /*
         int gBlock = 64;
         int ig = 0, ig0;
         for (int igstop = 0; ; igstop += gBlock)  {
            if (igstop > ng) igstop = ng;
            ig0 = ig;
            for (int ipl = 0; ipl < npl; ++ipl)  {
               const SxVecRef<PrecCoeffG> &phiRef = proj(is)(ipl);
               PrecCoeffG *srcPtr = &phaseNl(ig0, ipl);
               for (ig = ig0; ig < igstop; ++ig)
                  //res(ig) += phiRef(ig) * phaseNl(ig, ipl);
                  res(ig) += phiRef(ig) * *srcPtr++;
            }
            if (ig == ng) break;
         }
         */
#endif
      }
   }

   SX_VALIDATE_VECTOR(res);

   res.auxData = psiAO.auxData;
   res.auxData.ik = ik;
   res.setBasis (gPtr);
   return res;
}

SxAOBasis::TPsi
SxAOBasis::fromPWBasis (const SxVecRef<SxComplex16> &psiG) const
{
   // no projection of empty vectors
   SX_CHECK (psiG.getSize () > 0);

   // --- find k-point
   SX_CHECK (cacheRefOrb == CacheAll);
   int ik = psiG.auxData.ik;
   int nk = int(refOrbitals.getSize ());
   int is = orbitalMap(0).is; // first species, may be > 0
   if (ik >= 0 && ik < nk)  {
      SX_CHECK (psiG.getBasisPtr () == refOrbitals(ik)(is).getBasisPtr ());
   } else if (cacheRefOrb == CacheAll)  {
      SX_CHECK (is>=0, is);
      for (ik = 0; ik < nk; ++ik)  {
         if (refOrbitals(ik)(is).getBasisPtr () == psiG.getBasisPtr ())
            break;
      }
      if (ik == nk)  {
         cout << "Projection from unregistered G-basis in SxAOBasis\n";
         SX_EXIT;
      }
   }
   return SPtr ? fromPWBasis (SPtr->apply(psiG), ik)
               : fromPWBasis (psiG             , ik);
}

SxAOBasis::TPsi
SxAOBasis::fromPWBasis (const SxVecRef<PrecCoeffG> &psiS, int ik) const
{
   SX_CLOCK (Timer::AoTotalTime);
   SX_CLOCK (Timer::AoProjection);

#ifndef USE_SXGEMMM
   /* This is the BLAS algorithm
      It collects the projectors and uses zgemm calls to do
      the data reduction.
      sxgemmm is an alternative direct algorithm, see below.
    */
   const int SX_PWPROJ_BLOCK_MINPSI = 1;
   if (psiS.getNCols () > SX_PWPROJ_BLOCK_MINPSI)  {
      // --- use projector blocking
      SxAOBasisProj aop(*this, ik);

      TPsi result = aop.getProjectionFromExtended(psiS);

      SX_VALIDATE_VECTOR (result);
      result.auxData = psiS.auxData;
      result.setBasis (this);
      if (extraPhase.getSize () > 0)  {
         const SxVector<SxComplex16> &phase = extraPhase(ik);
         SX_CHECK (phase.getSize () == getNElements (),
                   phase.getSize (), getNElements ());
         for (ssize_t iState = 0, iRes = 0; iState < result.getNCols (); ++iState)
            for (ssize_t iOrb = 0; iOrb < getNElements (); ++iOrb, ++iRes)
               result(iRes) *= phase(iOrb).conj ();
      }
      return result;
   } else if (psiS.getNCols () > 1)  {
      int nPsi = (int)psiS.getNCols ();
      // use algorithm below for each psi
      SxVector<SxComplex16> res (getNOrb (), nPsi);

      for (int i = 0; i < nPsi; ++i)
         res.colRef(i) <<= fromPWBasis (psiS.colRef(i));
      res.auxData = psiS.auxData;
      res.setBasis (this);
      return res;
   }
#endif

   /* --- Use atom-blocking (BLAS algorithm)
      For a single psi, we follow this strategy:

      First note that
      <p(is,ipl) * T(is,ia)|psi> = <p(is,ipl) | T*(is,ia) * psi>

      We now perform the right-hand formula as matrix operation for
      each species:
      <T * p|psi> (npl * nAtoms) = p(npl x ng) ^ psiT(ng x nAtoms)

      If the number of atoms becomes large (>64), we use blocking
      to limit the size of psiT.

      --- Direct algorithm (sxgemmm algorithm)
      We also have a direct blocked algorithm that calculates
      p(ig,ipl) * T(ig,ia).conj () * psi(ig, iState)
      exploiting CPU vectorization and level caches.
      Its performance appears superior on modern CPUs.
      
   */
   const SxArray<PsiG> &proj = refOrbitals(ik);
   const SxGBasis *gk = dynamic_cast<const SxGBasis*>(psiS.getBasisPtr ());
   int nOrb = getNOrb ();
   ssize_t nStates = psiS.getNCols ();
   TPsi res(nOrb * nStates);
   res.reshape (nOrb, nStates);
   // --- loop over species
   for (int iOrb = 0; iOrb < nOrb; /* empty */) {
      int is = orbitalMap(iOrb).is;
      if (!gk) {
         // get G basis from PAW basis
         gk = dynamic_cast<const SxGBasis*> (proj(is).getBasisPtr ());
         SX_CHECK (gk);
      }
      SX_CHECK (gk == proj(is).getBasisPtr ());

      // number of projectors for this species
      int npl = (int)proj(is).getNCols ();

      // --- get number of atoms
      int nAtoms, iOrbStart = iOrb;
      for (nAtoms = 0; iOrb < nOrb; iOrb += npl, ++nAtoms)  {
         if (orbitalMap(iOrb).is != is) break;
         SX_CHECK (orbitalMap(iOrb).io == orbitalMap(iOrbStart).io);
      }
      iOrb = iOrbStart;

#ifndef USE_SXGEMMM
      // --- shortcut for nAtoms = 1
      if (nAtoms == 1)  {
         PsiG psiTrans = gk->getPhaseFactors(is, orbitalMap(iOrb).ia).conj()
                       * psiS;
         res (SxIdx(iOrb,iOrb+npl-1)) <<= proj(is).overlap (psiTrans);
         iOrb += npl;
         continue;
      }
#else
      SxVector<double> projReal = proj(is).real ();
#endif

      // --- perform projections (with blocking on atoms)
      SxVector<PrecCoeffG> atomBlock;
      int ng = gk->ng;
      for (int ia = 0, iab = 0; ia < nAtoms; ++ia)  {
         if (iab == 0)  {
            // start new (psi * T) block
            atomBlock.reformat (gk->ng, min(nAtoms-ia, blockSize));
         }
         {
            // collect into atomBlock:
            // - for BLAS-based algorithm: psi * T.conj ()
            // - for sxgemmm-based case: T.conj ();
            //psiTrans.colRef (iab++)<<= psi * gk.getPhaseFactors(is, ia).conj ();
            // get the real atom index
            int iAtom = orbitalMap(iOrbStart + npl * ia).ia;
            SxVecRef<PrecCoeffG> bcol = atomBlock.colRef(iab), T = gk->getPhaseFactors(is, iAtom);
#ifdef USE_OPENMP
#pragma omp parallel for if (ng > sxChunkSize)
#endif
            for (int ig = 0; ig < ng; ig++)
#ifdef USE_SXGEMMM
               bcol(ig) = T(ig).conj ();
#else
               bcol(ig) = psiS(ig) * T(ig).conj ();
#endif
            iab++;
         }

         if (iab == atomBlock.getNCols ())  {
            // project current block
            SX_CHECK (proj(is).imag ().normSqr () < 1e-20);
#ifdef USE_SXGEMMM
            sxgemmm(nStates, iab, npl, ng,
                    psiS.elements, psiS.getNRows (), atomBlock.elements,
                    projReal.elements, res.elements + iOrb, npl, nOrb);
#else
            SxIdx currentOrbs(iOrb, iOrb + iab*npl - 1);
            res(currentOrbs) <<= proj(is).overlap (atomBlock);
#endif
            iOrb += iab * npl;
            iab = 0;
         }
      }
      SX_CHECK (iOrbStart + nAtoms * npl == iOrb,
                iOrb - iOrbStart, nAtoms, npl);
   }
   SX_VALIDATE_VECTOR(res);
   res.auxData = psiS.auxData;
   res.setBasis (this);
   if (extraPhase.getSize () > 0)  {
      const SxVector<SxComplex16> &phase = extraPhase(ik);
      SX_CHECK (phase.getSize () == getNElements (),
                phase.getSize (), getNElements ());
      for (ssize_t iState = 0, iRes = 0; iState < res.getNCols (); ++iState)
         for (ssize_t iOrb = 0; iOrb < getNElements (); ++iOrb, ++iRes)
            res(iRes) *= phase(iOrb).conj ();
   }
   return res;
}

int SxAOBasis::getLMax () const
{
   int result = 0;
   int nOrbs = getNOrb ();
   for (int iOrb = 0; iOrb < nOrbs; iOrb++)  {
      int is = orbitalMap(iOrb).is;
      int io = orbitalMap(iOrb).io;
      int l = refOrbMap(is)(io).l;
      if (result < l) result = l;
   }

   return result;
}

PsiG SxAOBasis::getAOinG (int ik) const
{
   SX_CHECK (cacheRefOrb == CacheAll);
   if (cacheRefOrb == CacheCurrentK)  {
      SX_CHECK (refOrbCachedK == ik, refOrbCachedK, ik);
      ik = 0;
   }
   SX_CHECK (ik >= 0 && ik < refOrbitals.getSize (),
             ik, refOrbitals.getSize ());
   const SxGBasis *gPtr 
      = dynamic_cast<const SxGBasis *> (refOrbitals(ik)(0).getBasisPtr());
   ssize_t ng = gPtr->g2.getSize ();
   ssize_t nOrbitals = orbitalMap.getSize ();
   PsiG result(ng,nOrbitals);
   for (int iOrb = 0; iOrb < nOrbitals; iOrb++)  {
     result.colRef(iOrb) <<= getAOinG (ik,iOrb);
   }
   result.setBasis(gPtr);
   result.auxData.ik = ik;

   return result; 
}

PsiG SxAOBasis::getAOinG (int ik, int iOrb) const
{
   SX_CHECK (cacheRefOrb == CacheAll);
   int trueK = ik;
   if (cacheRefOrb == CacheCurrentK)  {
      SX_CHECK (refOrbCachedK == ik, refOrbCachedK, ik);
      ik = 0;
   }
   SX_CHECK (ik >= 0 && ik < refOrbitals.getSize (),
             ik, refOrbitals.getSize ());
   SX_CHECK (iOrb >= 0 && iOrb < orbitalMap.getSize (),
             iOrb, orbitalMap.getSize ());
   const OrbitalIndex &idx = orbitalMap(iOrb);
   const PsiG &refOrb = refOrbitals(ik)(idx.is);

   const SxGBasis *gPtr
      = dynamic_cast<const SxGBasis *> (refOrb.getBasisPtr ());
   SX_CHECK (gPtr);
   PsiG res;
   res = gPtr->getPhaseFactors(idx.is,idx.ia) * refOrb.colRef(idx.io);
   if (extraPhase.getSize () > 0) res *= extraPhase(ik)(iOrb);
   const AoIndex &nlm = refOrbMap(idx.is)(idx.io);
   res.auxData.ik = trueK;
   res.auxData.is = idx.is;
   res.auxData.ia = idx.ia;
   res.auxData.n = char(nlm.n);
   res.auxData.l = char(nlm.l);
   res.auxData.m = char(nlm.m);

   SX_VALIDATE_VECTOR (res);
   return res;
}

SxVecRef<SxComplex16> SxAOBasis::getOverlap (int ik) const
{
   // --- single k caching
   if (cacheOverlap == CacheCurrentK)  {
      if (ik != overlapCachedK)  {
         overlapCachedK = ik;
         overlap(0).unref ();
         overlap(0) = calculateOverlap(ik);
      }
      return overlap(0);
   }
   // --- all k caching
   if (cacheOverlap == CacheAll)  {
      SX_CHECK (ik >= 0 && ik < overlap.getSize (),
                ik, overlap.getSize ());
      if (overlap(ik).getSize () == 0)
         overlap(ik) = calculateOverlap(ik);
      return overlap(ik);
   }
   // no caching
   SX_CHECK (cacheOverlap == Unknown || cacheOverlap == Recompute);

   if (cacheOverlap == Unknown)  {
      cout << "Warning: undefined caching behaviour for overlap matrices\n"
              "         in SxAOBasis::getOverlap. Recomputing overlap...\n";
   }

   return calculateOverlap(ik);
}

SxVecRef<SxComplex16> SxAOBasis::getInverseOverlap (int ik) const
{
   SX_CLOCK (Timer::AoTotalTime);
   // --- single k caching
   if (cacheInverse == CacheCurrentK)  {
      if (ik != invOverlapCachedK)  {
         invOverlapCachedK = ik;
         invOverlap(0).unref ();
         // get the overlap matrix
         const SxVecRef<SxComplex16> &S = getOverlap(ik);
         SX_CLOCK (Timer::AoSInversion);
         invOverlap(0) = SxVector<SxComplex16> (S).inverse ();
      }
      return invOverlap(0);
   }
   // --- all k caching
   if (cacheInverse == CacheAll)  {
      SX_CHECK (ik >= 0 && ik < invOverlap.getSize (),
                ik, overlap.getSize ());
      if (invOverlap(ik).getSize () == 0)  {
         SX_CLOCK (Timer::AoSInversion);
         invOverlap(ik) = SxVector<SxComplex16> (getOverlap(ik)).inverse ();
      }
      return invOverlap(ik);
   }
   // no caching
   SX_CHECK (cacheInverse == Unknown || cacheInverse == Recompute);

   if (cacheInverse == Unknown)  {
      cout << "Warning: undefined caching behaviour for overlap matrices\n"
              "         in SxAOBasis::getOverlap. Recomputing overlap...\n";
   }

   SxVector<SxComplex16> res = getOverlap(ik);
   SX_CLOCK (Timer::AoSInversion);
   return std::move(res).inverse ();
}

void SxAOBasis::setOverlapCaching (enum Caching mode, int nk)
{
   SX_CHECK (mode == Recompute ||
             mode == CacheCurrentK ||
             mode == CacheAll);
   if (mode == Recompute)  {
      overlap.resize (0);
      overlapCachedK = -1;
   }
   if (mode == CacheCurrentK)  {
      overlap.resize (1);
      overlapCachedK = -1;
   }
   if (mode == CacheAll)  {
      if (nk <= 0) nk = int(refOrbitals.getSize ());
      if (cacheOverlap == CacheCurrentK && overlapCachedK != -1)  {
         SxVecRef<SxComplex16> ovlp = overlap(0);
         overlap.resize (nk);
         overlap(overlapCachedK).unref (); // needed for nk=1
         overlap(overlapCachedK) = ovlp;
      }
      overlapCachedK = -1;
   }
   cacheOverlap = mode;
}

void SxAOBasis::setInvOverlapCaching (enum Caching mode, int nk)
{
   SX_CHECK (mode == Recompute ||
             mode == CacheCurrentK ||
             mode == CacheAll);
   if (mode == Recompute)  {
      invOverlap.resize (0);
      invOverlapCachedK = -1;
   }
   if (mode == CacheCurrentK)  {
      invOverlap.resize (1);
      invOverlapCachedK = -1;
   }
   if (mode == CacheAll)  {
      if (nk <= 0) nk = int(refOrbitals.getSize ());
      if (cacheInverse == CacheCurrentK && invOverlapCachedK != -1)  {
         SxVecRef<SxComplex16> ovlp = overlap(0);
         invOverlap.resize (nk);
         invOverlap(invOverlapCachedK).unref (); // needed for nk=1
         invOverlap(invOverlapCachedK) = ovlp;
      } else {
         overlap.resize (nk);
      }
      invOverlapCachedK = -1;
   }
   cacheInverse = mode;
}

class SxAOBasisProjGrad
: public SxProjMatrix<PrecCoeffG>::SaveProjections
{
   public:
      const SxArray<SxAOBasis::OrbitalIndex> &orbitalMap;
      const SxArray<PsiG> &refOrbitals;
      mutable SxVector<PrecCoeffG> T; // cached phase factors ("atom translation")

      const SxGBasis &gk;

      SxAOBasisProjGrad (const SxAOBasis &ao, int ik)
         : SxProjMatrix<PrecCoeffG>::SaveProjections (
              3 * ao.getNOrb (),
              ao.blockSize,
              (int)ao.refOrbitals(ik)(ao.orbitalMap(0).is).getNRows ()),
         orbitalMap(ao.orbitalMap),
         refOrbitals(ao.refOrbitals(ik)),
         gk(refOrbitals(ao.orbitalMap(0).is).getBasis<SxGBasis> ())
      {
         // empty
      }
      virtual ~SxAOBasisProjGrad () {}

      HAS_TARGET_GETPROJECTOR;
      NO_GETFACTOR;
      virtual void getProjector(int i, SxVecRef<PrecCoeffG> *target) const
      {
         int iOrb = i / 3;
         int iDir = i - 3 * iOrb;
         SX_CHECK (iOrb >= 0 && iOrb < orbitalMap.getSize (),
                   iOrb, orbitalMap.getSize ());
         const SxAOBasis::OrbitalIndex &idx = orbitalMap(iOrb);
         const PsiRef &refOrb = refOrbitals(idx.is).colRef (idx.io);

         if (T.getSize () <= 0
             || T.auxData.is != idx.is
             || T.auxData.ia != idx.ia)
         {
            SxVecRef<PrecCoeffG> phi = refOrb.getBasis<SxGBasis> ().getPhaseFactors (idx.is, idx.ia);
            T.reformat (nElements, 3);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for (int jDir = 0; jDir < 3; ++jDir)  {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
               for (int ig = 0; ig < nElements; ++ig)  {
                  T(ig + jDir * nElements) = I * phi(ig) 
                                           * gk.gVec(ig + jDir * nElements);
               }
            }
            T.auxData.is = idx.is;
            T.auxData.ia = idx.ia;
         }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
         for (int ig = 0; ig < nElements; ++ig)
            (*target)(ig) =  T(ig + iDir * nElements) * refOrb(ig);
      }
};

SxArray<SxVector<SxComplex16> >
SxAOBasis::gradProject (const PsiRef &psi) const
{
   //SX_CLOCK (Timer::GradProject);
   SX_CHECK(psi.getSize () >0);

   // --- find k-point
   SX_CHECK (cacheRefOrb == CacheAll);
   int ik = psi.auxData.ik;
   int nk = int(refOrbitals.getSize ());
   PsiRef psiS = SPtr->apply (psi);
   if (ik >= 0 && ik < nk)  {
      SX_CHECK (psiS.getBasisPtr () == refOrbitals(ik)(0).getBasisPtr ());
   } else if (cacheRefOrb == CacheAll)  {
      for (ik = 0; ik < nk; ++ik)  {
         if (refOrbitals(ik)(0).getBasisPtr () == psiS.getBasisPtr ())
            break;
      }
      if (ik == nk)  {
         cout << "Projection from unregistered G-basis in SxAOBasis\n";
         SX_EXIT;
      }
   }
   return gradProject(psiS, ik);
}

SxArray<SxVector<SxComplex16> >
SxAOBasis::gradProject (const PsiRef &psiS, int ik) const
{
   SX_CLOCK (Timer::AoTotalTime);
   SX_CLOCK (Timer::AoGradProject);
   int nPsi = (int)psiS.getNCols ();
   if (nPsi==0) nPsi = 1;
   int nOrb = getNOrb ();

   SxArray<SxVector<SxComplex16> > result(3);
   SxVecRef<SxComplex16> allProj;
   allProj = SxAOBasisProjGrad(*this, ik).getProjectionFromExtended (psiS);
   allProj.reshape (3, nOrb * nPsi);
   for (int iDir = 0; iDir < 3; ++iDir)  {
      result(iDir) = allProj.rowRef (iDir); // performs a copy!
      result(iDir).reshape (nOrb, nPsi);
      if (extraPhase.getSize () > 0)  {
         const SxVector<SxComplex16> &phase = extraPhase(ik);
         SX_CHECK (phase.getSize () == getNElements (),
                   phase.getSize (), getNElements ());
         for (ssize_t iState = 0, iRes = 0; iState < nPsi; ++iState)
            for (ssize_t iOrb = 0; iOrb < getNElements (); ++iOrb, ++iRes)
               result(iDir)(iRes) *= phase(iOrb).conj ();
      }
   }
   return result;
}



