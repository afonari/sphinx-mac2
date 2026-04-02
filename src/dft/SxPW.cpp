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

#include <SxPW.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <SxMathLib.h>
//#include <SFHIngX.h>
#include <SxLoopMPI.h>  // LoopMPI
#include <SxAllocCache.h>
#include <SxPWOverlap.h>
#include <SxEigensystem.h>
#include <SxDiagMat.h>

#ifndef WIN32
#   include <unistd.h>
#endif /* WIN32 */

//Ref1: Rev. Mod. Phys. (64), 1045-1097 (1992)

SxPW::SxPW () : SxPWSet ()
{
   nStates = nSpin = nComp = nkPoints = 0;
   loadedK = -1;
   loadedSpin = -1;
   loadedI = -1;
   stored = Unknown;
   registerMemoryObservers ();
#ifdef SXPW_FINDLOOP
   findLoopLastIk = findLoopLastIspin = -1;
#endif
}


SxPW::SxPW (int nStatesIn, int nSpinIn, SxPtr<SxGkBasis> gkBasisPtrIn,
            int nCompIn,
            const SxString &tmpDirIn)
   : SxPWSet ()
{
   SX_CHECK (gkBasisPtrIn.getPtr () != NULL);
   SX_CHECK (gkBasisPtrIn->gBasisList(0));
   nStates    = nStatesIn;
   nSpin      = nSpinIn;
   tmpDir     = tmpDirIn;
   gkBasisPtr = gkBasisPtrIn;
   int ng, ik, nk, iSpin;
   nk = gkBasisPtr->nk;
   nComp = nCompIn;
   SX_CHECK ((char)nComp == nComp, nComp); // must be within char range
   SxGkBasis &gkBasis = *gkBasisPtr;

   SX_CHECK  (nk > 0, nk);

   nkPoints   = nk;
   if ( tmpDir.getSize () == 0 )  {
      // --- keep waves in memory
      stored = InMemory;
      waves.resize (nk);
      for (ik=0; ik < nk; ik++)
      {
         ng = gkBasis(ik).ng;
         waves(ik).resize (nSpin);
         nStatesPerK.append (nStates);
         nGIPerK.append (nStates*ng);
         for (iSpin=0; iSpin < nSpin; iSpin++)
         {
            SX_MPI_LEVEL("waves-k");
            if (SxLoopMPI::myWork(ik))
            {  // LoopMPI
               waves(ik)(iSpin).reformat (ng * nComp, nStates);
               waves(ik)(iSpin).auxData.nComp = (char)nComp;
               waves(ik)(iSpin).setBasis (&gkBasis(ik));
            }  // LoopMPI
         }
      }
   }  else  {
      // --- keep waves on disk
      stored = KeepOnDisk;
      createScratchFile ();
      loadWaves (0,0);
   }
   registerMemoryObservers ();
#ifdef SXPW_FINDLOOP
   findLoopLastIk = findLoopLastIspin = -1;
#endif
}

SxPW::SxPW (const SxPW &in)
   : SxPWSet ()
{
   this->gkBasisPtr = in.gkBasisPtr;
   SxGkBasis &gkBasis = *gkBasisPtr;
   int ik, nk, iSpin, nGI, nRows, nCols;
   nk = in.getNk();
   nSpin = in.getNSpin ();
   nComp = in.nComp;
   stored = in.stored;
   tmpDir     = in.tmpDir;
   if (stored != InMemory)  {
      SX_EXIT;
      stored = KeepOnDisk;
      createScratchFile ();
   }

   nkPoints   = nk;
   nStates = in.nStates;
   if ( stored == InMemory )  {
      waves.resize (nk);
      for (ik=0; ik < nk; ik++)  {
         waves(ik).resize (nSpin);
         nStatesPerK.append (in.getNStates (ik));
         if (nStates >= 0)  {
            if (nStates != nStatesPerK.last ())
               nStates = -1;
         }
         nGIPerK.append ((int)in.waves(ik)(0).getSize());
         for (iSpin=0; iSpin < nSpin; iSpin++)  {
            nRows = (int)in.waves(ik)(iSpin).getNRows ();
            nCols = (int)in.waves(ik)(iSpin).getNCols ();
            nGI   = (int)in.waves(ik)(iSpin).getSize();  // {ng x i}
            waves(ik)(iSpin).resize (nGI);
            waves(ik)(iSpin).reshape (nRows, nCols);
            waves(ik)(iSpin) <<= in.waves(ik)(iSpin);
            waves(ik)(iSpin).setBasis( &gkBasis(ik) );
         }
      }
   }  else  {  // keep waves on disk
      SX_EXIT;
   }

   registerMemoryObservers ();
#ifdef SXPW_FINDLOOP
   findLoopLastIk = findLoopLastIspin = -1;
#endif
}


SxPW::SxPW (int nSpinIn, const SxList<int> &ngPerK,
            int nCompIn,
            const SxString &tmpDirIn)
   : SxPWSet ()
{
   SX_EXIT;
   gkBasisPtr = SxPtr<SxGkBasis>::create();  // ???
   tmpDir = tmpDirIn;
   int ng, ik, nk, iSpin;
   nk = int(ngPerK.getSize());
   bool keepOnDisc = (tmpDir != "");
   nkPoints   = nk;

   if (keepOnDisc)  {
      SX_EXIT;
   }

   SX_CHECK  (nk > 0, nk);
   SX_CHECK  (nSpin == 1 || nSpin == 2, nSpin);
   nSpin = nSpinIn;
   nComp = nCompIn;
   SX_CHECK ((char)nComp == nComp, nComp); // must be within char range


   if (!keepOnDisc)  {
      stored = InMemory;
      waves.resize (nk);
      for (ik=0; ik < nk; ik++)  {
         ng = ngPerK(ik);
         waves(ik).resize (nSpin);
         nStatesPerK.append (nStates);
         nGIPerK.append (nStates*ng);
         for (iSpin=0; iSpin < nSpin; iSpin++)  {
            waves(ik)(iSpin).resize (nStates * ng * nComp);
            waves(ik)(iSpin).reshape (ng * nComp, nStates);
            waves(ik)(iSpin).auxData.nComp = (char)nComp;
         }
      }
   }  else  {
      // --- keep waves on disk
      stored = KeepOnDisk;
      createScratchFile ();
      loadWaves (0,0);
   }
   registerMemoryObservers ();
#ifdef SXPW_FINDLOOP
   findLoopLastIk = findLoopLastIspin = -1;
#endif
}


SxPW::SxPW (const SxString filename, enum StorageModel how, const SxString &tmpDirIn)
   : SxPWSet (),
     nkPoints(-1),
     nStates(-1),
     nSpin(-1),
     nComp (-1),
     stored(how)
{
   SX_CHECK (how == InMemory || how == ReadOnDemand ||
             how == ReadOneByOne || how == KeepOnDisk);
   try {
      wavesFile.open (filename, SxBinIO::BINARY_READ_ONLY);
      if (how == InMemory)  {
         read (wavesFile);
         wavesFile.close ();
      } else if (how == KeepOnDisk) {
         if (tmpDirIn.isEmpty ()) tmpDir =".";
         else tmpDir = tmpDirIn;
         stored = KeepOnDisk;
         createScratchFile ();
         read (wavesFile);
         wavesFile.close ();
      } else {
         nkPoints = wavesFile.getDimension ("nk");
         nSpin    = wavesFile.getDimension ("nSpin");
         if (wavesFile.containsDim("nComp"))
            nComp = wavesFile.getDimension ("nComp");
         else
            nComp = 1;
         SX_CHECK ((char)nComp == nComp, nComp); // must be within char range
         nStatesPerK.resize (nkPoints);
         wavesFile.read ("nPerK", &nStatesPerK, nkPoints);

         // get nGIPerK
         SxVector<int> nGk(nkPoints);
         wavesFile.read ("nGk", &nGk, nkPoints);
         nGIPerK.resize (nkPoints);
         for (int ik = 0; ik < nkPoints; ++ik)
            nGIPerK(ik) = nGk(ik) * nStatesPerK(ik);

         loadedK = -1;
         loadedSpin = -1;
         loadedI = -1;
         waves.resize(1);
         waves(0).resize (1);
         nStates = nStatesPerK(0); // ad hoc
      }
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   registerMemoryObservers ();
#ifdef SXPW_FINDLOOP
   findLoopLastIk = findLoopLastIspin = -1;
#endif
}

SxPW::~SxPW ()
{
   if (stored == KeepOnDisk)  {
      SxString tmpFilename = wavesFile.filename;
      try {
         wavesFile.close ();
      } catch (const SxException &e)  {
         e.print ();
         SX_EXIT;
      }
      unlink (tmpFilename.ascii());
   }
   if (stored == ReadOnDemand || stored == ReadOneByOne)  {
      try {
         wavesFile.close ();
      } catch (const SxException &e)  {
         e.print ();
         SX_EXIT;
      }
   }
}

void SxPW::operator= (const SxPW &in)
{
   SxPWSet::operator= (in);
   //cout << "call to SxPW::operator=" << endl; // identify hidden calls
   SX_CHECK (stored != KeepOnDisk); // would need closing of current file
   if (stored == ReadOnDemand || stored == ReadOneByOne)
      wavesFile.close ();
   this->gkBasisPtr = in.gkBasisPtr;
   nStates     = in.nStates;
   nSpin       = in.nSpin;
   nComp       = in.nComp;
   SX_CHECK ((char)nComp == nComp, nComp); // must be within char range
   nkPoints    = in.getNk();
   nStatesPerK = in.nStatesPerK;
   nGIPerK     = in.nGIPerK;

   stored     = in.stored;
   if (stored == KeepOnDisk)  {
      // take over file
      tmpDir      = in.tmpDir;
      loadedK     = in.loadedK;
      loadedSpin = in.loadedSpin;
      wavesFile.fp          = in.wavesFile.fp;
      in.wavesFile.fp       = NULL;
      wavesFile.mode        = SxBinIO::SCRATCH_READ_WRITE;
      in.wavesFile.mode     = SxBinIO::UNKNOWN;
      wavesFile.filename    = in.wavesFile.filename;
      in.wavesFile.filename = "";
      wavesFile.isOpen      = true;
      in.wavesFile.isOpen   = false;
      in.stored             = Unknown;
   }
   if (stored == ReadOnDemand)  {
      loadedK = in.loadedK;
      loadedSpin = in.loadedSpin;
      try  {
         wavesFile.open (in.wavesFile.filename, SxBinIO::BINARY_READ_ONLY);
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   }
   waves = in.waves;
   UPDATE_MEMORY (waves);
}

SxPtr<SxPWSet> SxPW::getNew () const
{
   return SxPtr<SxPW>::create (getNStates (), getNSpin (), gkBasisPtr,
                               nComp, tmpDir);
}

PsiG& SxPW::operator() (int iSpin, int ik)
{
#ifdef SXPW_FINDLOOP
   findLoopUpdate (iSpin, ik);
#endif
   SX_CHECK (stored == InMemory || stored == KeepOnDisk || stored == ReadOnDemand);
   // no ReadOneByOne (since one by one implies only one state at a time)
   if ( stored == KeepOnDisk || stored == ReadOnDemand)  {
      if (ik != loadedK || iSpin != loadedSpin)  {
         flushWaves ();
         loadWaves (ik,iSpin);
      }
      return waves(0)(0);
   }  else  {
      waves(ik)(iSpin).auxData.iSpin = char(iSpin);
      waves(ik)(iSpin).auxData.ik    = ik;
      return waves(ik)(iSpin);
   }
}


const PsiG& SxPW::operator() (int iSpin, int ik) const
{
#ifdef SXPW_FINDLOOP
   findLoopUpdate (iSpin, ik);
#endif
   SX_CHECK (stored == InMemory || stored == KeepOnDisk || stored == ReadOnDemand);
   if ( stored == KeepOnDisk || stored == ReadOnDemand)  {
      if (ik != loadedK || iSpin != loadedSpin)  {
         if (loadedK >= 0)  {
            SX_CHECK  (iSpin >= 0, iSpin);
            flushWaves ();
         }
         loadWaves (ik,iSpin);
      }
      return waves(0)(0);
   }  else  {
      waves(ik)(iSpin).auxData.iSpin = char(iSpin);
      waves(ik)(iSpin).auxData.ik    = ik;
      return waves(ik)(iSpin);
   }
}

PsiRef SxPW::operator() (const SxIdx &idx, int iSpin, int ik)
{
#ifdef SXPW_FINDLOOP
   findLoopUpdate (iSpin, ik);
#endif
   SX_CHECK (ik >= 0 && ik < getNk (), ik, getNk ());
   SX_CHECK (iSpin >= 0 && iSpin < getNSpin (), iSpin, getNSpin ());
   SX_CHECK (idx.start >= 0, idx.start);
   SX_CHECK (idx.end <= getNStates(ik), idx.end);
   PsiG &allStates = operator() (iSpin, ik);
   PsiRef res;
   int ng = (int)allStates.getNRows ();
   res = allStates(SxIdx (idx.start * ng, idx.end * ng + ng - 1));
   res.reshape (ng, idx.end - idx.start + 1);
   res.auxData = allStates.auxData;
   return res;
}

const PsiRef SxPW::operator() (const SxIdx &idx, int iSpin, int ik) const
{
#ifdef SXPW_FINDLOOP
   findLoopUpdate (iSpin, ik);
#endif
   SX_CHECK (ik >= 0 && ik < getNk (), ik, getNk ());
   SX_CHECK (iSpin >= 0 && iSpin < getNSpin (), iSpin, getNSpin ());
   SX_CHECK (idx.start >= 0, idx.start);
   SX_CHECK (idx.end <= getNStates(ik), idx.end);
   const PsiG &allStates = operator() (iSpin, ik);
   PsiRef res;
   int ng = (int)allStates.getNRows ();
   res = allStates(SxIdx (idx.start * ng, idx.end * ng + ng - 1));
   res.reshape (ng, idx.end - idx.start + 1);
   res.auxData = allStates.auxData;
   return res;
}

PsiRef SxPW::operator() (int i, int iSpin, int ik)
{
#ifdef SXPW_FINDLOOP
   findLoopUpdate (iSpin, ik);
#endif
   SX_CHECK (stored == InMemory || stored == KeepOnDisk || stored == ReadOnDemand || stored == ReadOneByOne);
   PsiRef psiOut;
   if ( stored == KeepOnDisk || stored == ReadOnDemand)  {
      if (ik != loadedK || iSpin != loadedSpin)  {
         flushWaves ();
         loadWaves (ik,iSpin);
      }

      psiOut = waves(0)(0).colRef(i);
   } else if (stored == ReadOneByOne) {
      if (   waves(0)(0).getSize () == 0
          || loadedI != i
          || waves(0)(0).auxData.ik != ik)
         loadWaves(i, iSpin, ik);
      psiOut = waves(0)(0);
   }  else  {
      psiOut = waves(ik)(iSpin).colRef(i);
   }
   psiOut.auxData.iSpin = char(iSpin);
   psiOut.auxData.ik    = ik;
   return psiOut;
}


const PsiRef SxPW::operator() (int i, int iSpin, int ik) const
{
#ifdef SXPW_FINDLOOP
   findLoopUpdate (iSpin, ik);
#endif
   SX_CHECK (stored == InMemory || stored == KeepOnDisk
             || stored == ReadOnDemand || stored == ReadOneByOne);
   PsiRef psiOut;
   if ( stored == KeepOnDisk || stored == ReadOnDemand)  {
      if (ik != loadedK || iSpin != loadedSpin)  {
         flushWaves ();
         loadWaves (ik,iSpin);
      }
      psiOut = waves(0)(0).colRef(i);
   } else if (stored == ReadOneByOne) {
      if (   waves(0)(0).getSize () == 0
          || loadedI != i
          || waves(0)(0).auxData.ik != ik)
         loadWaves(i, iSpin, ik);
      psiOut = waves(0)(0);
   }  else if (stored == InMemory)  {
      psiOut = waves(ik)(iSpin).colRef(i);
   }
   psiOut.auxData.iSpin = char(iSpin);
   psiOut.auxData.ik    = ik;
   return psiOut;
}


int SxPW::getNStates (int ik) const
{
   if (ik == -1) return nStates;
   SX_CHECK (ik >= 0 && ik < getNk (), ik, getNk ());
   if (stored == InMemory && waves(ik).getSize () > 0
       && waves(ik)(0).getSize () > 0)
   {
      SX_CHECK (nStatesPerK(ik) == waves(ik)(0).getNCols (),
                nStatesPerK(ik), waves(ik)(0).getNCols ());
      return (int)waves(ik)(0).getNCols ();
   }
   return nStatesPerK(ik);
}


void SxPW::changeNStates (int newN)
{
   SX_CHECK (newN > getNStates() || stored == ReadOnDemand ||
             stored == ReadOneByOne, newN, getNStates());

   int i, nOrig, iSpin, ng, ik, nk = getNk ();
   nkPoints = nk;
   nSpin = getNSpin();

   if ( stored == KeepOnDisk )  {
      // --- resize temporary file
      SX_EXIT;
      SxString tmpFilename = wavesFile.filename;
      wavesFile.close ();
      unlink (tmpFilename.ascii());  // remove old file
      createScratchFile ();
   }
   if (stored == ReadOnDemand || stored == ReadOneByOne)  {
      for (ik = 0; ik < nk; ++ik)  {
         if (nStatesPerK(ik) < newN)  {
            cout << "Too few states in '" << wavesFile.filename;
            cout << "' for SxPW::changeNStates to " << newN << endl;
            SX_EXIT;
         }
         nGIPerK(ik) /= nStatesPerK(ik);
         nStatesPerK(ik) = newN;
         nGIPerK(ik) *= newN;
      }
      loadedK = -1;
      loadedSpin = -1;
      loadedI = -1;
      return;
   }

   SxPW origWaves ( *this );

   for (ik=0; ik < nk; ik++)  {
      //waves.resize (nk);
      nOrig = origWaves.getNStates (ik);
      //waves(ik).resize (nSpin);
      nStatesPerK(ik) = newN;
      ng          = (int)origWaves.waves(ik)(0).getNRows ();
      nGIPerK(ik) = (int)origWaves.waves(ik)(0).getNRows () * newN;
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         SX_CHECK (nComp > 0, nComp);
         SX_CHECK ((char)nComp == nComp, nComp); // must be within char range
         waves(ik)(iSpin).resize  (ng * newN * nComp);
         waves(ik)(iSpin).reshape (ng * nComp, newN);
         waves(ik)(iSpin).auxData.nComp = (char)nComp;
         for (i=0; i < nOrig; i++)
            waves(ik)(iSpin).colRef(i) <<= origWaves(i,iSpin,ik);
         for (i=nOrig; i < newN; i++)
            waves(ik)(iSpin).colRef(i).randomize();
      }
   }
   orthonormalize ();
}



void SxPW::randomize ()
{
   int ik, nk = getNk(), iSpin;
   nkPoints = nk;
   SX_MPI_LEVEL("waves-k");
   for (ik=0; ik < nk; ik++)  {
      if (SxLoopMPI::myWork(ik))  {
         for (iSpin=0; iSpin < nSpin; iSpin++)
            (*this)(iSpin,ik).randomize();
      }
   }

//   orthonormalize ();
}



SxPW &SxPW::normalize ()
{
   int iSpin, ik, nk = int(waves.getSize());
   nkPoints = nk;
   for (ik=0; ik < nk; ik++)
      for (iSpin=0; iSpin < nSpin; iSpin++)
         waves(ik)(iSpin).normalize ();
   return *this;
}


SxPW &SxPW::orthogonalize (enum SxOrthoMethod)
{
   SX_EXIT;
   return *this;
}


SxPW &SxPW::orthogonalize (SxPW &)
{
   SX_EXIT;
   return *this;
}

void SxPW::setOrthogonal (SxVecRef<PrecCoeffG> *psiPtr,
                          int firstN, int iSpin, int ik,
                          enum NormMethod normal,
                          double threshold)
{
   // by C. Freysoldt
   // This routine is heavily commented on technical details
   // This is because I'm implementing it now, but the XPress library
   // will need changes here, but maybe not by me
   SX_CLOCK (Timer::ortho);

   SX_CHECK (iSpin >= 0 && iSpin < getNSpin (), iSpin, getNSpin ());
   SX_CHECK (ik >= 0 && ik < getNk (), ik, getNk ());
   SX_CHECK (firstN >= 0 && firstN <= getNStates (ik),
             firstN, getNStates (ik));
   SX_CHECK (threshold >= 0., threshold);
   SX_CHECK (   normal == DONT_NORMALIZE
             || normal == NORMALIZE);

   SX_CHECK (psiPtr);
   if (firstN == 0)  {
      if (normal == NORMALIZE) psiPtr->normalize ();
      return;
   }
   PsiRef &psi = *psiPtr;
   int nPsi = (int)psi.getNCols ();
   if (nPsi == 0) nPsi = 1;

   int j;
   PrecCoeffG scp;
   if (firstN > 150 || nPsi > 1)  {
      // for low n, the direct (convential) BLAS1 approach is faster
      // because there is a certain overhead, and BLAS2 routines may
      // not be fastet for BLAS1-like tasks
      // the cross-over point should be tested
      // this should be the zdotc / gemv crossover for nPsi == 1
      // this should be the gemv / gemm crossover for nPsi > 1
      SX_CHECK (stored == InMemory || loadedK == ik, loadedK, ik);
      SX_CHECK (stored == InMemory || loadedK >= 0, loadedK);
      SX_CHECK (stored == InMemory || loadedSpin >= 0, loadedSpin);

      PsiG *allStates
         = (stored == InMemory) ? &(waves(ik)(iSpin))
                                : &(waves(0)(0));

      // --- get matrix |j> for j < firstN
      int ng = (int)allStates->getNRows ();
      PsiRef lowerStates
         = (*allStates)(SxIdx(0, ng * firstN - 1));
      lowerStates.reshape(ng, firstN);
      // --- get <j|psi> for all 0 <= j < firstN

      SxVector<PrecCoeffG> scps = lowerStates.overlap (psi);

      // --- subtract |j><j|psi>
      // bool loopSubtr = (firstN < 500); // axpy / gemv crossover, to be tested
      bool loopSubtr = false;
      /*
      if (! loopSubtr && nPsi == 1)  {
         // count negligible scps
         int scpNegligible = firstN / 2; // if fewer negligible scps, use BLAS2
         SxVector<PrecCoeffG>::Iterator it = scps.begin ();
         for (j = 0; (j < firstN) && (scpNegligible > 0); j++, ++it)
            if ((*it).absSqr () > threshold) scpNegligible--;
         loopSubtr = (scpNegligible > 0);
      }
      */
      if (nPsi == 1 && loopSubtr)  {
         // BLAS1 version
         SxVector<PrecCoeffG>::Iterator it = scps.begin ();
         for (j = 0; j < firstN; j++, ++it)
            if ((*it).absSqr () > threshold)
               psi.plus_assign_ax (-(*it), lowerStates.colRef(j));
      } else {
         // BLAS2 (nPsi == 1) or BLAS3 version
         psi -= (lowerStates ^ scps);
      }

   } else {
      // --- conventional BLAS1 approach
      double scp2;
      for (j=0; j < firstN; j++)  {
         PsiRef psiJ = operator() (j,iSpin,ik);
         scp = dot (psiJ, psi); // <j|psi>
         scp2 = scp.absSqr ();
//       psi -= psiJ * scp;
         if (scp2 > threshold)
            psi.plus_assign_ax (-scp, psiJ); // psi -= |j><j|psi>
      }
   }

   if (normal == NORMALIZE)  {
      if (nPsi == 1)
         psi.normalize ();
      else
         for (int i = 0; i < nPsi; ++i) psi.colRef(i).normalize ();
   }
}

SxPW &SxPW::orthonormalize (enum SxOrthoMethod method,
                            SxArray<SxArray<SxVector<PrecCoeffG> > > *uPtr)
{
   SxPW &psiSet = *this;
//   timer->start (TIME_ORTHO);

   if (method == GramSchmidt)  {
      SxPWOverlap S;
      SX_MPI_LEVEL("waves-k");
      for (int ik = 0; ik < getNk (); ik++)  {
         if (!SxLoopMPI::myWork(ik)) continue;
         for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
            PsiG &psi = psiSet (iSpin, ik);
            S.orthonormalize (&psi, method);
         }
      }
   }  else  {
      // --- diagonalization, according to LOEWDIN orthogonalization scheme
      SX_CLOCK (Timer::ortho);

      SxVector<PrecCoeffG> S, U;

      SX_MPI_LEVEL("waves-k");
      for (int ik=0; ik < getNk (); ik++)  {
         if (!SxLoopMPI::myWork(ik)) continue;
         for (int iSpin=0; iSpin < nSpin; iSpin++)  {
            PsiG &psi = psiSet(iSpin, ik);

            // --- get overlap matrix S
            {
               SX_CLOCK (Timer::loewdinS);
               S = psi.overlap (psi);
            }

            // --- diagonalize S
            SxSymEigensystem<PrecCoeffG> eig (S);

            // --- set up U = S^{-1/2} from eigenbasis
            SxDiagMat<double> Ieps(1. / sqrt(eig.vals));
            U = eig.vecs ^ Ieps ^ eig.vecs.adjoint();

            // rotate wavefunctions as psi' = psi ^ U;
            psi.rotate (U);

            // --- store U, if required
            if (uPtr)  (*uPtr)(ik)(iSpin) = std::move (U);
         }
      }
   }

   return *this;
}


void SxPW::setGkBasisPtr (SxPtr<SxGkBasis> gkBasisPtrIn)
{
   gkBasisPtr = gkBasisPtrIn;
   SxGkBasis &gkBasis = *gkBasisPtr;
   int nk = getNk ();
   // Gk basis must fit unless this is uninitialized
   SX_CHECK (nk <= 0 || gkBasisPtr->getNk () == nk, nk, gkBasisPtr->getNk ());

   // set basis in waves
   if (stored == InMemory)  {
      for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
         for (int ik = 0; ik < nk; ik++)  {
           waves(ik)(iSpin).setBasis( &gkBasis(ik) );
         }
      }
   }   else if (loadedK >= 0)  {
      waves(0)(0).setBasis( &gkBasis(loadedK) );
   }
}


void SxPW::setNStates (const SxArray<int> &nPerK)
{
   SX_CHECK (waves.getSize() == nPerK.getSize(),
             waves.getSize(),   nPerK.getSize());
   SX_CHECK      (gkBasisPtr);

   int i, n, iSpin, ik, nGk, sizeOld, nOld;
   int nk = getNk ();
   const SxGkBasis &Gk = *gkBasisPtr;

   nStates = nPerK(0);
   if (stored == InMemory)  {

      for (ik=0; ik < nk; ik++)  {
         n = nPerK(ik);
         if (n < nStates) nStates = n; // get min. no. of states
         nGk = Gk(ik).ng;
         if ( n > nGk)  n = nGk;
//         SX_CHECK (n == nGk, n, nGk);
         nGIPerK(ik) = (int)waves(ik)(0).getNRows () * nPerK(ik);
         for (iSpin=0;  iSpin < nSpin;  iSpin++)  {
            sizeOld = (int)waves(ik)(iSpin).getSize ();
            nOld = (int)waves(ik)(iSpin).getNCols ();
            SX_CHECK (waves(ik)(iSpin).auxData.nComp == nComp,
                      (int)waves(ik)(iSpin).auxData.nComp, nComp);
            waves(ik)(iSpin).reshape (sizeOld, 1);
            waves(ik)(iSpin).resize  (nGk*n*nComp, true);
            waves(ik)(iSpin).reshape (nGk*nComp, n);

            // --- fill rest of wavefunctions with random numbers
            if (nOld < n)  {
               for (i = nOld; i < n; i++)  {
                  waves(ik)(iSpin).colRef(i).randomize ();
//                  waves(ik)(iSpin).normalize ();  // maybe, that gets necess.
               }
            }
         }
      }
   }  else  {
      // TODO: setNStates(nPerK) for case 'keepOnDisk' not yet implemented
      SX_EXIT;
   }
   UPDATE_MEMORY (waves);
   TRACK_MALLOC (*this, 1);
}


void SxPW::createScratchFile ()
{
   SX_CHECK (stored == KeepOnDisk);
   SX_CHECK (!wavesFile.isOpen);
   SX_CHECK (gkBasisPtr);
   const SxGkBasis &gkBasis = *gkBasisPtr;
   try  {
      wavesFile.open (SxBinIO::createTempName(tmpDir), SxBinIO::SCRATCH_READ_WRITE);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   waves.resize (1);  // only 1 k-point is kept in memory
   waves(0).resize (1); // resize by 1 spin
   size_t nElem, byteLen = sizeof (waves(0)(0)(0));
   int ik, iSpin;
   int nk = nkPoints, ng;
   nStatesPerK.resize (0);
   nGIPerK.resize (0);
   SX_CHECK (nComp > 0, nComp);
   SX_CHECK ((char)nComp == nComp, nComp); // must be within char range
   for (ik=0; ik < nk; ik++)  {
      ng = gkBasis(ik).ng;
      nElem = nStates*ng;
      nStatesPerK.append (nStates);
      nGIPerK.append ((int)nElem);
      // --- fill file buffer
      waves(0)(0).reformat(ng * nComp, nStates);
      waves(0)(0).auxData.nComp = (char)nComp;
      for (iSpin=0; iSpin < nSpin; iSpin++)  {
         fwrite (waves(0)(0).elements, byteLen, nElem * nComp, wavesFile.fp);
      }
   }
   loadedK    = -1;
   loadedSpin = -1;
}

void SxPW::flushWaves () const
{
   if (stored == ReadOnDemand || stored == ReadOneByOne) return;
   SX_CHECK (stored == KeepOnDisk);
   if (loadedK == -1) return;
   SX_CHECK (loadedK <= nkPoints, loadedK, nkPoints);
   int i;
   size_t pos, byteLen = sizeof (waves(0)(0)(0));
#  ifndef NDEBUG
   cout << "Flushing k-point " << (loadedK+1) << endl;
#  endif /* NDEBUG */
   for (pos=0, i=0; i < loadedK; i++)  pos += nGIPerK (i)*nSpin * nComp;
   pos += nGIPerK (loadedK)*loadedSpin * nComp;
   fseek (wavesFile.fp, pos*byteLen, SEEK_SET);
   fwrite (waves(0)(0).elements, byteLen, nGIPerK(loadedK) * nComp,
              wavesFile.fp);
}


void SxPW::loadWaves (int ik, int iSpin) const
{
   SX_CHECK (stored == KeepOnDisk || stored == ReadOnDemand);
   int i, ng;
   size_t pos=0, byteLen = sizeof (PrecCoeffG);
   SxVector<int> nPerK;
   /* FFT mapping not in use: assume G-basis fits, because it can be read
      from file.  This allows to use ReadOnDemand with modified FFT meshes.

   SxVector<int> fftIdx;

   */
#  ifndef NDEBUG
   cout << "Loading k-point " << (ik+1) << endl;
#  endif /* NDEBUG */

   bool readKwise = false;

   // --- get starting point
   if (stored == KeepOnDisk)  {
      for (i=0; i < ik; i++) pos += nGIPerK (i)*nSpin * nComp;
      pos += nGIPerK (ik)*iSpin * nComp;
      fseek (wavesFile.fp, pos*byteLen, SEEK_SET);
      ng = nGIPerK(ik) / nStatesPerK(ik);
   } else {
      if (nComp <= 0)  {
         if (wavesFile.containsDim("nComp"))
            const_cast<SxPW*>(this)->nComp = wavesFile.getDimension ("nComp");
         else
            const_cast<SxPW*>(this)->nComp = 1;
      }
      SX_CHECK ((char)nComp == nComp, nComp); // must be within char range
      nPerK.resize (nkPoints);
      wavesFile.read ("nPerK", &nPerK, nkPoints);
      SxVector<int> nGk(nkPoints);
      wavesFile.read ("nGk", &nGk, nkPoints);
      int fftOffset = 0;
      if (wavesFile.containsDim("nCoeff"))  {
         readKwise = false;
         for (i=0; i < ik; i++)  {
            pos += nPerK(i) * nGk(i) * nSpin * nComp;
            fftOffset += nGk(i);
         }
      } else {
         readKwise = true;
         pos = 0;
      }
      ng = nGk(ik);
      pos += nPerK(ik) * nGk(ik) * iSpin * nComp;
      /* FFT mapping not in use, see above
      if (gkBasisPtr)  {
         fftIdx.resize (ng);
         wavesFile.read ("fftIdx", &fftIdx, ng, fftOffset);
      }
      */
   }

   PsiG psiTmp;
   size_t offset = pos; // for ReadOnDemand
   if ((stored == ReadOnDemand) && (!gkBasisPtr))
      psiTmp.resize (ng * nComp);

   SxString psiVarName = readKwise ? SxString("psi-") + SxString(ik+1)
                                   : SxString("psi");

   waves(0)(0) = PsiG ();
   waves(0)(0).reformat(ng * nComp, nStatesPerK(ik));
   if (stored == KeepOnDisk)  {
      size_t nRead = fread (waves(0)(0).elements, byteLen, nGIPerK(ik) * nComp,
                            wavesFile.fp);
      if (nRead != byteLen * nGIPerK(ik) * nComp) { SX_EXIT; }
   } else  {
        for (i = 0; i < nStatesPerK(ik); ++i)  {
            PsiRef psi = waves(0)(0).colRef (i);
            /* FFT mapping not in use, see above

            if (gkBasisPtr)  {
               wavesFile.read (psiVarName, &psiTmp, ng, offset);
               psi <<= (*gkBasisPtr)(ik).mapToFFT (psiTmp, fftIdx);
            } else
            */
            {
               wavesFile.readVec (psiVarName, &psi, ng * nComp, (int)offset);
            }
            offset += ng * nComp;
         }
         offset += ng * nComp * (nPerK(ik) - nStatesPerK(ik));
      }
   if (gkBasisPtr)
      waves(0)(0).setBasis (& (*gkBasisPtr)(ik));
   waves(0)(0).auxData.iSpin = char(iSpin);
   waves(0)(0).auxData.ik    = ik;
   waves(0)(0).auxData.nComp = (char)nComp;
   loadedK = ik;
   loadedSpin = iSpin;
}

void SxPW::loadWaves (int i, int iSpin, int ik) const
{
   SX_CHECK (stored == ReadOneByOne);
   waves(0)(0) = PsiG ();
   waves(0)(0) = readPsi(wavesFile, i, iSpin, ik);
   if (gkBasisPtr)
      waves(0)(0).setBasis (& (*gkBasisPtr)(ik));
   loadedI = i;
}

PsiG SxPW::readPsi (const SxBinIO &io, int i, int iSpin, int ik)
{
   // SX_CHECK (nSpin == 1,nSpin);

   int ng;
   size_t pos=0;
   /* FFT mapping not in use: assume G-basis fits, because it can be read
      from file.  This allows to use ReadOnDemand with modified FFT meshes.

   SxVector<int> fftIdx;

   */

   // --- get starting point
   int nk, nSpin, nComp;
   try {
      nk    = io.getDimension ("nk");
      nSpin = io.getDimension ("nSpin");
      if (io.containsDim("nComp"))
         nComp = io.getDimension ("nComp");
      else
         nComp = 1;
      SX_CHECK ((char)nComp == nComp, nComp); // must be within char range
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   SxVector<int> nPerK(nk);
   SxVector<int> nGk(nk);
   try {
      io.read ("nPerK", &nPerK, nk);
      io.read ("nGk", &nGk, nk);
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   size_t fftOffset = 0;
   SxString psiVarName;
   if (io.containsDim("nCoeff"))  {
      psiVarName = "psi";
      for (i=0; i < ik; i++)  {
         pos += nPerK(i) * nGk(i);
         fftOffset += nGk(i);
      }
   } else {
      psiVarName = "psi-" + SxString(ik+1);
      pos = 0;
   }
   ng = nGk(ik);
   /* FFT mapping not in use, see above
   if (gkBasisPtr)  {
      fftIdx.resize (ng);
      io.read ("fftIdx", &fftIdx, ng, fftOffset);
   }
   */

   size_t offset = nSpin * nComp * pos + size_t(i) * ng * nComp;

   PsiG psi(ng * nComp);

   try {
      /* FFT mapping not in use, see above
      if (gkBasisPtr)  {
         PsiG psiTmp(ng);
         io.read (psiVarName, &psiTmp, ng, offset);
         psi <<= (*gkBasisPtr)(ik).mapToFFT (psiTmp, fftIdx);
      } else
      */
      io.readVec (psiVarName, &psi, ng, (int)offset);
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }
   SxAuxData &auxData = psi.auxData;
   auxData.iSpin = char(iSpin);
   auxData.ik    = ik;
   auxData.nComp = (char)nComp;
   return psi;
}

void SxPW::memMinimize ()
{
   if (stored == InMemory) return; // nothing to do
   // write to disk (if necessary)
   if (stored == KeepOnDisk) flushWaves ();
   // unload waves
   loadedK = -1;
   loadedSpin = -1;
   loadedI = -1;
   waves(0)(0) = PsiG ();
}

void SxPW::read (const SxBinIO &io, int mode)
{
   SX_CHECK (stored != ReadOnDemand && stored != ReadOneByOne);
   if (stored == Unknown) stored = InMemory;
   SX_CHECK (stored == InMemory || mode & KeepNStates); // to be implemented

   bool hasBasis = (gkBasisPtr.getPtr() != NULL);
   bool needsReorder = hasBasis;

   cout << "SxPW::read: |G+k> basis ";
   if (!hasBasis) cout << "not ";
   cout << "present." << endl;

   // --- variable for k-wise read (allows larger files)
   bool readKwise;
   SxString psiVarName;

   try  {
      int ik, nk, iSpin, nSpinFile, ng, i, iOffset;
#ifndef NDEBUG
      int nCoeff = -1;
#endif
      if(io.containsDim("nComp"))  {
         int nCompFile = io.getDimension("nComp");
         if (nComp > 0 && nComp != nCompFile)  {
            cout << "Inconsistency while reading " << io.filename << endl;
            cout << "Expected nComp=" << nComp << ", but found "
                 << nCompFile << " components in file" << endl;
            SX_QUIT;
         }
         nComp = nCompFile;
         SX_CHECK ((char)nComp == nComp, nComp); // must be within char range
      } else {
         nComp = 1;
      }
      nk        = io.getDimension ("nk");
      nSpinFile = io.getDimension ("nSpin");
      if (nSpinFile != nSpin && nSpin > 0)  {
         // can happen in non-DEBUG cases, so SX_CHECK not sufficient
         cout << "Inconsistency: this->nSpin = " << this->nSpin;
         cout << "; nSpin = " << nSpinFile << " in " << io.filename << endl;
         SX_EXIT;
      } else {
         nSpin = nSpinFile;
      }
      readKwise = !io.containsDim("nCoeff");
      if (readKwise)  {
         cout << "Reading waves k-wise" << endl;
      } else {
#ifndef NDEBUG
         nCoeff   = io.getDimension ("nCoeff");
#endif
         psiVarName = "psi";
      }
      sxprintf ("Reading wavefunction:\n");
      sxprintf ("---------------------\n");
      sxprintf ("nk: %d, nSpin: %d\n", nk, nSpin);

      SxString name;
      int offset = 0;
      int nAllGk = io.getDimension ("nAllGk");
      SxVector<int> nPerK (nk), nGk (nk), fftIdxIn(nAllGk);

      io.read ("nPerK",    &nPerK,    nk);
      io.read ("nGk",      &nGk,      nk);
      io.read ("fftIdx",   &fftIdxIn, nAllGk);
      // --- check if #k, k and #G(k) are correct
      if (hasBasis && nk != gkBasisPtr->nk)  {
         sxprintf ("Error: Input file is corrupt.\n");
         sxprintf ("       Number of k-points: %d (should be %d)\n", nk,
                 gkBasisPtr->nk);
         SX_EXIT;
      }
      if (hasBasis)  {
         SxVector<PrecG> kVecs(nk,3); // need matrix for reading
         Coord kVec;
         io.read ("kVec", &kVecs, nk, 3);
         for (ik = 0; ik < nk; ik++)  {
            kVec = kVecs.rowRef(ik).toVector3 ();
            if ( (kVec - gkBasisPtr->getK(ik)).normSqr () > 1e-10)  {
               cout << endl;
               cout << "Consistency check failed for k-point " << (ik+1);
               cout << ".\n Expected " << gkBasisPtr->getK(ik) << " but found ";
               cout << kVec<< " in '" << io.filename << "'." << endl;
               SX_EXIT; // may indicate changes in the k-point setup!
               // if this occurs, we should remap the k-points
            }
         }
      }
      for (ik=0; ik < nk; ik++)  {
         if (hasBasis && nGk(ik) != (*gkBasisPtr)(ik).ng)  {
            sxprintf ("Error: Input file is corrupt.\n");
            sxprintf ("       Number of |G+k> vectors: %d (should be %d)\n",
                     nGk(ik), (*gkBasisPtr)(ik).ng);
            SX_EXIT;
         }
      }

      SxMesh3D mesh;
      io.read("meshDim", &mesh);

      // --- resize to kpoints
      if (stored == InMemory)  {
         waves.resize (nk);
      }
      if (!(mode & KeepNStates))  nStatesPerK.resize (nk);
      SX_CHECK (nStatesPerK.getSize () == nk, nStatesPerK.getSize (), nk);

      nGIPerK.resize (nk);

      iOffset = 0;
      for (ik=0; ik < nk; ik++)  {
         if (readKwise)  {
            offset = 0;
            psiVarName = "psi-" + SxString(ik+1);
         }
         if (stored == InMemory)
            waves(ik).resize (nSpin);
         if (mode & KeepNStates)
            nStates = nStatesPerK(ik);
         else
            nStatesPerK(ik) = nStates = nPerK(ik);
         ng = nGk(ik);
         int n = nStatesPerK(ik);
         nGIPerK(ik) = n * ng;
         ssize_t ngc = ng * nComp;
         SxVecRef<PrecFFTIdx> fftIdx = fftIdxIn(SxIdx(iOffset, iOffset+int(ngc)-1));
         iOffset += int(ngc);
         if (hasBasis)  {
            if (mesh == (*gkBasisPtr)(ik).fft3d(0).mesh)  {
               // --- check n123 vs fftIdx
               SxVector<PrecFFTIdx>::Iterator n123It, fftIdxIt;
               n123It   = (*gkBasisPtr)(ik).n123(0).begin ();
               fftIdxIt = fftIdx.begin ();
               int ig;
               for (ig = 0; ig < ng; ++ig)
                  if (*fftIdxIt++ != *n123It++) break;
               needsReorder = (ig < ng);
            } else {
               needsReorder = true;
               SxVector3<int> vec;
               const SxFFT3d &fft = (*gkBasisPtr)(ik).fft3d(0);
               // --- check fftIdx fits into mesh
               for (int ig = 0; ig < ngc; ++ig)  {
                  vec = mesh.getMeshVec (fftIdx(ig), SxMesh3D::Origin);
                  if (   2 * abs(vec(0)) > fft.mesh(0)
                      || 2 * abs(vec(1)) > fft.mesh(1)
                      || 2 * abs(vec(2)) > fft.mesh(2))  {
                     cout << "Mesh incompatibility between file and basis."
                          << endl;
                     SX_EXIT;
                  }
               }
            }
         }

         for (iSpin=0; iSpin < nSpin; iSpin++)  {

            sxprintf ("reading state (%d,%d)...%dx%d elements\n",
                  iSpin, ik, int(ngc), n < nPerK(ik) ? n : nPerK(ik));
            fflush (stdout);

            // resize waves
            if (stored == InMemory)  {
               waves(ik)(iSpin).reformat (ngc, n);
               if (gkBasisPtr.getPtr() != NULL)
                  waves(ik)(iSpin).setBasis (&(*gkBasisPtr)(ik));  // UGLY
            } else  {
               // stored == KeepOnDisk
               waves(0)(0).reformat (ngc, n);
               loadedK = ik;
               loadedSpin = iSpin;
            }
            for (i=0; i < nPerK(ik); i++)  {
               if (i >= n)  {
                  // read fewer states than there are in file
                  // jump other states
                  offset += ng * (nPerK(ik) - n);
                  break;
               }

               PsiRef psi = (*this)(i,iSpin,ik);
               if (hasBasis && needsReorder)  {
                  PsiG psiTmp(ngc);
                  io.readVec (psiVarName, &psiTmp, int(ngc), offset);
                  SX_VALIDATE_VECTOR (psiTmp);
                  psi <<= (*gkBasisPtr)(ik).mapToFFT (psiTmp, fftIdx, mesh);
               } else {
                  io.readVec (psiVarName, &psi, int(ngc), offset);
                  SX_VALIDATE_VECTOR (psi);
               }
               offset += int(ngc);
            }

            // randomize rest of waves
            for (i = nPerK(ik); i < n; i++)  {
               (*this)(i, iSpin, ik).randomize ();
               (*this)(i, iSpin, ik).normalize ();
            }
            if (stored == KeepOnDisk)
               flushWaves ();

         } // iSpin

         // at the end, all coefficients must be read
         SX_CHECK ( (!readKwise) || offset == ngc * nPerK(ik) * nSpin,
                   offset, ngc * nPerK(ik) * nSpin);
      }
      // at the end, all coefficients must be read
      SX_CHECK (readKwise || offset == nCoeff, offset, nCoeff);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   UPDATE_MEMORY (waves);
   TRACK_MALLOC (*this, 1);

   //--- if GkBasis should be read from file create it from file and set
   if(!(mode & KeepGkBasis)) {
      //TODO CHECK MACRO for debugging only, overwrite GK basis is a bad idea but who knows when it is needed.
      SX_CHECK(gkBasisPtr.getPtr() == NULL);
      if (gkBasisPtr.getPtr() != NULL)  {
         cout << SX_SEPARATOR;
         cout << "WARNING: Overwrite GkBasisPtr in SxPW!" << endl;
         cout << SX_SEPARATOR;
      }
      if (mode & SaveMemory) gkBasisPtr = SxPtr<SxGkBasis>::create(io, true, true);
      else gkBasisPtr = SxPtr<SxGkBasis>::create(io, true, false);
      setGkBasisPtr(gkBasisPtr);
   }
}

void SxPW::write (SxBinIO &io) const
{
   int nk = getNk ();
   // for large numbers of kPoints reaching ncdump limit
   // set writekwise kPoint dependent
   bool writeKwise = nk < 1000;
   if (stored == KeepOnDisk && io.ncMode == SxBinIO::WRITE_DATA)
      flushWaves ();
   try  {
      int iSpin;
      int ik;
      int nCoeff = 0, i, n, ng;

      SxVector<int> nPerK(nk), nGk(nk);
      for (ik=0; ik < nk; ik++)  {
         nPerK(ik) = n  = getNStates (ik);
         nGk(ik)   = ng = nGIPerK(ik) / n;
         nCoeff += n * ng * nSpin * nComp;
      }
      int offset = 0;

      // --- create dimensions
      io.addDimension ("nk",        nk);
      if (writeKwise) {
         for (ik = 0; ik < nk; ++ik)
            io.addDimension ("nCoeff-"+SxString(ik+1), nGIPerK(ik) * nSpin * nComp);
      } else {
         io.addDimension ("nCoeff",    nCoeff);
      }
      io.addDimension ("nSpin",     nSpin);
      io.addDimension ("nComp",     nComp);

      // --- write data
      if (!io.contains ("nGk") || (io.ncMode == SxBinIO::WRITE_DATA))
         io.write ("nGk",   nGk,   "nk");
      if (!io.contains ("nPerK") || (io.ncMode == SxBinIO::WRITE_DATA))
         io.write ("nPerK", nPerK, "nk");
      SxString psiVarName = "psi",
               coeffDimName = "nCoeff";
      SX_NEW_LOOP (*this);
      for (ik=0; ik < nk; ik++)  {
         if (writeKwise)  {
            offset = 0;
            SxString sK(ik+1);
            psiVarName = "psi-" + sK;
            coeffDimName = "nCoeff-" + sK;
         }
         for (iSpin=0; iSpin < nSpin; iSpin++)  {

            if (   stored != InMemory && stored != Unknown
                && io.ncMode == SxBinIO::WRITE_DATA)
            {
               loadWaves (ik,iSpin);
            }
            SX_CHECK (   io.ncMode != SxBinIO::WRITE_DATA
                      || (*this)(iSpin,ik).getSize () == nGk(ik) * nPerK(ik) * nComp,
                      (*this)(iSpin,ik).getSize (), nGk(ik) * nPerK(ik) * nComp);

            for (i=0; i < nPerK (ik); i++)  {
               io.writeVec (psiVarName,
                            /* real psi only upon real writing
                               to avoid superfluous io when waves are
                               in scratch file on disk */
                            (io.ncMode == SxBinIO::WRITE_DATA)
                               ? (*this)(i, iSpin, ik)
                               : SxVecRef<PrecCoeffG> (),
                            coeffDimName, offset);
               offset += nGk(ik) * nComp;
            }
         }
         // at the end, all coefficients must be written
         SX_CHECK ((!writeKwise) || offset == nGIPerK(ik) * nSpin * nComp,
                   offset, nGIPerK(ik) * nSpin * nComp);
      }

      // at the end, all coefficients must be written
      SX_CHECK (writeKwise || offset == nCoeff, offset, nCoeff);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}


size_t SxPW::getNBytes () const
{
   return ::getNBytes (waves)
        + ::getNBytes (nStates)
        + ::getNBytes (nGIPerK);
}

void SxPW::registerMemoryObservers ()
{
   TRACK_MEMORY (waves);
   TRACK_MEMORY (nStates);
   TRACK_MEMORY (nGIPerK);
   TRACK_MALLOC (*this, 1);
}

void SxPW::setZero()
{
   SX_CHECK(stored == InMemory);
   for (int ik = 0; ik < getNk(); ik++)
      for (int iSpin = 0; iSpin < getNSpin(); iSpin++)
            waves(ik)(iSpin).set(0.);

}

#ifdef SXPW_FINDLOOP
void SxPW::findLoopUpdate (int iSpin, int ik) const
{
   if (iSpin == findLoopLastIk && ik == findLoopLastIk) return;
   if (iSpin == 0 && ik == 0)  {
      cout << "NEW waves loop!" << endl;
   }
   findLoopLastIk = ik; findLoopLastIspin = iSpin;
}

void SxPW::findLoopReset () const
{
   findLoopLastIk = findLoopLastIspin = 0;
}
#endif

