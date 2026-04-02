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

#include <SxRBasis.h>
#include <SxGBasis.h>

SxRBasis::SxRBasis ()
   : dOmega (0.),
     fft3d (),
     gBasisPtr(NULL)
{
   // empty
}

SxRBasis::SxRBasis (int nx, int ny, int nz, const SxCell &cellIn) 
   : fft3d (SxFFT3d::Reverse, nx, ny, nz, cellIn.volume, true),
     cell(cellIn),
     gBasisPtr(NULL)
{
   fft3d.autoscale = false;
  
   dOmega     = cell.volume / (nx * ny * nz);
}

void SxRBasis::set (int nx, int ny, int nz, const SxCell &cellIn)
{
   fft3d.symmetric = true;
   fft3d.dir = SxFFT3d::Reverse;
   fft3d.setMesh (nx, ny, nz, cellIn.volume);
   fft3d.autoscale = false;
   cell = cellIn;
   dOmega     = cell.volume / (nx * ny * nz);
}

void SxRBasis::set (const SxVector3<int> &mesh, const SxCell &cellIn)
{
   set (mesh(0), mesh(1), mesh(2), cellIn);
}

SxRBasis::SxRBasis (const SxVector3<int> &mesh, 
                    const SxCell &cellIn)
   : fft3d (SxFFT3d::Reverse, mesh(0), mesh(1), mesh(2), cellIn.volume),
     cell(cellIn),
     gBasisPtr(NULL)
{
   fft3d.autoscale = false;
   dOmega     = cell.volume / fft3d.meshSize;
}
      

SxRBasis::~SxRBasis ()
{
   deregisterAll ();
}

void SxRBasis::deregister (const SxBasis *basis) const
{
   if (gBasisPtr == basis) gBasisPtr = NULL;
   deregisterBasis (basis);
}

Real8 SxRBasis::tr (const SxVecRef<double> &vec) const
{
   SX_CHECK (vec.getBasisPtr() == this);
   return vec.sum() * dOmega;
}

Real8 SxRBasis::tr (const SxVecRef<SxComplex16> &vec) const
{
   SX_CHECK (vec.getBasisPtr() == this);
   SxComplex16 sum = vec.sum();
   SX_CHECK (fabs(sum.im) < 1e-12 || fabs(sum.im) < 1e-12 * fabs(sum.re),
             sum.re, sum.im);
   return sum.re * dOmega;
}


SxVecRef<SxRBasis::CoeffType>
SxRBasis::identity (const SxRBasis *rBasis,
                    const SxVecRef<CoeffType> &vec) const
{
   // --- verify that this is really the identity projector, i.e.
   //     that the provided basis is the same as the vector's basis
   SX_CHECK (vec.getBasisPtr() == this);
   SX_CHECK (rBasis == this);
   return const_cast<SxVecRef<SxRBasis::CoeffType>&>(vec);
}


SxVector<SxRBasis::CoeffType>
SxRBasis::toGSpace (const SxGBasis *gBasis,
                    const SxVecRef<CoeffType> &psiR) const
{  SX_CLOCK (Timer::RtoG);
   SX_CHECK (gBasis);
   // "handshake"
   gBasis->registerBasis(*this); // does nothing if *this is already registered
   int nComp = psiR.auxData.nComp;
   SX_CHECK (psiR.getSize() == fft3d.meshSize * nComp,
             psiR.getSize(), fft3d.meshSize, nComp);

   // --- set up fft mesh
   int ng = gBasis->ng;
   int idx = gBasis->getBasisId (this);
   int size = (int)psiR.getSize() / nComp;
   SX_CHECK (idx >=0 && idx < gBasis->n123.getSize (),
             idx, gBasis->n123.getSize ());
   SX_CHECK (fft3d.mesh(0) == gBasis->fft3d(idx).mesh(0),
             fft3d.mesh(0), gBasis->fft3d(idx).mesh(0));
   SX_CHECK (fft3d.mesh(1) == gBasis->fft3d(idx).mesh(1),
             fft3d.mesh(1), gBasis->fft3d(idx).mesh(1));
   SX_CHECK (fft3d.mesh(2) == gBasis->fft3d(idx).mesh(2),
             fft3d.mesh(2), gBasis->fft3d(idx).mesh(2));
   SxVector<CoeffType> psiOut (ng*nComp);
   psiOut.auxData = psiR.auxData;
   psiOut.setBasis (gBasis);

   if (gBasis->hasMixedFFT (fft3d.mesh))  {
      SxFFT2d1d &mixedFFT = *(gBasis->getMixedFFT ());
      SxComplex16* mesh = mixedFFT.getMesh ();
      for (int c = 0; c < nComp; c++)  {
         // do the FT
         {  SX_CLOCK (Timer::RtoG_FFT);
            mixedFFT.fftBackward (psiR.elements + c * fft3d.meshSize);
         }
         SX_FFT_REAL scale = fft3d.scaleRevFFT;
#ifdef USE_OPENMP
#pragma omp parallel for if (ng > sxChunkSize)
#endif
         // get data from mesh
         for (ssize_t ig = 0; ig < ng; ++ig)
            psiOut(ig + c * ng) = scale * mesh[gBasis->n231(ig)];

      }
      mixedFFT.freeMesh ();
      SX_VALIDATE_VECTOR (psiOut);
      return psiOut;
   }

   for (int c = 0; c < nComp; c++)  {
      // --- set up single component
      const SxVecRef<CoeffType> &psiRC = psiR.getRef(size, c * size);

      const SxVecRef<SxComplex16> psiIn  (const_cast<SxVecRef<CoeffType>&>(psiRC));
      PrecFFTIdx *n123Ptr = gBasis->n123(idx).elements;
      SxVector<SxComplex16> psiWork (fft3d.meshSize);

      {  SX_CLOCK (Timer::RtoG_FFT);
         fft3d.fftReverse (fft3d.meshSize, psiIn.elements,psiWork.elements);
      }
      PrecCoeffG *destPtr = psiOut.elements + c * ng;
      //#  ifdef USE_OPENMP
      //      if (ng > chunkSize)  {
      //#        pragma omp parallel for 
      //         for (ig=0; ig < ng; ig++)
      //            destPtr[ig] = psiWork.elements[n123Ptr[ig]] * fft3d.scaleRevFFT;
      //      } else
      //#  endif /* USE_OPENMP */
      {
         for (ssize_t ig=0; ig < ng; ig++)
            *destPtr++ = psiWork.elements[*n123Ptr++] * fft3d.scaleRevFFT;
      }
   }

   // cout << "--> toGSpace called!  nComp = " << nComp << endl;

   SX_VALIDATE_VECTOR (psiOut);
   return psiOut;
}



void SxRBasis::writeMesh3d (const SxString &file, 
                            const TRho &mesh3d) const
{
   try  {
      SxBinIO io (file, SxBinIO::BINARY_WRITE_ONLY);
      io.writeMesh (mesh3d, cell, fft3d.mesh);
      io.setMode (SxBinIO::WRITE_DATA);
      io.writeMesh (mesh3d, cell, fft3d.mesh);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}


RhoR SxRBasis::readMesh3d (const SxString &file) const
{
   RhoR meshes;
   try  {
      SxBinIO io (file, SxBinIO::BINARY_READ_ONLY);
      meshes = io.readMesh ();
      io.close ();

      sxprintf ("SxRBasis::readMesh3d  -  consistency checks not implemented\n");

   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   return meshes;
}



void SxRBasis::print () const
{
   sxprintf ("SxRBasis:\n");
}
