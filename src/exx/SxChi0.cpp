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

#include <SxChi0.h>
#include <SxEigensystem.h>
#include <SxIdx.h>

SxChi0::SxChi0 ()
   : nG      (-1),
     GPtr    (NULL),
     gkPtr   (NULL),
     nk      (-1)
{
   // --- set uninitialized variables to nonsense value
   nv = nc = -1;
   fftFactorRtoG = 0.;
   sizeG = -1;
}

SxChi0::SxChi0 (const SxGBasis &G, const SxGkBasis &gk, int nG_)
   : nG      (nG_),
     GPtr    (&G),
     gkPtr   (&gk),
     nk      (gkPtr->getNk())
{
   // --- set uninitialized variables to nonsense values
   nv = nc = -1;

   // --- initialise R basis
   R = G.getRBasis ();

   // --- reformat Chi0
   Chi0.reformat (nG-1, nG-1);  // no (G = 0) components considered

   // --- fft prefactor
   //   ( To check the validity of this factor, use SxChi0::checkFFTfactor! )
   fftFactorRtoG = 1. / R.fft3d.scaleFor;

   // --- size of the compressed G mesh the vectors are stored in
   sizeG = G.ng;
}

SxChi0::SxChi0 (SxGBasis &G, SxGkBasis &gk, int nG_, const SxCell &cell)
   : nG      (nG_),
     GPtr    (&G),
     gkPtr   (&gk),
     R       (G.fft3d(0).mesh, cell),
     nk      (gkPtr->getNk())
{
   // --- set uninitialized variables to nonsense value
   nv = nc = -1;

   // --- reformat Chi0
   Chi0.reformat (nG-1, nG-1);  // no (G = 0) components considered

   // --- fft prefactor
   //   ( To check the validity of this factor, use SxChi0::checkFFTfactor! )
   fftFactorRtoG = 1. / G.fft3d(0).scaleFor;

   // --- size of the compressed G mesh the vectors are stored in
   sizeG = (int)G.n123.getSize();
}

SxChi0::SxChi0 (SxGBasis &G, SxGkBasis &gk, const SxCell &cell)
   : GPtr    (&G),
     gkPtr   (&gk),
     R       (G.fft3d(0).mesh, cell),
     nk      (gkPtr->getNk())
{
   // --- set uninitialized variables to nonsense values
   nG = nv = nc = -1;

   // --- fft prefactor
   //   ( To check the validity of this factor, use SxChi0::checkFFTfactor! )
   fftFactorRtoG = 1. / R.fft3d.scaleFor;

   // size of the compressed G mesh the vectors are stored in
   sizeG = (int)G.n123.getSize();
}

void SxChi0::checkFFTfactor (const SxPW &waves)
{
   SX_CHECK (GPtr);

   const SxGBasis    &G    = *GPtr;
   const PsiRef &psiG = waves (0,0,0);
   PrecCoeffG         g0;
   double             g0abs, eps = 1.e-10;

   PsiR psiR  = ( R | psiG );
   SxMeshG rhoVV = ( G | (psiR * psiR.conj()) );
   g0    = rhoVV(0) * fftFactorRtoG;
   g0abs = sqrt (g0.absSqr ());

   if (fabs(g0abs - 1.) < eps)  {
      cout << "SxChi0::testFFTfactor --> (G=0)-component: " << g0
           << " ... ok" << endl;
   }  else  {
      cout << "SxChi0::testFFTfactor --> (G=0)-component: " << g0
           << " ... failed" << endl;
      SX_QUIT;
   }
}
      
void SxChi0::compute (const SxPW &waves, const SxFermi &fermi)
{
   SX_CHECK  (waves.getNSpin() == 1, waves.getNSpin());
               // spin pol. not yet implemented
   SX_CHECK (waves.getNk() == nk, waves.getNk(), nk);
   SX_CHECK      (GPtr);

   int              ik;
   int              iG, iGp;
   int              iv, ic, vIdx, cIdx;
   PrecEps          vEps, cEps, vcEps;
   const SxGBasis  &G = *GPtr;
   PrecWeights      kWeight;

   // - - - - - - - - - - - - - - - - - - - - - - - - - -
   nv = fermi.getNValenceBands(0,0);
   nc = fermi.getNConductionBands(0,0);

   // TODO: The following lines should go to SxFermi.
   // --- check for equality of nv at all k-points
   for (ik = 0; ik < nk; ik++)  {
      if (fermi.getNValenceBands(0,ik) != nv)  {
         sxprintf ("Number of valence states is not the same for all k-points. ");
         sxprintf ("Sorry.\n");
         SX_QUIT;
      }
   }

   // --- find the smallest number of conduction bands
   for (ik = 0; ik < nk; ik++)  {
      if (fermi.getNConductionBands(0,ik) < nc)  {
         nc = fermi.getNConductionBands(0,ik);
      }
   }
   // - - - - - - - - - - - - - - - - - - - - - - - - - -

   sxprintf ("SxChi0::compute ... using %d G vectors "
           "and %d conduction bands\n", nG, nc);
   sxprintf ("SxChi0::compute ... not using k <=> -k\n");
   cout << "SxChi0::compute ... starting\n";  cout.flush();

   Chi0.set (0.);

   for (ik = 0; ik < nk; ik++)  {
      kWeight = gkPtr->weights(ik);

      for (iv = 0; iv < nv; iv++)  {
         vIdx  = fermi.getValenceBandIdx (iv,0,ik);
         vEps  = fermi.eps (vIdx,0,ik);
         PsiR vPsiR = R | waves (vIdx,0,ik);

         for (ic = 0; ic < nc; ic++)  {
            cIdx  = fermi.getConductionBandIdx (ic,0,ik);
            cEps  = fermi.eps (cIdx,0,ik);
            PsiR cPsiR = R | waves (cIdx,0,ik);
            // Abdullah version
            SxMeshG rhoVC = ( G | (vPsiR * cPsiR.conj()) );
            SxMeshG rhoCV = ( G | (cPsiR * vPsiR.conj()) );
//cout << SX_SEPARATOR;
//cout << "rhoVC = " << rhoVC << endl;
//cout << SX_SEPARATOR;
//cout << "rhoCV = " << rhoCV << endl;
//EXIT;
            vcEps = 2. * kWeight / (vEps - cEps);

            // TODO: this is a symmetric rank 1 update
            // Chi0(G,G') -> Chi0(G,G') + alpha x(G)x(G')*
            // with alpha = vcEps and x = (rhoVC * rhoCV).conj ()
            // => use a special call for this
            for (iG = 1; iG < nG; iG++)  {
               for (iGp = iG; iGp < nG; iGp++)  {
                  Chi0(iG-1,iGp-1) +=     vcEps
                                      * ( rhoVC(iG).conj() * rhoVC(iGp)
                                      +   rhoCV(iG).conj() * rhoCV(iGp) );
/*
cout << SX_SEPARATOR;
cout << "new (iG = " << iG << ", iGp = " << iGp << ") = "
     << rhoVC(iG).conj() * rhoVC(iGp) + rhoCV(iG).conj() * rhoCV(iGp) << endl;
cout << "old (iG = " << iG << ", iGp = " << iGp << ") = "
     << 2. * ( rhoVC(iG).conj() * rhoVC(iGp) ) << endl;
cout << "rhoVC(iG)  = " << rhoVC(iG) << endl;
cout << "rhoVC(iGp) = " << rhoVC(iGp) << endl;
cout << "rhoCV(iG)  = " << rhoCV(iG) << endl;
cout << "rhoCV(iGp) = " << rhoCV(iGp) << endl;
*/
               }  // :iGp
            }  // :iG
         }  // :ic
         sxprintf ("k=%d, v=%d ... done\n", ik, iv);  cout.flush();
      }  // :iv
   }  // :ik

   // --- fill lower triangle of matrix
   for (iG = 1; iG < nG; iG++)  {
      for (iGp = iG; iGp < nG; iGp++)  {
         Chi0(iGp-1,iG-1) = Chi0(iG-1,iGp-1).conj();
                            // conj() just to assure hermitecity even in case
                            // of numerical inaccuracy
      }  // :iGp
   }  // :iG

   // --- multiply with FFT factor
   // 
   //   ( Within the upper calculation the following scheme is applied:
   //     <G|psi_i> --> <R|psi_i> --> <G|psi_i*psi_j> =: rhoVC without
   //     considering any FFT prefactor. However, in the second trafo
   //     occurs a (pointwise) product of two wavefunctions. For one
   //     function, missing factors cancel due to the consecutive
   //     performance of backward and forward trafo. For the other
   //     wavefunction, however, one needs to care about the FFT factors. )
   Chi0 *= fftFactorRtoG * fftFactorRtoG;

   // --- divide by omega
//   Chi0 /= R.cell.volume;

// --- conjucate chi0, do NOT adjoint it
//Chi0 = Chi0.conj ();
      
   cout << "SxChi0::compute ... done\n";  cout.flush();
}

void SxChi0::computeOld (const SxPW &waves, const SxFermi &fermi)
{
   SX_CHECK  (waves.getNSpin() == 1, waves.getNSpin());
               // spin pol. not yet implemented
   SX_CHECK (waves.getNk() == nk, waves.getNk(), nk);
   SX_CHECK      (GPtr);

   int              ik;
   int              iG, iGp;
   int              iv, ic, vIdx, cIdx;
   PrecEps          vEps, cEps, epsVC;
   const SxGBasis  &G = *GPtr;
   PrecWeights      kWeight;
sxprintf ("within ::computeOld -> %p\n", (const void *)GPtr);

   // ---------------------------------------------------
const PsiRef &psiV_G = waves(0,0,0);
cout << "waves(0,0,0).sum () = " << psiV_G.sum() << endl;

   // - - - - - - - - - - - - - - - - - - - - - - - - - -
   nv = fermi.getNValenceBands(0,0);
   nc = fermi.getNConductionBands(0,0);

   // TODO: The following lines should go to SxFermi.
   // --- check for equality of nv at all k-points
   for (ik = 0; ik < nk; ik++)  {
      if (fermi.getNValenceBands(0,ik) != nv)  {
         sxprintf ("Number of valence states is not the same for all k-points. ");
         sxprintf ("Sorry.\n");
         SX_QUIT;
      }
   }

   // --- find the smallest number of conduction bands
   for (ik = 0; ik < nk; ik++)  {
      if (fermi.getNConductionBands(0,ik) < nc)  {
         nc = fermi.getNConductionBands(0,ik);
      }
   }
   // - - - - - - - - - - - - - - - - - - - - - - - - - -

   sxprintf ("SxChi0::computeOld ... using %d G vectors "
           "and %d conduction bands\n", nG, nc);
   cout << "SxChi0::computeOld ... starting\n";  cout.flush();

   Chi0.set (0.);

   for (ik = 0; ik < nk; ik++)  {
      kWeight = gkPtr->weights(ik);

      for (iv = 0; iv < nv; iv++)  {
         vIdx  = fermi.getValenceBandIdx (iv,0,ik);
         vEps  = fermi.eps (vIdx,0,ik);
         PsiR vPsiR = R | waves (vIdx,0,ik);
// - - - fac
//cout << "im G+k Raum: " << (vPsiG ^ vPsiG).chop() << endl;
//cout << "      L�ge: " << vPsiG.getSize() << endl;
//cout << "im r   Raum: " << (vPsiR ^ vPsiR).chop() << endl;
//cout << "      L�ge: " << vPsiR.getSize() << endl;
//
//SxMeshG vPsiGG = ( G | vPsiR );
//cout << "im G   Raum: " << (vPsiGG ^ vPsiGG).chop() << endl;
//cout << "      L�ge: " << vPsiGG.getSize() << endl;
//
//SxMeshG rhoVV = ( G | (vPsiR * vPsiR.conj()) );
//cout << "rhoVV(G=0) = " << rhoVV(0) << endl;
// - - - fac

         for (ic = 0; ic < nc; ic++)  {
            cIdx  = fermi.getConductionBandIdx (ic,0,ik);
            cEps  = fermi.eps (cIdx,0,ik);
            PsiR cPsiR = R | waves (cIdx,0,ik);

            // Abdullah version
            SxMeshG rhoVC = ( G | (vPsiR * cPsiR.conj()) );
//            rhoVC = ( G | (vPsiR.conj() * cPsiR) );  // --> footnote (1)
            epsVC = 4. * kWeight / (vEps - cEps);

            for (iG = 1; iG < nG; iG++)  {
               for (iGp = iG; iGp < nG; iGp++)  {
                  Chi0(iG-1,iGp-1) += epsVC * rhoVC(iG).conj() * rhoVC(iGp);
               }  // :iGp
            }  // :iG
         }  // :ic
         sxprintf ("k=%d, v=%d ... done\n", ik, iv);  cout.flush();
      }  // :iv
   }  // :ik

   // --- fill lower triangle of matrix
   for (iG = 1; iG < nG; iG++)  {
      for (iGp = iG; iGp < nG; iGp++)  {
         Chi0(iGp-1,iG-1) = Chi0(iG-1,iGp-1).conj();
                            // conj() just to assure hermitecity even in case
                            // of numerical inaccuracy
      }  // :iGp
   }  // :iG

   // --- multiply with FFT factor
   // 
   //   ( Within the upper calculation the following scheme is applied:
   //     <G|psi_i> --> <R|psi_i> --> <G|psi_ipsi_j> =: rhoVC without
   //     considering any FFT prefactor. However, in the second trafo
   //     occurs a (pointwise) product of two wavefunctions. For one
   //     function, missing factors cancel due to the consecutive
   //     performance of backward and forward trafo. For the other
   //     wavefunction, however, one needs to care about the FFT factors. )
   Chi0 *= fftFactorRtoG * fftFactorRtoG;

   // --- divide by omega
//   Chi0 /= R.cell.volume;  // test 4th July

// --- conjucate chi0, do NOT adjoint it
//Chi0 = Chi0.conj ();
      
   cout << "SxChi0::computeOld ... done\n";  cout.flush();
}

void SxChi0::computeEigVals ()
{
   eigVals = symEigenvalues (Chi0);
}

void SxChi0::computeEigSys ()
{
   SxSymEigensystem<PrecCoeffG> eigSys(Chi0);
   eigVals = std::move(eigSys.vals);
   eigVecs = std::move(eigSys.vecs);
}

PsiG SxChi0::getEigVecG (const int iG)
{
   SX_CHECK (iG >= 0 && iG < nG-1, iG, nG);
   SX_CHECK  (eigVecs.getSize() > 0, eigVecs.getSize());

   PsiG out(*GPtr);
   out(0) = 0.;
   out( SxIdx(1,nG-1) ) <<= eigVecs.colRef(iG);
   return out;
}

PsiR SxChi0::getEigVecR (const int iG)
{
   SX_CHECK (iG >= 0 && iG < nG-1, iG, nG);
   SX_CHECK  (eigVecs.getSize() > 0, eigVecs.getSize());

   PsiR eigVecR = ( R | getEigVecG(iG) );

   // The following corresponds to a division by the phasefactor
   // exp{ip} (with the phase p), so that the vector is completely
   // real, and a subsequent muliplication by cos p.
   for (int d = 0; d < eigVecR.getSize(); d++)  eigVecR(d).im = 0.;

   return eigVecR;
}

PsiG SxChi0::operator* (const SxVecRef<PrecCoeffG> &v)
{
   int nGin = (int)v.getSize();

   SX_CHECK (nGin > nG, nGin, nG);

   PsiG out(nGin);
   out(0) = 0.;
   out(SxIdx(1,nG-1)) = Chi0 ^ v.getRef(nG-1,1);
   out(SxIdx(nG,nGin-1)).set (0.);
   return out;
}

// --- I/O
void SxChi0::write (const SxString &filename) const
{
   SX_CHECK  (eigVecs.getSize() > 0, eigVecs.getSize());
   SX_CHECK (eigVals.getSize() == nG-1, eigVals.getSize(), nG);

   int dim = nG-1;

   try  {
      SxBinIO io;
      io.open (filename, SxBinIO::BINARY_WRITE_ONLY);

      // --- create dimensions
      io.addDimension ("dim", dim);  // chi_0 is a (dim x dim) matrix
      io.addDimension ("nv", nv);    // number of valence bands used
      io.addDimension ("nc", nc);    // number of conduction bands used

      // --- write data
      io.write ("chi0", Chi0, "dim", "dim");
      io.write ("eigVals", eigVals, "dim");
      io.write ("eigVecs", eigVecs, "dim", "dim");
      io.setMode (SxBinIO::WRITE_DATA);
      io.write ("chi0", Chi0, "dim", "dim");
      io.write ("eigVals", eigVals, "dim");
      io.write ("eigVecs", eigVecs, "dim", "dim");

      io.close ();
   }  catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

void SxChi0::read (const SxString &filename)
{
   try  {
      int dim;
      SxVector<double> eigValsDouble;

      SxBinIO io;
      io.open (filename, SxBinIO::BINARY_READ_ONLY);

      // --- get dimensions
      dim = io.getDimension ("dim");
      nG  = dim + 1;  // nG considers (G=0)-component, too
      nv  = io.getDimension ("nv");
      nc  = io.getDimension ("nc");

      // --- resize vectors/matrices
      Chi0.reformat (dim, dim);
      eigValsDouble.resize (dim);
      eigVecs.reformat (dim, dim);

      // --- get data
      io.read ("chi0", &Chi0, dim, dim);
      io.read ("eigVals", &eigValsDouble, dim);
      io.read ("eigVecs", &eigVecs, dim, dim);

      // --- cast EWs to complex vector
      eigVals = SxVector<PrecCoeffG> (eigValsDouble);

      io.close ();
   }  catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

// ==========================================================================
//    fotenotes
// ==========================================================================
// (1)
// Is the original version of Abdullah replaced by the line commented out
// ("to comment out" really exists), all eigenvalues remain the same. When
// representing the eigenvectors in real space, their real part also is in-
// variant against swapping those lines. Their strange imaginary part, however,
// is replaced by its negative. Shortly, eigenvectors in real space become
// complex conjugated.
// Representing the eigenvectors in reciprocal space, nothing can be intuitivly
// seen.
