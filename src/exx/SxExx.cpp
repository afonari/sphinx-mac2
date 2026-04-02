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

#include <SxExx.h>
#include <SxConfig.h>
#include <SxChi0.h>
#include <SxRhoIndG.h>
#include <SxFockGk.h>
#include <SxSymFockGk.h>

SxExx::SxExx ()
   : SxXC(),
     nGChi(-1)
{
   // empty
}

SxExx::SxExx (const SxXC &xc)
   : SxXC(xc),
     nGChi(-1)
{
   vExUnsym.resize (nSpin); vSlaterUnsym.resize (nSpin);  // TODO: strike out
   vExLoopG.resize (nSpin);
}

SxExx::~SxExx ()
{
   // empty
}

void SxExx::read (const SxSymbolTable *table, const SxGBasis &G)
{
   try  {
      const SxSymbolTable *top = table->topLevel();
      const SxSymbolTable *basis = top->getGroup("basis");

      // --- get number of G vectors
      if (basis->contains("eCutChi"))  {
         PrecEnergy eCutChi = basis->get("eCutChi")->toReal();
         int ig = 0;
         while (G.g2(ig) < eCutChi && ig < G.ng-1)  ig++;
         nGChi = ig;
         cout << "eCutChi = " << eCutChi << " --> nGChi = " << nGChi << endl;
      }

   }  catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
} 

void SxExx::updateXC (const RhoR &rhoR, const SxPWSet *wavesPtr,
                      const SxFermi *fermiPtr)
{
   const SxRBasis *rBasisPtr 
      = dynamic_cast<const SxRBasis *>(rhoR(0).getBasisPtr());
   SX_CHECK (rBasisPtr);

   int iSpin;

   // --- add core charge density in case of NLCC
   for (iSpin=0; iSpin < nSpin; iSpin++)  {
      if (nlcc)
         rhoXcR(iSpin) = rhoR(iSpin) + rhoCoreR;
      else
         rhoXcR(iSpin) = rhoR(iSpin);
   }

   // --- compute total charge density: sum up spin and core
   SxMeshR rhoTlR;
   if (nSpin == 1)
      rhoTlR = rhoXcR(0);
   else  {
      rhoTlR = SxMeshR (0., *rBasisPtr);
      for (iSpin=0; iSpin < nSpin; iSpin++)  rhoTlR += rhoXcR(iSpin);
   }

   // --- XC functional
   switch (xcFunctional)  {
      case SxXC::EXX : SX_CHECK (wavesPtr);
                        SX_CHECK (fermiPtr);
                        if (nSpin == 1)  {
                           computeEXX_noCorr (*wavesPtr, *fermiPtr);
                        }  else  {
                           sxprintf ("Spin polarization for EXX not yet "
                                     "implemented. Sorry.\n");
                           SX_QUIT;
                        }
                        vCor(0) = SxMeshR (0., *rBasisPtr);
                        eCorrelation = eC = 0.;
                        break;
      case SxXC::EXX_LDA :
                        SX_CHECK (wavesPtr);
                        SX_CHECK (fermiPtr);
                        if (nSpin == 1)  {
                           computeEXX_LDA (rhoTlR, *wavesPtr, *fermiPtr);
                        }  else  {
                           sxprintf ("Spin polarization for EXX not yet "
                                     "implemented. Sorry.\n");
                           SX_QUIT;
                        }
                        break;
      case SxXC::EXX_PBE :
                        SX_CHECK (wavesPtr);
                        SX_CHECK (fermiPtr);
                        if (nSpin == 1)  {
                           computeEXX_PBE (rhoTlR, *wavesPtr, *fermiPtr);
                        }  else  {
                           sxprintf ("Spin polarization for EXX not yet "
                                     "implemented. Sorry.\n");
                           SX_QUIT;
                        }
                        break;
      case SxXC::EXX_PBE_WB :
                        SX_CHECK (wavesPtr);
                        SX_CHECK (fermiPtr);
                        if (nSpin == 1)  {
                           computeEXX_PBE_WB (rhoTlR, *wavesPtr, *fermiPtr);
                        }  else  {
                           sxprintf ("Spin polarization for EXX not yet "
                                     "implemented. Sorry.\n");
                           SX_QUIT;
                        }
                        break;
      case SxXC::Slater :
                        SX_CHECK (wavesPtr);
                        SX_CHECK (fermiPtr);
                        vEx(0)  = SxMeshR (*rBasisPtr);
                        vCor(0) = SxMeshR (0., *rBasisPtr);
                        eCorrelation = eC = 0.;
                        if (nSpin == 1)  {
                           computeSlater (rhoTlR, *wavesPtr, *fermiPtr);
                        }  else  {
                           sxprintf ("Spin polarisation for Slater potential "
                                     "not yet implemented. Sorry.\n");
                           SX_QUIT;
                        }
                        break;
      case SxXC::KLI  : SX_CHECK (wavesPtr);
                        SX_CHECK (fermiPtr);
                        vEx(0)  = SxMeshR (*rBasisPtr);
                        vCor(0) = SxMeshR (0., *rBasisPtr);
                        eCorrelation = eC = 0.;
                        if (nSpin == 1)  {
                           computeKLI_noCorr (rhoTlR, *wavesPtr, *fermiPtr);
                        }  else  {
                           sxprintf ("Spin polarisation for KLI not yet "
                                     "implemented. Sorry.\n");
                           SX_QUIT;
                        }
                        break;
      default         : sxprintf ("Unknown XC functional.\n");
                        SX_EXIT;
   }

   double eVEx, eVCor;

   // --- compute exchange-correlation potential energy
   for (iSpin=0, eVxc=0.; iSpin < nSpin; iSpin++)  {

      if (!calcX)  {
         vXc(iSpin) -= vEx(iSpin);
         vEx(iSpin).set (0.);
         eExchange   = 0.;
      }
      if (!calcC)  {
         vXc(iSpin)  -= vCor(iSpin);
         vCor(iSpin).set (0.);
         eCorrelation = 0.;
      }
      
      eVxc += tr(rhoR(iSpin) * vXc(iSpin));
      if (xcFunctional != SxXC::READ_VXC)  {
         eVEx = tr(rhoR(iSpin) * vEx(iSpin));
         eVCor = tr(rhoR(iSpin) * vCor(iSpin));
         sxprintf("spin: %d  eVEx: %.15f \n", iSpin, eVEx); 
         sxprintf("spin: %d  eVCor: %.15f \n", iSpin, eVCor);
      }
   }

}


void SxExx::computeEXX_noCorr (const SxPWSet &waves, const SxFermi &fermi)
{
   cout << "SxExx::computeEXX_noCorr -- Hallo!\n" << endl;
   SX_CHECK (fockPtr);
   const SxRBasis *rBasisPtr = dynamic_cast<const SxRBasis *>
                              (rhoXcR(0).getBasisPtr ());
   SX_CHECK (rBasisPtr);

   const SxRBasis &R = *rBasisPtr;
   const SxGBasis &G = R.getGBasis ();
//sxprintf ("outside of ::computeOld -> %p\n", pwBasis);
//cout << "fermi -> " << fermi.eps(0,0) << endl;
//cout << "fermi.eps.getSize() : " << fermi.eps(0,0).getSize() << endl;

   const SxPW *wavesPtr = dynamic_cast<const SxPW *> (&waves);
   SX_CHECK (wavesPtr);

   // --- initialise chi_0 and the Fock-induced charge density
   SxChi0    chi0 (G, waves.getGkBasis(), nGChi);
   SxRhoIndG rhoInd (G, waves.getGkBasis(), nGChi);

   if (useCorrectChi)  chi0.compute (*wavesPtr, fermi);
   else                chi0.computeOld (*wavesPtr, fermi);  // Vogl version
//cout << "chi0.diag --> " << chi0.getDiag() * 100. << endl;
//cout << "chi0 (neu) = " << chi0.getChi0 () << endl;
//chi0.computeOld (*wavesPtr, fermi);
//cout << "chi0 (alt) = " << chi0.getChi0 () << endl;
//SX_EXIT;

   rhoInd.compute (*wavesPtr, fermi, *fockPtr);
//cout << "E --> " << rhoInd.getRhoIndG () * 100. << endl;
//SX_EXIT;

   SxMeshG vExG(G.ng, 0.);
   vExG.setBasis (&G);
   vExG( SxIdx(1,nGChi-1) ) = rhoInd.getRhoIndG () ^ chi0.getInverse ();
//cout << "vExG :  " << vExG * 100. << endl;
//SX_EXIT;

SX_FFT_REAL fftFactorGtoR = R.fft3d.scaleFor;  // = 1/sqrt(Omega)
//cout << "fftFactorGtoR = " << fftFactorGtoR << endl;
//vExG(0) = 1.;

   vExUnsym(0) = SxVector<PrecCoeffG> (vExG);

compareRhoInd (rhoInd.getRhoIndG(), vExG, *wavesPtr, fermi);
//SX_EXIT;

   vEx(0) = (R | vExG).real() 
                 * fftFactorGtoR;  // Why? Compare with the probably wrong
                                   // fft factor added at the end of
                                   // rhoIndG::compute!
//cout << "vEx(0) :  " << vEx(0) << endl;
//SX_EXIT;
   vEx(0) = R.symmetrize (vEx(0));

//vExG = (G | vEx(0)) / fftFactorGtoR;
//cout << "vExG :  " << vExG * 100. << endl;
//SX_EXIT;

   eXc = eExchange = eX = rhoInd.getXEnergy ();

   vXc(0) <<= vEx(0);
}

void SxExx::computeEXX_LDA (const SxMeshR &rho,
                           const SxPWSet &waves, const SxFermi &fermi)
{
   bool calcXo = calcX;
   cout << "SxExx::computeEXX_LDA -- Hallo!\n" << endl;
   calcX = false;
   rhoXcR(0) = rho;
   computeLDA (); // correlation within LDA
   calcX = calcXo;
   computeEXX_noCorr (waves, fermi);  // overwrites x and xc of LDA computation
                                      //    with x = xc = EXX
   vXc(0) += vCor(0);   // xc = x(EXX) + c(LDA)
   eXc    += eC;

   shiftXC ();
}

void SxExx::computeEXX_PBE (/*const*/ SxMeshR &/*rho*/,
                           const SxPWSet &waves, const SxFermi &fermi)
{
   cout << "SxExx::computeEXX_PBE -- Hallo!\n" << endl;
//cout << "rho = " << rho << endl;
//sxprintf ("sum = %20.15f\n", rho.sum());
   SX_EXIT; // vCor/vEx not set up in updateGGA_PBE_WB
   //updateGGA_PBE_WB (rho);            // x, c, xc within PBE (WB)
   computeEXX_noCorr (waves, fermi);  // overwrites x and xc of PBE computation
                                      //    with x = xc = EXX
   vXc(0) += vCor(0);   // xc = x(EXX) + c(PBE)
//cout << "vCor = " << vCor(0) << endl;
//sxprintf ("sum = %20.15f\n", vCor(0).sum());
//SX_EXIT;
   eXc    += eCorrelation; // that's correct
//   eXc    += eC;                    // that's wrong, since eC is local energy
}

void SxExx::computeEXX_PBE_WB (/*const*/ SxMeshR &/*rho*/,
                              const SxPWSet &waves, const SxFermi &fermi)
{
   cout << "SxExx::computeEXX_PBE_WB -- Hallo!\n" << endl;
   SX_EXIT; // vCor/vEx not set up in updateGGA_PBE_WB
   //updateGGA_PBE_WB (rho);            // x, c, xc within PBE (WB)
   computeEXX_noCorr (waves, fermi);  // overwrites x and xc of PBE (WB)
                                      //    computation with x = xc = EXX
   vXc(0) += vCor(0);                 // xc = x(EXX) + c(PBE(WB))
   eXc    += eCorrelation;
}

void SxExx::computeKLI_noCorr (const SxMeshR &rho,
                              const SxPWSet &waves, const SxFermi &fermi)
{
   cout << "SxExx::computeKLI_noCorr -- Hallo!\n" << endl;

   int               ik, nk = waves.getNk();
   int               iv, nv, vIdx;
   SxMeshR           vkRho, q;
   PrecWeights       kWeight;
   const SxGkBasis  &gkSet = waves.getGkBasis ();  // set of |G+k> basises
   const SxRBasis *rBasisPtr 
      = dynamic_cast<const SxRBasis *>((rho.getBasisPtr ()));
   SX_CHECK (rBasisPtr);
   const SxRBasis   &R = *rBasisPtr;
   PrecEnergy        vkPotOld;

   // - - - - - - - - - - - - - - - - - - - - - - - - - -
   nv = fermi.getNValenceBands(0,0);

//TODO: test norm of psi!!

   // --- store old x-potential
   VxcR vXcOld(nSpin);
   for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
      vXcOld(iSpin).resize (R.getMeshSize());
      vXcOld(iSpin) <<= vXc(iSpin);
   }

   // --- get the Slater potential
   VxcR vSlater (getVSlater(rho, waves, fermi, true));

   // --- calculate the correction term
   q.resize (R.getMeshSize());
   q.setBasis (&R);
   q.set (0.);
   for (ik = 0; ik < nk; ik++)  {
      kWeight = gkSet.weights(ik);

      for (iv = 0; iv < nv; iv++)  {
         vIdx   = fermi.getValenceBandIdx (iv,0,ik);
         PsiR vkPsiR = R | waves (vIdx,0,ik);
         vkRho  = vkPsiR.conj() * vkPsiR;

         // --- get expectation value of old vXc
         vkPotOld = tr(vkRho * vXcOld(0));

         // --- get contribution of vk
         q += kWeight * vkRho * (vkPotOld - slaterExpect(ik)(iv));
      }
   }

   // --- division by rho and multiplication by 2 for spin
   q = 2. * q / rho;

   // --- symmetrisation
   q = R.symmetrize (q);

   // --- add correction term to Slater potential
   vEx(0) = vSlater(0) + q;
   vXc(0) <<= vEx(0);

   // --- set (G=0)-component to 0 for comparison
   shiftXC ();
}

void SxExx::computeSlater (const SxMeshR &rho,
                          const SxPWSet &waves, const SxFermi &fermi)
{
   cout << "SxExx::computeSlater -- Hallo!\n" << endl;

   VxcR vSlater (getVSlater(rho, waves, fermi));
   vEx(0) <<= vSlater(0);
   vXc(0) <<= vSlater(0);

   // --- set (G=0)-component to 0 for comparison
   shiftXC ();
}

VxcR SxExx::getVSlater (const SxMeshR &rho,
                       const SxPWSet &waves, const SxFermi &fermi, bool kli)
{
   cout << "SxExx::slater -- Hallo!\n" << endl;

//   int               iSpin;
   int               ik, nk = waves.getNk();
   int               iv, nv, vIdx;
   VxcR              vSlater(nSpin);
   PrecWeights       kWeight;
   SxFockGk         &FockOp = *fockPtr;
   const SxGkBasis  &gkSet = waves.getGkBasis ();  // set of |G+k> basises
   const SxRBasis *rBasisPtr 
      = dynamic_cast<const SxRBasis *>((rho.getBasisPtr ()));
   SX_CHECK (rBasisPtr);
   const SxRBasis   &R = *rBasisPtr;
//const SxGBasis   &G = *pwBasis;
   const SxPW       *wavesPtr = dynamic_cast<const SxPW *> (&waves);
   SX_CHECK            (wavesPtr);
//SX_FFT_REAL       fftFactorGtoR = R.fft3d.scaleFor;
//cout << "FFT factor: " << fftFactorGtoR << endl;

   // - - - - - - - - - - - - - - - - - - - - - - - - - -
   nv = fermi.getNValenceBands(0,0);

   // TODO: The following lines should go into SxFermi.
   // --- check the equality of nv for all k-points
   for (ik = 0; ik < nk; ik++)  {
      if (fermi.getNValenceBands(0,ik) != nv)  {
         sxprintf ("Number of valence states is not the same for all ");
         sxprintf ("k-points. Sorry.\n");
         SX_QUIT;
      }
   }

   eX = 0.;
   if (kli)  {
      slaterExpect.resize (nk);
      for (ik = 0; ik < nk; ik++)  slaterExpect(ik).resize (nv);
   }

   // --- compute Slater potential
   PsiR p(0., R);
   for (ik = 0; ik < nk; ik++)  {
      kWeight = gkSet.weights(ik);
      FockOp.compute (ik, *wavesPtr, fermi);

      for (iv = 0; iv < nv; iv++)  {
         vIdx   = fermi.getValenceBandIdx (iv,0,ik);
         const PsiRef &vPsiG = waves (vIdx,0,ik);
         PsiR vPsiR  = R | vPsiG;
         PsiG vFockG = FockOp * vPsiG;
         PsiR vFockR = R | vFockG;

         // --- counter of Slater potential
         PsiR pvk  = vPsiR.conj() * vFockR;
         p.plus_assign_ax (kWeight, pvk);

         // --- expectation value of pvk
         if (kli)  slaterExpect(ik)(iv) = tr(pvk);

         // --- exchange contribution to the total energy
         eX += kWeight * dot(vFockG, vPsiG).re;
      }
   }

   // --- exchange contribution to total energy
   eXc = eExchange = eX;

   // --- division by rho and multiplication by 2 for spin
   vSlater(0) = 2. * p / rho;

// --- for exx debugging only
const SxGBasis &G = R.getGBasis ();
vSlaterUnsym(0) = ( G | vSlater(0) );

   // --- symmetrisation
   vSlater(0) = R.symmetrize (vSlater(0));

   return vSlater;
}

//-----------------------------------------------------------------------------
//   gsEXX  <-->  stEXX
//-----------------------------------------------------------------------------
PsiG
SxExx::greensfunction (const PsiRef &h1psi0, double epsV,
                       const SxPW &waves, const SxFermi &fermi)
{
   const SxGkBasis &Gk = waves.getGkBasis();

   SX_CHECK (h1psi0.getBasisPtr() == &(Gk(h1psi0.auxData.ik)));

   int    iSpin  = h1psi0.auxData.iSpin;
   int    ik     = h1psi0.auxData.ik;
   int    nc     = fermi.getNConductionBands(0,0);
   PsiG psi1(0., Gk(ik));
   SxComplex16  elem;

   cout << "nc = " << nc << endl;

   // --- initialise |psi1>
   psi1.auxData.iSpin = char(iSpin);
   psi1.auxData.ik    = ik;

   // --- Greens function
   for (int ic = 0; ic < nc; ic++)  {

      int cIdx = fermi.getConductionBandIdx (ic,iSpin,ik);
      const PsiRef &psiC = waves(cIdx,iSpin,ik);
      double epsC = fermi.eps(cIdx,iSpin,ik);

      // < psi0 | H1 | psi1>
      elem = dot(psiC, h1psi0);

      // G H1 | psi1 >
      psi1 += psiC * elem / (epsV - epsC);
   }

   return psi1;
}

PsiG SxExx::probier (const SxVecRef<PrecCoeffG> &vExG, const SxPW &waves, const SxFermi &fermi)
{
   SX_CHECK (dynamic_cast<const SxGBasis *>(vExG.getBasisPtr ()));
   const SxGBasis  &G  = *dynamic_cast<const SxGBasis *>(vExG.getBasisPtr ());
   const SxRBasis  &R  = G.getRBasis ();
   const SxGkBasis &Gk = waves.getGkBasis();

   int  nValStates = fermi.getNValenceBands(0,0);
   int  nc         = fermi.getNConductionBands(0,0);
   cout << "nValBands = " << nValStates << endl;
   cout << "nc        = " << nc << endl;

   int  iSpin = 0;
   int  ik = 0;
   int  mu, muIdx;

   PsiR vExR = ( R | vExG );
   PsiR resR(0., R);

   for (mu = 0; mu < nValStates; mu++)  {
      muIdx = fermi.getValenceBandIdx(mu,iSpin,ik);
      const PsiRef &psi0G = waves(muIdx,iSpin,ik);
      PsiR psi0R = ( R | psi0G );

      PsiG h1psi0G = ( Gk(ik) | (vExR * psi0R ) );
      // TODO: ugly
      double epsV = fermi.eps(muIdx,iSpin,ik);
      h1psi0G.auxData.iSpin = (char)iSpin;
      h1psi0G.auxData.ik    = ik;

      PsiR psi1GreensR
         = R | greensfunction (h1psi0G, epsV, waves, fermi);

      resR += psi0R.conj() * psi1GreensR;
   }

   return 4. * R.cell.volume * ( G | resR );
}

void SxExx::compareRhoInd (const SxVecRef<PrecCoeffG> &rhoIndG, const SxVecRef<PrecCoeffG> &vExG,
                          const SxPW &waves, const SxFermi &fermi)
{
   cout << "Compare rho_ind ..." << endl;

   SX_CHECK (dynamic_cast<const SxGBasis *>(vExG.getBasisPtr ()));
   const SxGBasis  &G  = *dynamic_cast<const SxGBasis *>(vExG.getBasisPtr ());
   const SxRBasis  &R  = G.getRBasis ();
   const SxGkBasis &Gk = waves.getGkBasis();

   int   ik = 0;
   int   iSpin = 0;
   int   i = 0;
   const PsiRef &psi0G = waves(i,iSpin,ik);
   PsiR  psi0R = ( R | psi0G );
   int   ng = G.ng;

   PsiR  vExR = ( R | vExG );
   PsiG  h1psi0 = ( Gk(ik) | (vExR * psi0R) );
   // TODO: ugly
   h1psi0.auxData.iSpin = char(iSpin);
   h1psi0.auxData.ik    = ik;

//   PsiG psi1Greens = greensfunction (h1psi0, waves, fermi);
//   cout << "psi1Greens = " << psi1Greens << endl;

   cout << SX_SEPARATOR;

   PsiG res = probier (vExG, waves, fermi);
   cout << "probier = " << res << endl;

   cout << SX_SEPARATOR;

   PsiG rhoInd(ng);
   rhoInd(0) = 0.;
   rhoInd(SxIdx(1,nGChi-1)) = rhoIndG;
   rhoInd(SxIdx(nGChi,ng-1)).set (0.);
   cout << "rhoInd = " << rhoInd << endl;

   cout << SX_SEPARATOR;

   PsiG diff = res - rhoInd;
   cout << "diff = " << diff << endl;
   cout << "|diff| = " << diff.norm () << endl;

   cout << SX_SEPARATOR;

   PsiRef bier = res( SxIdx(1,nGChi-1) );
   PsiG diffNGChi = bier - rhoIndG;
   cout << "diff (nGChi) = " << diffNGChi << endl;
   cout << "|diff (nGChi)| = " << diffNGChi.norm () << endl;
}

void SxExx::computeVExUnsym (const SxPWSet &waves, const SxFermi &fermi,
                             const SxGBasis &G)
{
   cout << "SxExx::computeVExUnsym -- Hallo!\n" << endl;
   SX_CHECK (fockPtr);

   const SxPW *wavesPtr = dynamic_cast<const SxPW *> (&waves);
   SX_CHECK (wavesPtr);

   // --- initialise chi_0 and the Fock-induced charge density
   SxChi0    chi0 (G, waves.getGkBasis(), nGChi);
   SxRhoIndG rhoInd (G, waves.getGkBasis(), nGChi);

   chi0.computeOld (*wavesPtr, fermi);
   rhoInd.compute (*wavesPtr, fermi, *fockPtr);

   SxMeshG vExG(G.ng, 0.);
   vExG.setBasis (&G);
   vExG( SxIdx(1,nGChi-1) ) = rhoInd.getRhoIndG () ^ chi0.getInverse ();

cout << "rhoInd (stEXX) = " << rhoInd.getRhoIndG () << endl;

   vExUnsym(0) = PsiG (vExG);
}

void SxExx::computeFockEnergy (const SxPW &waves, const SxFermi &fermi)
{
   SX_CHECK (fock.getSize() > 0, fock.getSize());

   int               iv, vIdx, iSpin, ik;
   int               nValStates = fermi.getNValenceBands(0,0);
   int               nk         = waves.getNk();
   SX_CHECK (nSpin == waves.getNSpin(), nSpin, waves.getNSpin());
   const SxGkBasis  &Gk         = waves.getGkBasis();
   PrecWeights       kWeight;

   SX_CHECK (nSpin == 1, nSpin);  // spin pol. not yet implemented

   eX = 0.;

   for (ik = 0; ik < nk; ik++)  {
      kWeight = Gk.weights(ik);

      for (iSpin = 0; iSpin < nSpin; iSpin++)  {
         for (iv = 0; iv < nValStates; iv++)  {
            vIdx     = fermi.getValenceBandIdx (iv,iSpin,ik);
            const PsiRef &psi = waves (vIdx,iSpin,ik);
            PsiG fockPsi  = *(fock(ik)) * psi;
            eX += kWeight * dot(fockPsi, psi).re;
         }
      }
   }

   eXc = eExchange = eX;
}

void SxExx::computeSymFockEnergy (const SxPW &waves, const SxFermi &fermi)
{
   SX_CHECK (symFock.getSize() > 0, symFock.getSize());

   int               iv, vIdx, iSpin, ik;
   int               nValStates = fermi.getNValenceBands(0,0);
   int               nk         = waves.getNk();
   SX_CHECK (nSpin == waves.getNSpin(), nSpin, waves.getNSpin());
   const SxGkBasis  &Gk         = waves.getGkBasis();
   PrecWeights       kWeight;

   SX_CHECK (nSpin == 1, nSpin);  // spin pol. not yet implemented

   eX = 0.;

   for (ik = 0; ik < nk; ik++)  {
      kWeight = Gk.weights(ik);

      for (iSpin = 0; iSpin < nSpin; iSpin++)  {
         for (iv = 0; iv < nValStates; iv++)  {
            vIdx     = fermi.getValenceBandIdx (iv,iSpin,ik);
            const PsiRef &psi = waves (vIdx,iSpin,ik);
            PsiG fockPsi  = *(symFock(ik)) * psi;
            eX += kWeight * dot(fockPsi, psi).re;
         }
      }
   }

   eXc = eExchange = eX;
}


