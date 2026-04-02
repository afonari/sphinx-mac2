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

#include <SxPartialWaveBasis.h>
#include <SxPAWPot.h>
#include <SxPAWBasis.h>

SxPartialWaveBasis::SxPartialWaveBasis (const SxConstPtr<SxPAWPot> &pot,
                                        const SxAtomicStructure &str)
   : nElements(0)
{
   SX_CHECK (pot);
   potPtr = pot;

   SX_CHECK (potPtr->getNSpecies () == str.getNSpecies (),
             potPtr->getNSpecies (), str.getNSpecies ());
               
   for (int is = 0; is < str.getNSpecies (); ++is)
      nElements += potPtr->getNProj (is) * str.getNAtoms (is);

   SX_CHECK (nElements > 0);
}

void SxPartialWaveBasis::print () const
{
   cout << SX_SEPARATOR;
   cout << "| Partial wave basis |p>" << endl;
   cout << "| Number of partials:      " << nElements << endl;
   cout << "| Projections to Gk-basis: " 
        << (projectors ? "enabled" : "disabled") << endl;
   cout << SX_SEPARATOR;
}

SxPtr<SxAOBasis> 
SxPartialWaveBasis::createProjBasis (const SxGkBasis &gk,
                                     const SxPAWPot &pawPot)
{
   SxArray<SxVector<double> > projRad(pawPot.getNSpecies ());
   for (int is = 0; is < pawPot.getNSpecies (); ++is)  {
      // enforce referencing in assignments below...
      SxVecRef<double> &projRef = projRad(is);
      if (    pawPot.pPsFine.getSize () > 0
           && pawPot.pPsFine(is).getSize () > 0)
      {
         projRef = pawPot.pPsFine(is);
      } else {
         projRef = pawPot.pPS(is);
      }
   }
   return SxPtr<SxAOBasis>::create (gk, projRad, pawPot.lPhi);
}

PsiRef
SxPartialWaveBasis::projectFrom (const SxPAWBasis *pawIn,
                                 const SxVecRef<SxComplex16> &psi) const
{
   SX_CHECK (pawIn->gBasis);
   SX_CHECK (this == pawIn->pBasis.getPtr ());
   int nPsi = (int)psi.getNCols ();
   SX_CHECK (psi.getSize () > 0);
   PsiG res (getNElements (), nPsi);
   res.auxData = psi.auxData;
   res.setBasis (this);
   // --- copy partial wave part into new vector
   res <<= psi.getRef<SubMatrix> (pawIn->gBasis->ng, getNElements (),
                                  psi.getNCols ());
   return res;
}

PsiRef SxPartialWaveBasis::projectFrom (const SxGBasis *,
                                        const PsiRef &psi) const
{
   SX_CHECK (projectors);
   SxVector<PrecCoeffG> res = projectors->fromPWBasis (psi);
   res.setBasis (this);
   res.auxData.ik = psi.auxData.ik;
   return res;
}


