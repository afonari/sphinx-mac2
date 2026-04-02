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

#include <SxSpeciesRef.h>

SxAtomicStructure SxSpeciesRef::operator* (const SxAtomicStructure &in) const
{
   SX_CHECK (in.atomInfo);
   int nAtom = in.getNAtoms(iSpec);

   SxAtomicStructure res;
   res.cell = in.cell;

   // --- create reference to coordinates
   res.set (const_cast<SxVecRef<PrecTauR>&>(in.coords) // cast away constness
              .getRef<Compact>(in.atomInfo->offset(iSpec) * 3, 3, nAtom, 1));

   // --- create new atomInfo
   SxAtomInfo::Ptr aInfo;
   if (mode == MapToFiltered)  {
      aInfo = SxAtomInfo::derive (in.atomInfo);
   }  else  {
      SX_CHECK (in.atomInfo->getParent ());
      aInfo = SxAtomInfo::derive (in.atomInfo->getParent ());
   }

   aInfo->nSpecies = in.atomInfo->nSpecies;
   aInfo->nAtoms.resize (aInfo->nSpecies);
   aInfo->nAtoms.set (0);
   aInfo->nAtoms(iSpec) = nAtom;
   aInfo->setupOffset ();

   // --- create map to parent
   if (mode == MapToFiltered)  {
      aInfo->parentMap.resize (nAtom);
      SxVector<int>::Iterator mapIt = aInfo->parentMap.begin ();
      int idx = in.atomInfo->offset(iSpec);
      for ( ; nAtom; --nAtom)
         *mapIt++ = idx++;
   } else {
      SX_CHECK (MapToParentOfFiltered);
      aInfo->parentMap = in.atomInfo->parentMap(in.getRange (iSpec));
      // TODO: this won't be possible in XPress. Remove this mode.
   }

   res.atomInfo = aInfo;
   
   return res;

}

