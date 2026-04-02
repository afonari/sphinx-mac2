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

#include <SxForstData.h>

SxForstData::SxForstData () : SxSpeciesData ()
{
   // empty
}


SxForstData::SxForstData (const SxSymbolTable *table) : SxSpeciesData (table)
{
   try  {
      // someValue = table->get("someIdentifier")->toDouble ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}


SxForstData::~SxForstData ()
{
   // empty
}



