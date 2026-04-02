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
#include <SxSimpleParser.h>


double SxParserGet::operator|| (double defaultValue) const
{
   double res = defaultValue;
   if (symbol) {
      try {
         res = symbol->toReal ();
      } catch (SxException e)  {
         theParser.error ();
         e.print ();
         SX_EXIT;
      }
   }
   return res;
}

/// get string from symbol or default value
SxString SxParserGet::operator|| (const char *defaultValue) const
{
   if (symbol) {
      SxString res;
      try {
         res = symbol->toString ();
      } catch (SxException e)  {
         theParser.error ();
         e.print ();
         SX_EXIT;
      }
      return res;
   }
   return defaultValue;
}

/// get string from symbol or default value
SxString SxParserGet::operator|| (const SxString &defaultValue) const
{
   if (symbol) {
      SxString res;
      try {
         res = symbol->toString ();
      } catch (SxException e)  {
         theParser.error ();
         e.print ();
         SX_EXIT;
      }
      return res;
   }
   return defaultValue;
}

/// get int from symbol or default value
int SxParserGet::operator|| (int defaultValue) const
{
   int res = defaultValue;
   if (symbol) {
      try {
         res = symbol->toInt ();
      } catch (SxException e)  {
         theParser.error ();
         e.print ();
         SX_EXIT;
      }
   }
   return res;
}
/// get int from symbol
SxParserGet::operator int () const
{
   try {
      if (symbol) return symbol->toInt ();
      theParser.current()->get (name); // raise error
   } catch (SxException e)  {
      theParser.error ();
      e.print ();
      SX_EXIT;
   }
   return -1;
}
/// get string from symbol
SxParserGet::operator SxString () const
{
   try {
      if (symbol) return symbol->toString ();
      theParser.current()->get (name); // raise error
   } catch (SxException e)  {
      theParser.error ();
      e.print ();
      SX_EXIT;
   }
   return -1;
}
/// get double from symbol
SxParserGet::operator double () const
{
   try {
      if (symbol) return symbol->toReal ();
      theParser.current()->get (name); // raise error
   } catch (SxException e)  {
      theParser.error ();
      e.print ();
      SX_EXIT;
   }
   return -1.;
}

/// get int array from symbol
SxParserGet::operator SxArray<int> () const
{
   try {
      if (symbol) return symbol->toIntList ();
      theParser.current()->get (name); // raise error
   } catch (SxException e)  {
      theParser.error ();
      e.print ();
      SX_EXIT;
   }
   return SxArray<int> ();
}

/// Direct access to the SxSymbol
const SxSymbol* SxParserGet::operator-> () const
{
   if (symbol) return symbol;
   try {
      return theParser.current()->get (name);
   } catch (SxException e)  {
      theParser.error ();
      e.print ();
      SX_EXIT;
   }
}

SxParserPush::SxParserPush (const char *groupName,
                            SxSimpleParser &parser,
                            const char* sourceFile,
                            int sourceLine)
   : theParser(parser)
{
   const SxSymbolTable *validGroup = theParser.currentGroup;
   if (!validGroup)  {
      // this may happen in constructs like
      // SYMBOLGROUP ("A") { ... } else SYMBOLGROUP ("B") { ... }
      // when symbolgroup "A" is not present
      int nInvalid = 0;
      // search for a valid group in the outer groups
      SxStack<const SxSymbolTable*> &outer = theParser.outerGroups;
      while (outer.getSize () > 0 && !(validGroup = outer.top ())) {
         outer.pop ();
         nInvalid++;
      }
      SX_CHECK (validGroup);
      // push back the invalid groups
      while (nInvalid--) outer.push (NULL);
   }
   // set location for error messages
   theParser.where (sourceFile, sourceLine);
   // put current group on stack
   theParser.outerGroups.push (theParser.currentGroup);
   // --- now try to get the requested group
   SxString name = groupName;
   theParser.currentGroup = validGroup->containsGroup (name)
                          ? validGroup->getGroup (name)
                          : NULL;
   // fallback: the first group may match the table we started with
   if (!theParser.currentGroup && theParser.outerGroups.getSize () == 1
       && theParser.outerGroups.top ()->name == name)
   {
      theParser.currentGroup = theParser.outerGroups.top ();
   }
}

SxParserPush::~SxParserPush ()  {
   theParser.currentGroup = theParser.outerGroups.pop ();
}
