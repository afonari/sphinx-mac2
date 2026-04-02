// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------
#include <SxUUIDv4.h>
#include <SxHashFunction.h>

SxUUIDv4::SxUUIDv4 ()
{
   ssize_t i;
   const ssize_t nBytes = 16;
   data.resize (nBytes);
   for (i=0; i < 16; ++i)  {  // generate 128 random bits
      data(i) = getRndByte ();
   }
   // --- RFC 4122, 4.4
   data(6) = 0x40 | (data(6) & 0x0f);
   data(8) = 0x80 | (data(8) & 0x3f);
}

SxUUIDv4::~SxUUIDv4 ()
{
   // empty
}

unsigned char SxUUIDv4::getRndByte () const
{
   std::uniform_int_distribution<unsigned int> dis (0, 255);
   return (unsigned char)dis (SxUtil::getGlobalObj ().mtEngine);
}

SxString SxUUIDv4::getStr () const
{
   SxString res; // "XXXXXXXX-XXXX-4XXX-XXXX-XXXXXXXXXXXX";
   for (ssize_t i = 0; i < 16; ++i) {
      res += SxString::sprintf ("%02X", data(i));
      if (i == 3 || i == 5 || i == 7 || i == 9) {
         res += "-";
      }
   }
   return res;
}

uint64_t SxUUIDv4::getHash () const
{
   return SxHashFunction::murmur64 ((const void *)data.elements,
                                    (int)data.getSize (), 12345);
}
