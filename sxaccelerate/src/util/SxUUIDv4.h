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

#ifndef _SX_UUID_V4_H_
#define _SX_UUID_V4_H_

#include <SxUtil.h>
#include <SxString.h>

class SX_EXPORT_UTIL SxUUIDv4
{
   public:
      SxUUIDv4 ();
     ~SxUUIDv4 ();

      SxString getStr () const;
      uint64_t getHash () const;

   protected:
      SxArray<unsigned char> data;

      unsigned char getRndByte () const;
};

#endif /* _SX_UUID_V4_H_ */
