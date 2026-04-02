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

#include <SxVector.h>
#include <string.h>

SxAllocMem* SxAllocMem::create (ssize_t n, int memTypeIn)
{
   SX_CHECK (n > 0, n);
   // allocate handle (separate cacheline)
   void *handle = SxAllocation::getAligned (roundCL(sizeof(SxAllocMem)),
                                            roundCL(1));
   //SX_CHECK (handle);
   // construct SxAllocMem on this memory
   return new (handle) SxAllocMem((size_t)n, memTypeIn);
}

inline
SxAllocMem::SxAllocMem (size_t n, int memTypeIn)
   : allocSize (n), memType(memTypeIn), refCounter (1)
{
   SX_CHECK (memType == 0, memType);
   // --- memtype = 0  => normal RAM
   mem = SxAllocation::get (n);
   //SX_CHECK (mem);
#  ifndef NDEBUG
      // fill with nan
      ::memset (mem, 0xff, allocSize);
#  endif
}
