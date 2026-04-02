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

#ifndef _SX_CACHE_OBLIVIOUS_H_
#define _SX_CACHE_OBLIVIOUS_H_

#include <SxMath.h>

/** \brief Small helper class for cache-oblivious algorithms

  The idea of a cache-oblivious algorithm is to recursively
  split a 2D problem such that there is a high chance that at
  some point, the problem fits entirely into a cache level.
  By doing it recursively, we do not need to know the cache
  size in advance.

  This class splits a larger problem into smaller bits
  of a certain maximum size (which should fit into the smallest
  cache). Since we use a finite internal stack, the total problem
  must first be split into large squares with a fixed chunk size,
  which are at most 2^15 = 32768x larger than the work size.

  Usage:
\code
   size_t blockSize = ... ; // less then 65535
   // --- outer loop: split nX,nY into chunks within uint16_t range
   for (size_t j = 0; j < nY; j+=blockSize)  {
      for (size_t i = 0; i < nX; i+=blockSize)  {
         uint16_t nx = uint16_t(max(blockSize, nX - i));
         uint16_t ny = uint16_t(max(blockSize, nY - j));

         SxCacheOblivious task(n, m, maxSize);
         while (task.splitProblem ())  {
            nx = task.nX;
            ny = task.nY;
            size_t offX = i + task.getOffX ();
            size_t offY = j + task.getOffY ();
            // work on problem nx x ny at offset (offX, offY)
         }
      }
   }
\endcode

    \author C. Freysoldt, freysoldt@mpie.de
*/
class SX_EXPORT_MATH SxCacheOblivious
{
   public:
      // These variables are set by splitProblem
      uint16_t nX, nY, offX, offY;
   protected:
      /// stack position (points to next empty slot)
      int top;
      /// max. work size
      int maxWork;
      /// The internal stack of unsolved problems
      struct { uint16_t offX, nX, offY, nY; } stack[16];
   public:
      /// Constructor
      SxCacheOblivious (uint16_t nX_, uint16_t nY_, int maxSize)
         : nX(0), nY(0), offX(0), offY(0),
           top (1), maxWork(maxSize)
      {
         SX_CHECK (((int(nX_) * int(nY_))>>15) <= maxSize, (int(nX_)*int(nY_)) >> 15, maxSize);
         stack[0].nX = uint16_t(nX_);
         stack[0].nY = uint16_t(nY_);
         stack[0].offX = stack[0].offY = 0;
      }
      /** \brief Set nX,nY,offX,offY to next problem
        @return Return false if done
        */
      inline bool splitProblem ();
};

bool SxCacheOblivious::splitProblem ()
{
   if (--top < 0) return false;
   // pop from stack
   nX = stack[top].nX;
   nY = stack[top].nY;
   offX = stack[top].offX;
   offY = stack[top].offY;
   do {
      if (int(nX) * int(nY) <= maxWork) {
         return true;
      } else {
         // --- split on longer dimension
         if (nX > nY) {
            // split on nX
            uint16_t newNX = uint16_t(nX>>1);
            // put upper half on stack
            stack[top].offX = uint16_t(offX + newNX);
            stack[top].nX   = uint16_t(nX - newNX);
            stack[top].offY = offY;
            stack[top].nY   = nY;
            // do lower half next
            nX = newNX;
         } else {
            // split on nY
            uint16_t newNY = uint16_t(nY>>1);
            // put right half on stack
            stack[top].offX = offX;
            stack[top].nX   = nX;
            stack[top].offY = uint16_t(offY + newNY);
            stack[top].nY   = uint16_t(nY - newNY);
            // do left half next
            nY = newNY;
         }
         top++;
      }
   } while (top < 16);
   if (int(nX) * int(nY) <= maxWork) return true;
   SX_EXIT;
}

#endif /* _SX_CACHE_OBLIVIOUS_H_ */
