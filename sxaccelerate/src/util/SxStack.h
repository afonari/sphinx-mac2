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

#ifndef _SX_STACK_H_
#define _SX_STACK_H_
#include <SxError.h>
#include <string.h>

#ifdef __GNUC__
#define forced_inline __attribute__((always_inline))
#else
#define forced_inline
#endif

/** \brief Stack (push, pop) with blocked memory allocation

    \b SxClass = S/PHI/nX stack class

    This class provides a simple stack. You can push data onto
    the stack, and then pop up the most recently pushed data.

    This is much faster than SxList, but less flexible.
    The stack can be copied, but this is usually a bad idea.

    The stack can also be exported to an SxArray or SxVector.
    The ordering will be "as pushed", i.e. the first pushed
    data will be the first of the array. Once exported, the stack
    is empty again. This is done by
    \code
SxStack<double> stack;
...
// fill stack
...
SxArray<double> array (stack);
// now stack is empty
    \endcode

    The template parameter must have the copy constructor, for
    exporting to SxArray or SxVector also the assignment operator.

    For a simple class SxSimpleClass (no constructor,
    no destructor, no pointer members) one can declare
\code
    SXSTACK_SIMPLE(SxSimpleType) // no semicolon
\endcode
    to request that fast memcpy is used to export the stack.

    \note One block is always allocated, even if the stack
    is never used.

    \author Christoph Freysoldt, freysoldt@mpie.de */

template <class T>
class SxStack
{
   protected:
      /// Number of items to allocate per block
      size_t blockSize;

      /** \brief Memory block
        \note Singly linked list (to previous blocks)
       */
      class Block  {
         public:
            /// The actual data (managed by Block)
            T *elements;
            /// Block below current one
            Block *prev;
            /** Constructor
              @param size    size of the block
              @param current current top block, becomes prev
             */
            inline Block (size_t size, Block *current = NULL);
            /// Destructor
            ~Block ();
      };

      /// Top block of the stack
      Block *topBlock;
      /// Index in block
      size_t iElem;

      /** \brief Delete topmost block */
      inline void deleteTopBlock ()
      {
         Block *oldTopBlock = topBlock;
         topBlock = topBlock->prev;
         SX_CHECK (topBlock);
         delete oldTopBlock;
      }

      /// \name Auxiliary functions that may be specialized for simple types
      /// @{
      /// \brief Destroy elements of the top block
      inline void destroyElements ();
      ///@}

   public:
      /** \brief Constructor
        @param blockSizeIn the stack block size to use
      */
      SxStack (size_t blockSizeIn = sizeof(T) > 4096 ? 1 : (4096 / sizeof(T)) )
         : blockSize(blockSizeIn), iElem(0)
      {
         SX_CHECK(blockSize > 0);
         topBlock = new Block(blockSize);
      }

      /// Copy stack (usually a bad idea)
      void copy (const SxStack &in);

      /// Copy constructor (usually a bad idea)
      SxStack (const SxStack<T> &) { SX_EXIT; }
      /// Assignment operator (usually a bad idea)
      void operator= (const SxStack<T> &) { SX_EXIT; }

      /// removeAll stack
      void removeAll ();

      /** Destructor */
      ~SxStack () { removeAll (); delete topBlock; }

      /// \brief Copy-push data onto stack
      inline void push (const T& in) forced_inline
      {
         if (iElem == blockSize) {
            topBlock = new Block(blockSize, topBlock);
            iElem = 0;
         }
         // use placement new to copyconstruct in
         new (topBlock->elements + iElem++) T(in);
      }

      /// \brief Copy-push data onto stack
      inline SxStack<T>& operator<< (const T& in) forced_inline
      {
         push(in);
         return *this;
      }

      /// \brief Move-push data onto stack
      inline void push (T&& in) forced_inline
      {
         if (iElem == blockSize) {
            topBlock = new Block(blockSize, topBlock);
            iElem = 0;
         }
         // use placement new to move-construct in
         new (topBlock->elements + iElem++) T(std::move(in));
      }

      /// \brief Move-push data onto stack
      inline SxStack<T>&  operator<< (T&& in) forced_inline
      {
         push (std::move(in));
         return *this;
      }

      /// Pop up most recently pushed data
      inline T pop () forced_inline
      {
         SX_CHECK (iElem > 0 || topBlock->prev);
         if (iElem == 0)  {
            deleteTopBlock ();
            iElem = blockSize;
         }
         T* ptr = topBlock->elements + --iElem;
         // move element (this temporary is optimized away)
         T res (std::move(*ptr));
         // call destructor for popped element (cf. push)
         ptr-> ~T ();
         return res;
      }

      /// Access most recently pushed data
      inline T& top ()
      {
         SX_CHECK (iElem > 0 || topBlock->prev);
         if (iElem == 0)
            return topBlock->prev->elements[blockSize-1];
         return topBlock->elements[iElem - 1];
      }

      /** \brief Export stack data
        \param target Memory area where to copy the data to. Must be large
                      enough to hold the stack
        \param n      size of the stack
        Ordering: the first data pushed onto the stack will be
        first in the array.
        */
      void exportStack (T *target, size_t n) const&;

      /** \brief Export stack data
        \param target Memory area where to move the data to. Must be large
                      enough to hold the stack
        \param n      size of the stack
        Ordering: the first data pushed onto the stack will be
        first in the array.

        This function will move the elements in the stack.
        */
      void exportStack (T *target, size_t n) &&;

   public:

      /// \brief Return true if no data is on the stack
      inline bool isEmpty () const
      {
         return (! (iElem || topBlock->prev));
      }

      /** Compute the currently used stack size */
      size_t getSize () const;

      /** \brief Change block size (only empty stack!)
        @param newBlockSize the new block size
        */
      void changeBlockSize (size_t newBlockSize);
};

template<class T>
SxStack<T>::Block::Block (size_t size, Block *current)
: prev(current)
{
   SX_CHECK(size > 0);
   // allocate memory, but do not call any constructors.
   // use operator new rather than malloc to get correct alignment
   elements = static_cast<T*> (operator new (size * sizeof(T)));
}

template<class T>
SxStack<T>::Block::~Block ()
{
   // Deallocate memory consistently with constructor.
   operator delete (elements);
}

template <class T>
void SxStack<T>::removeAll ()
{
   while (true) {
      destroyElements ();
      if (!topBlock->prev) break; // this is the exit point.
      deleteTopBlock ();
      iElem = blockSize;
   }
}

// defined in SxArray.h
template<class T>
inline void sxCopyArray (T* dest, const T* src, ssize_t n);

template <class T>
void SxStack<T>::exportStack (T *target, size_t n) const&
{
   if (n == 0)  {
     SX_CHECK (iElem == 0 && !topBlock->prev);
     return;
   }
   SX_CHECK (n >= iElem, n, iElem);

   // --- copy data last (partial) block
   n -= iElem;
   sxCopyArray(target + n, topBlock->elements, iElem);

   // copy other blocks
   Block *block = topBlock;
   while (block->prev)  {
      block = block->prev;
      SX_CHECK (n >= blockSize, n, blockSize);
      n -= blockSize;
      sxCopyArray(target + n, block->elements, blockSize);
   }
   // make sure that stack has n elements
   SX_CHECK (n == 0, n);
}

template <class T>
void SxStack<T>::exportStack (T *target, size_t n) &&
{
   if (std::is_trivially_copyable<T>::value)  {
      static_cast<const SxStack<T>*>(this)->exportStack (target, n);
      removeAll ();
      return;
   }
   if (n == 0)  {
     SX_CHECK (this->isEmpty ());
     return;
   }
   // --- move data
   while (n)  {
      if (iElem == 0)  {
         deleteTopBlock ();
         iElem = blockSize;
         SX_CHECK(topBlock);
      }
      // move element
      target[--n] = std::move(topBlock->elements[--iElem]);
      // call destructor
      topBlock->elements[iElem].~T ();
   }
   // make sure that stack had n elements and that stack is empty now
   SX_CHECK(iElem == 0, iElem);
   SX_CHECK(!topBlock->prev);
}


template <class T>
inline void SxStack<T>::destroyElements ()
{
   if (std::is_trivially_destructible<T>::value)  {
      iElem = 0;
   } else {
      T* srcPtr = topBlock->elements;
      // count backwards, but destroy forward
      for (; iElem; --iElem, ++srcPtr) srcPtr->~T ();
   }
}

template <class T>
void SxStack<T>::copy (const SxStack<T> &in)
{
   size_t n = in.iElem;
   T *dest, *src;
   Block *destBlock = topBlock, *srcBlock = in.topBlock;
   for (; srcBlock; srcBlock = srcBlock->prev)  {
      src  = srcBlock->elements;
      dest = destBlock->elements;
      // use placement new to copyconstruct elements on stack
      for (; n; --n) new (src++) T(*dest++);

      if (srcBlock->prev)  {
         // allocate new Block before(!) current one
         destBlock->prev = new Block(blockSize);
         n = blockSize;
      }
   }
}

template <class T>
size_t SxStack<T>::getSize () const
{
   size_t size = iElem;
   for (Block *ptr = topBlock->prev; ptr; ptr = ptr->prev) size += blockSize;
   return size;
}

template <class T>
void SxStack<T>::changeBlockSize (size_t newBlockSize)
{
   SX_CHECK (newBlockSize > 0);
   SX_CHECK (iElem == 0 && !(topBlock->prev));
   delete topBlock;
   blockSize = newBlockSize;
   topBlock = new Block(blockSize);
}

#undef forced_inline

#endif /* _SX_STACK_H_ */
