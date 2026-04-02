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

#ifndef _SX_VEC_IT_H_
#define _SX_VEC_IT_H_

#include <SxMath.h>

// \brief Fundamental functionality of a vector iterator
// internal pointer that can be dereferenced and compared
template<class T>
class SxVecItBase
{
   protected:
      /// Pointer to element
      T* ptr;

      /// Empty constructor
      SxVecItBase() : ptr(NULL) {}
   public:
      /// Constructor
      SxVecItBase(T* in) : ptr(in) {}
      /// Dereference
      T& operator* () const forced_inline
      {
         SX_CHECK(ptr);
         return *ptr;
      }
      /// Comparison (equality)
      inline bool operator== (const SxVecItBase &it)  {
         return ( it.ptr == ptr );
      }
      /// Comparison (inequality)
      inline bool operator!= (const SxVecItBase &it)  {
         return ( it.ptr != ptr);
      }
};

template<class, StrideType> class SxVecConstIt;
template<class, StrideType> class SxVecRef;

// macro for implementing standard part of SxVecIt
// * friend declaration of SxVecRef<T,Layout>
// * Empty constructor
// * Constructor from SxVecIt base class
// * postfix operator
// * prefix operator
// * typecast operator to SxVecConstIt<T,Layout>&
// the last is just a declaration and absorbs the trailing ';'
#define SX_IT_IMPLEMENT(Layout) \
   SxVecIt () = default; \
   inline const SxVecIt<T,Layout> operator++ (int) forced_inline {\
      SxVecIt<T,Layout> copy (*this); \
      next (); \
      return copy; \
   } \
   inline SxVecIt<T,Layout> &operator++ () forced_inline { \
      next (); \
      return *this; \
   } \
   inline operator SxVecConstIt<T, Layout>& () \
   { return *reinterpret_cast<SxVecConstIt<T, Layout>*>(this); } \
   friend class SxVecRef<T,Layout>

// --- Iterator class (compact layout: Compact, PackedSymMatrix, BandedMatrix)
template<class T, StrideType Layout>
class SxVecIt : public SxVecItBase<T>
{
   protected:
      // How to increment
      void next () forced_inline { this->ptr++; }
      /// Constructor
      inline SxVecIt (T *from, const SxVecLayout<Layout> &)
         : SxVecItBase<T> (from) {}
      /// How to compute end
      static T* end (T *from, const SxVecLayout<Layout> &layout)
      {
         return from + layout.size;
      }
   public:
      /// Iterator from ptr
      inline SxVecIt (T *from) : SxVecItBase<T> (from) {}

      SX_IT_IMPLEMENT(Layout);
      /// Increment by n
      void operator+= (ssize_t n)  {
         this->ptr += n;
      }
};

// --- Iterator class (minimum functionality)
template<class T>
class SxVecIt<T, Strided> : public SxVecItBase<T>
{
   protected:
      ssize_t stride;
      // How to increment
      void next () forced_inline { this->ptr += stride; }
      // Constructor
      inline SxVecIt (T *from, const SxVecLayout<Strided> &layout)
         : SxVecItBase<T> (from), stride(layout.getRowStride ())
      {
         SX_CHECK (stride > 0, stride);
      }
      // How to compute end
      static T* end (T* start, const SxVecLayout<Strided> &layout)
      {
         return start + layout.size * layout.getRowStride ();
      }
      /// Constructor for end iterator
      inline SxVecIt (T *from) : SxVecItBase<T> (from) {}
   public:
      SX_IT_IMPLEMENT(Strided);
      /// Increment by n
      void operator+= (ssize_t n)  {
         this->ptr += n * stride;
      }
};

// --- Iterator class (minimum functionality)
template<class T>
class SxVecIt<T, SubMatrix> : public SxVecItBase<T>
{
   protected:
      ssize_t colStride;
      ssize_t nRows;
      T*      colEnd;
      // how to increment
      void next () forced_inline
      {
         this->ptr++;
         if (this->ptr == colEnd)  {
            this->ptr += colStride - nRows;
            colEnd += colStride;
         }
      }
      // How to compute end
      static T* end (T* start, const SxVecLayout<SubMatrix> &layout)
      {
         return start + layout.nCols * layout.getColStride ();
      }
      /// Constructor for end iterator
      inline SxVecIt (T *from) : SxVecItBase<T> (from) {}
   public:
      SX_IT_IMPLEMENT(SubMatrix);

      /// Initialize from ptr
      SxVecIt (T* in, const SxVecLayout<SubMatrix> &layout)
      : SxVecItBase<T>(in), colStride(layout.getColStride ()),
        nRows(layout.nRows), colEnd(in + nRows)
      {
         SX_CHECK (nRows > 0);
         SX_CHECK (nRows <= colStride);
      }
      /// Increment by n
      void operator+= (ssize_t n)  {
         // --- do column increment
         ssize_t nCols = n / nRows;
         this->ptr += nCols * colStride;
         colEnd += nCols * colStride;
         // --- do row increment
         n -= nCols * nRows;
         this->ptr += n;
         if (this->ptr >= colEnd)  {
            this->ptr += colStride - nRows;
            colEnd += colStride;
         }
      }
};

// --- Iterator class (minimum functionality)
template<class T>
class SxVecIt<T, GeneralMatrix> : public SxVecItBase<T>
{
   protected:
      ssize_t rowStride, colStride;
      ssize_t nRows;
      T* colEnd;
      // how to increment
      void next () forced_inline
      {
         this->ptr += rowStride;
         if (this->ptr == colEnd)  {
            this->ptr += colStride - nRows * rowStride;
            colEnd += colStride;
         }
      }
      // How to compute end
      static T* end (T* start, const SxVecLayout<GeneralMatrix> &layout)
      {
         return start + layout.nCols * layout.getColStride ();
      }
      /// Constructor for end iterator
      inline SxVecIt (T *from) : SxVecItBase<T> (from) {}
   public:
      SX_IT_IMPLEMENT(GeneralMatrix);

      /// Initialize from ptr
      SxVecIt (T* in, const SxVecLayout<GeneralMatrix> &layout)
      : SxVecItBase<T>(in),
        rowStride (layout.getRowStride ()), colStride(layout.getColStride ()),
        nRows(layout.nRows), colEnd(in + nRows * rowStride)
      {
         SX_CHECK (nRows > 0);
         SX_CHECK (rowStride > 0);
         SX_CHECK (colStride > 0);
         SX_CHECK (nRows * rowStride <= colStride);
      }

      /// Increment by n
      void operator+= (ssize_t n)  {
         // --- do column increment
         ssize_t nCols = n / nRows;
         this->ptr += nCols * colStride;
         colEnd += nCols * colStride;
         // --- do row increment
         n -= nCols * nRows;
         this->ptr += n * rowStride;
         if (this->ptr >= colEnd)  {
            this->ptr += colStride - nRows * rowStride;
            colEnd += colStride;
         }
      }

};

template<class T, StrideType Layout>
class SxVecConstIt : public SxVecIt<T, Layout>
{
   public:
      /// Add constness to dereferencing
      const T& operator* () forced_inline
      {
         return SxVecItBase<T>::operator* ();
      }
      /// Iterator from ptr
      inline SxVecConstIt (T *from) : SxVecIt<T, Layout> (from) {}
};

// ---


#endif /* _SX_VEC_IT_H_ */
