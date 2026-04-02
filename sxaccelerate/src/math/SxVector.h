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

#ifndef _SX_VECTOR_H_
#define _SX_VECTOR_H_

#include <SxMath.h>
#include <SxComplex.h>
#include <SxArray.h>
#include <SxIdx.h>
#include <SxBlasLib.h>
#include <atomic>
#include <SxAllocCache.h>
#include <SxVector3.h>
#include <SxMatrix3.h>
#include <SxRandom.h>

#include <SxTypeMapper.h>
#include <SxVecLayout.h>
#include <SxVecIt.h>
#include <SxVecCompute.h>

class SxBasis;
ssize_t getNElements (const SxBasis &b);

class SX_EXPORT_MATH SxAuxData
{
   public:
      const SxBasis* basisPtr;
      char iSpin, n, l, m;
      unsigned char nComp;
      int ik;
      int is;
      int ia;
      int iOther[4]; // can be extended to more, if needed

      /// Set from y if only y has basisPtr set, otherwise from x
      inline void setCombined (const SxAuxData &x, const SxAuxData &y)
      {
         if (x.basisPtr || !y.basisPtr)  {
            operator= (x);
         } else {
            operator= (y);
         }
      }

      /// Empty constructor
      SxAuxData ()
         : basisPtr (NULL), iSpin(-1), n(-1), l(-1), m(-1),
           nComp(1), ik(-1), is(-1), ia(-1)
      {
         iOther[0] = iOther[1] = iOther[2] = iOther[3] = -1;
      }

      /// Combine constructor
      SxAuxData (const SxAuxData aux1, const SxAuxData aux2)
      {
         setCombined (aux1, aux2);
      }

      /// Unset all data
      inline void unset ()
      {
         basisPtr = NULL;
         iSpin = n = l = m = -1;
         nComp = 1;
         ik = is = ia = -1;
         iOther[0] = iOther[1] = iOther[2] = iOther[3] = -1;
      }
};

// --- descriptor of allocated memory (incl. reference counting)
class SxAllocMem
{
   public:
      /// Allocated memory
      void *mem;
      /// Allocated size (bytes)
      size_t allocSize;
      /// Type of memory (currently unused)
      const int memType;
      /// Reference counter
      std::atomic<int> refCounter;

      /// Increase reference
      void ref () {
         SX_CHECK (mem);
         refCounter++;
      }

      /** \brief Decrease reference
          If true, this allocation can be destroyed
        */
      bool unref () {
         refCounter--;
         if (refCounter == 0)  {
            SX_CHECK (memType == 0);
            SxAllocation::retain (mem, allocSize);
            mem = NULL;
            return true;
         }
         return false;
      }

      template<class T>
      bool isInside (const T* start, ssize_t n) const
      {
         if (sizeof(T) * n > allocSize) return false;
         return (mem <= start && (start + n) <= (void*)((char*)mem + allocSize));
      }


      static constexpr size_t roundCL (size_t n)
      {
         const size_t cacheLine = 64;
         return ((n & (cacheLine - 1)) ? cacheLine : 0)
               + (n & ~(cacheLine - 1));
      }

      /// Create allocation (already referenced)
      static SxAllocMem *create (ssize_t n, int memTypeIn = 0);
      /// Destroy "created" allocation (return memory to pool)
      static void destroy (SxAllocMem *alloc)
      {
         SX_CHECK (alloc->refCounter == 0);
         SX_CHECK (alloc->mem == NULL);
         SxAllocation::retain (alloc, roundCL(sizeof(SxAllocMem)));
      }
   protected:
      /// Constructor (already referenced)
      SxAllocMem (size_t n, int memTypeIn = 0);
};

inline
std::ostream& operator<< (std::ostream& out, const SxAllocMem &allocMem)
{
   void *memEnd =(void*)((char*)allocMem.mem+allocMem.allocSize-1);
   out << "mem@" << allocMem.mem << ".." << memEnd
       << " (size " << allocMem.allocSize << ")";
   return out;
}

#ifdef NDEBUG
#   define SX_VALIDATE_VECTOR(expr) ((void)0);
#else
template<class T, StrideType Layout> class SxVecRef;
template<class T, StrideType Layout>
inline void SX_VALIDATE_VECTOR (const SxVecRef<T,Layout>&);
#endif /* NDEBUG */


// ---
template<class T> class SxVector;

template<bool Packed> class SxVectorReIm;

template<class T, StrideType Layout = Compact>
class SxVecRef : protected SxVecLayout<Layout>
{
   public:
      // provide scalar type as typedef
      typedef T ScalarType;

      // SxVector is a friend
      friend class SxVector<T>;
      // all SxVecRef's are friends
      template <class T2, StrideType Layout2>
      friend class SxVecRef;

   protected:
      /// Pointer to allocated memory handle
      SxAllocMem *allocMem;
   public:
      // --- data
      /// Pointer to first data elementsent
      T* elements;


   public:

      /// The auxiliary data
      SxAuxData auxData;

      /// Get the basis pointer
      inline const SxBasis *getBasisPtr () const
      {
         return auxData.basisPtr;
      }

      /// Get basis reference of specific type
      template<class B>
      const B& getBasis () const
      {
         const B* ptr = dynamic_cast<const B*>(auxData.basisPtr);
         SX_CHECK (ptr);
         return *ptr;
      }

      /// Set the basis pointer
      inline void setBasis (const SxBasis *basis)
      {
         auxData.basisPtr = basis;
      }

      /// Set the basis pointer
      inline void setBasis (const SxBasis &basis)
      {
         auxData.basisPtr = &basis;
      }

      // --- Constructors
      /// Empty constructor
      SxVecRef () : allocMem(NULL), elements(NULL) {}

      /// "Copy" constructor -> references data
      SxVecRef (SxVecRef<T,Layout> &);

      /// Avoid removing const-ness by creating a referencing non-const copy
      SxVecRef (const SxVecRef<T,Layout> &) = delete;

      /// Copy constructor from different type (copies)
      template<class T2>
      SxVecRef (const SxVecRef<T2, Layout> &);

      /// Layout cast constructor -> references data
      template<StrideType Layout2>
      SxVecRef (SxVecRef<T,Layout2> &);

      /// Move constructor
      SxVecRef (SxVecRef<T,Layout> &&);
   protected:
      /// Set size, but no memory
      SxVecRef (ssize_t nRowsIn, ssize_t rowStrideIn = 1,
                ssize_t nColsIn = 1, ssize_t colStrideIn = 1)
         : SxVecLayout<Layout> (nRowsIn, nColsIn, rowStrideIn, colStrideIn),
           allocMem(NULL), elements(NULL)
      { /* empty */ }

      /** \brief Allocate memory and set elements
          \note This is an low-level auxiliary function for - no automatic
                deallocation/no checks
        */
      inline void alloc (ssize_t n)
      {
         allocMem = SxAllocMem::create (sizeof(T) * n);
         elements = (T*)allocMem->mem;
      }

      /// Constructor for result from binary vector-vector operations
      template<StrideType Layout2>
      SxVecRef (SxAllocMem* mem, const SxVecLayout<Layout2> &layout,
                const SxAuxData &aux1, const SxAuxData &aux2)
         : SxVecLayout<Layout> (layout),
           allocMem(mem), auxData (aux1, aux2)
      {
         SX_CHECK(mem);
         elements = (T*)mem->mem;
         SX_VALIDATE_VECTOR (*this);
      }

      /// Lease association without clearing allocMem/elements/SxVecLayout
      inline void unref_ ()
      {
         // free memory (inside allocMem and for allocMem handle)
         if (allocMem && allocMem->unref ()) SxAllocMem::destroy (allocMem);
      }

   public:
      /// Destructor
      inline ~SxVecRef ()
      {
         unref_ ();
      }

      /// Typecast to SxMatrix3 (vector must be 3x3 matrix)
      inline operator SxMatrix3<T> ();

      /// Change to 3-vector
      SxVector3<T> toVector3 () const
      {
         SX_CHECK (getSize () == 3);
         return SxVector3<T> ((*this)(0), (*this)(1), (*this)(2));
      }

   public:
      /// Get number of elements
      ssize_t getSize () const { return this->size; }
      /// Get number of rows
      ssize_t getNRows () const { return this->nRows; }
      /// Get number of columns
      ssize_t getNCols () const { return this->nCols; }
      /// Get row stride
      using SxVecLayout<Layout>::getRowStride;
      /// Get column stride
      using SxVecLayout<Layout>::getColStride;
      /// Get memory offset
      using SxVecLayout<Layout>::getIdx;

   protected:
      // ensure visibility of nRows and nCols inside our functions
      using SxVecLayout<Layout>::nRows;
      using SxVecLayout<Layout>::nCols;

   public:
      // --- associating (use with care!)
      SxVecRef (SxAllocMem* allocIn,
                T*          data,
                ssize_t     nRowsIn,
                ssize_t     rowStrideIn = 1,
                ssize_t     nColsIn = 1,
                ssize_t     colStrideIn = 1);

   protected:
      // --- associating (use with care!)
      SxVecRef (SxAllocMem* allocIn,
                T*          data,
                ssize_t     nRowsIn,
                ssize_t     rowStrideIn,
                ssize_t     nColsIn,
                ssize_t     colStrideIn,
                const SxAuxData &auxIn);
   protected:
      SxVecRef (const SxVecLayout<Layout> &layout,
                const SxAuxData &auxIn)
         : SxVecLayout<Layout>(layout),
           allocMem (SxAllocMem::create (layout.size * sizeof(T))),
           elements ((T*)allocMem->mem),
           auxData(auxIn)
      {
         // empty
      }

   public:
      /// Invalidate this reference (and maybe free allocation)
      inline void unref ()  {
         unref_ ();
         allocMem = NULL;
         elements = NULL;
         SxVecLayout<Layout>::unset ();
         auxData.unset ();
      }

      /// True if no other references to memory exist
      inline bool isOnlyRef () const
      {
         return allocMem ? (allocMem->refCounter == 1) : true;
      }

      // --- assignment
      /** \brief Copy assignment (associates if unassociated, otherwise sets)
        \note Assigning an empty vector will remove the current association
        */
      inline SxVecRef<T,Layout>& operator= (const SxVecRef<T,Layout>& in);

      /// Copy assignment (sets)
      template<StrideType Layout2>
      inline SxVecRef<T,Layout>& operator<<= (const SxVecRef<T,Layout2>& in);

      /// Copy assignment (sets)
      inline SxVecRef<T,Layout>& operator<<= (const SxArray<T>& in);

      /// Assignment from scalar (sets elements)
      inline void operator= (const T &x)
      {
         SX_CHECK_NUM(x);
         set (x);
      }

      /// Assignment from different SxVecRef (copies)
      template<class T2, StrideType Layout2>
      inline SxVecRef<T,Layout>& operator= (const SxVecRef<T2, Layout2> &);

      /** \brief Assignment from different SxVecRef (passes on memory if possible)

          @note: This takes over the memory if the receiving layout is compatible
                 with the actual layout, the receiving reference is not associated,
                 and the incoming reference is not referenced multiple times
        */
      template<StrideType Layout2>
      inline SxVecRef<T,Layout>& operator= (SxVecRef<T, Layout2> &&);

      /// Set all values
      inline void set (const T&);
      /// Set from array
      inline void set (const SxArray<T>&);

      /// Reorder data (in-place)
      void sortByIdx (const SxArray<ssize_t> &sortIdx);

      // --- iterators
      typedef SxVecIt<T, Layout>      Iterator;
      typedef SxVecConstIt<T, Layout> ConstIterator;

      /// Get begin iterator
      Iterator      begin ()       { return SxVecIt<T,Layout> (elements, *this); }
      ConstIterator begin () const { return SxVecIt<T,Layout> (const_cast<T*>(elements), *this); }
      /// Get end iterator
      Iterator      end ()       { return SxVecIt<T,Layout>::end (elements, *this); }
      ConstIterator end () const { return SxVecIt<T,Layout>::end (elements, *this); }


      // --- access elements (non-const)
      inline T& operator() (ssize_t i)
      {
         SX_CHECK(i >= 0 && i < getSize (), i, getSize ());
         // only non-modified vectors can give direct access to elements
         return elements[SxVecLayout<Layout>::getIdx (i)];
      }
      inline T& operator() (SxAutoLoop &i) forced_inline
      {
         i.setLimit (getSize ());
         return operator() ((ssize_t)i);
      }
      /// Get matrix elementsent
      inline T& operator() (ssize_t r, ssize_t c)
      {
         SX_CHECK (r >= 0 && r < getNRows (), r, getNRows ());
         SX_CHECK (c >= 0 && c < getNCols (), c, getNCols ());
         return elements[SxVecLayout<Layout>::getIdx (r,c)];
      }

      // autoloop wrapper (I1, I2 can be int, ssize_t, or SxAutoLoop)
      template<class I1, class I2>
      inline T& operator() (const I1 &r, const I2 &c)
      {
         SxAutoLoop::setLimit (r, getNRows ());
         SxAutoLoop::setLimit (c, getNCols ());
         return operator() ((ssize_t)r, (ssize_t)c);
      }

      /// Reference to column
      inline SxVecRef<T, colStrideType(Layout)> colRef (ssize_t c);
      inline SxVecRef<T, colStrideType(Layout)> colRef (SxAutoLoop &c)
      {
         c.setLimit (getNCols ());
         return colRef((ssize_t)c);
      }

      /// Reference to several columns
      inline SxVecRef<T, Layout> colsRef (ssize_t c, ssize_t nColIn)
      {
         SX_CHECK (!isPacked(Layout));
         SX_CHECK (c >= 0 && c + nColIn <= getNCols (), c, nColIn, getNCols ());
         return SxVecRef<T, Layout> (allocMem, elements + c * getColStride (),
                                     nRows, getRowStride (),
                                     nColIn, getColStride (), auxData);
      }

      /// Reference to row
      inline SxVecRef<T, Strided> rowRef (ssize_t r);
      inline SxVecRef<T, Strided> rowRef (SxAutoLoop &r)
      {
         r.setLimit (getNRows ());
         return rowRef((ssize_t)r);
      }

      /// Reference to diagonal (of square matrix)
      inline SxVecRef<T, Strided> diag ();

      /// Reference to a component
      inline SxVecRef<T, subMatrixLayout(Layout)> compRef (int iComp);

      /** \brief Reference to subelements
        @param offset offset into current vector
        @param nRowsIn number of rows in new reference
        @param nColsIn number of cols in new reference
        @param rStride row stride (absolute, i.e. with respect to memory),
               Must be multiple of current row stride.
        @param cStride (if given and > 1): col stride (absolute, wrt to memory)
               Must be multiple of current column stride

        @note: This returns a modifiable reference to a (possibly) const object
        @note: Use with care!
        */
      template<StrideType Layout2>
      SxVecRef<T, Layout2> getRef (ssize_t offset,
                                   ssize_t nRowsIn, ssize_t nColsIn,
                                   ssize_t rStride = 1, ssize_t cStride = -1) const;

      /// Reference to subelements
      inline SxVecRef<T, Layout> getRef (ssize_t newSize, ssize_t offset = 0)
      {
         SX_CHECK (offset + newSize <= getSize (), offset, newSize, getSize ());
         // make sure that we can safely reference multiple columns, or that
         // the reference stays within one column
         SX_CHECK (getRowStride () * getNRows () == getColStride ()
                   || ((offset % getColStride ())
                       + newSize * getRowStride () <= getColStride ()));
         return getRef<Layout> (offset, newSize, 1,
                                getRowStride (), getColStride ());
      }

      /// Reference to subelements
      inline const SxVecRef<T, Layout>
      getRef (ssize_t newSize, ssize_t offset = 0) const
      {
         SX_CHECK (offset + newSize <= getSize (), offset, newSize, getSize ());
         // make sure that we can safely reference multiple columns, or that
         // the reference stays within one column
         SX_CHECK (getRowStride () * getNRows () == getColStride ()
                   || ((offset % getColStride ())
                       + newSize * getRowStride () <= getColStride ()));
         return getRef<Layout> (offset, newSize, 1,
                                getRowStride (), getColStride ());
      }

      // --- access elements (const)
      inline const T& operator() (ssize_t i) const
      {
         SX_CHECK(i >= 0 && i < getSize (), i, getSize ());
         return elements[SxVecLayout<Layout>::getIdx (i)];
      }
      inline const T& operator() (SxAutoLoop &i) const
      {
         i.setLimit (getSize ());
         return operator() ((ssize_t)i);
      }

      /// Get matrix element
      inline const T& operator() (ssize_t r, ssize_t c) const
      {
         SX_CHECK (r >= 0 && r < getNRows (), r, getNRows ());
         SX_CHECK (c >= 0 && c < getNCols (), c, getNCols ());
         return elements[SxVecLayout<Layout>::getIdx (r,c)];
      }

      // autoloop wrapper (I1, I2 can be int, ssize_t, or SxAutoLoop)
      template<class I1, class I2>
      inline const T& operator() (const I1 &r, const I2 &c) const
      {
         SxAutoLoop::setLimit (r, getNRows ());
         SxAutoLoop::setLimit (c, getNCols ());
         return operator() ((ssize_t)r, (ssize_t)c);
      }

      /// Reference to column
      inline const SxVecRef<T, colStrideType(Layout)> colRef (ssize_t c) const;
      inline const SxVecRef<T, colStrideType(Layout)> colRef (SxAutoLoop& c) const
      {
         c.setLimit (getNCols ());
         return const_cast<SxVecRef<T,Layout>*>(this)->colRef((ssize_t)c);
      }

      /// Reference to several columns
      inline const SxVecRef<T, Layout> colsRef (ssize_t c, ssize_t nColIn) const
      {
         SX_CHECK (!isPacked(Layout));
         SX_CHECK (c >= 0 && c + nColIn <= getNCols (), c, nColIn, getNCols ());
         return SxVecRef<T, Layout> (allocMem, elements + c * getColStride (),
                                     nRows, getRowStride (),
                                     nColIn, getColStride (), auxData);
      }

      /// Reference to row
      inline const SxVecRef<T, Strided> rowRef (ssize_t r) const;
      inline const SxVecRef<T, Strided> rowRef (SxAutoLoop& r) const
      {
         r.setLimit (getNRows ());
         return rowRef((ssize_t)r);
      }

      /// Reference to diagonal (of square matrix)
      inline const SxVecRef<T, Strided> diag () const;

      /// Reference to a component
      inline const SxVecRef<T, subMatrixLayout(Layout)>
      compRef (int iComp) const;

      typedef typename SxTypeMapper<T>::TReal TReal;
      friend class SxVectorReIm<false>;

      /// Real part (only for complex vectors!)
      inline const SxVecRef<TReal,realPartLayout(Layout)> real () const
      {
         return SxVectorReIm<isPacked(Layout)>::real (*this);
      }

      /// Imaginary part (only for complex vectors!)
      inline const SxVecRef<TReal,realPartLayout(Layout)> imag () const
      {
         return SxVectorReIm<isPacked(Layout)>::imag (*this);
      }

      // --- change shape
      SxVecRef<T,Layout>&
      reshape (ssize_t nRowsIn, ssize_t nColsIn = 1);

      // --- arithmetic operations
      template<SxVecAlgorithm> friend class SxVecCompute;

      // --- The next 136 lines were generated from snippets/SxVector.h snippet vec_op_vec
      // vector x + y (x, y have same layout)
      template<class T2>
      typename SxVecResult<decltype(T(0.) + T2(0.)), Layout>::VecType
      operator+ (const SxVecRef<T2, Layout> &x) const&
      {
         return SxVecCompute<algoSelect(std::is_same<T, T2>::value, Layout, Layout)>
                ::add (*this, x);
      }

      // vector x + y (x, y have different layout)
      template<class T2, StrideType Layout2>
      SxVector<decltype(T(0.) + T2(0.))>
      operator+ (const SxVecRef<T2, Layout2> &x) const
      {
         return SxVecCompute<algoSelect(std::is_same<T, T2>::value, Layout, Layout2)>
                ::add (*this, x);
      }

      // vector x + y (x, y have different layout, x is about to die)
      template<StrideType Layout2>
      typename SxVecResult<T, Layout>::VecType
      operator+ (const SxVecRef<T, Layout2> &y) &&
      {
         SX_CHECK (allocMem);
         if (isOnlyRef () && isContiguous(Layout))  {
            //cout << "MOVE add" << endl;
            SxVecCompute<algoInPlace(Layout, Layout2)>::addInPlace (*this, y);
            auxData.setCombined (auxData, y.auxData);
            return std::move(*this);
         } else {
            return SxVecCompute<algoSelect(true, Layout, Layout2)>::add (*this, y);
         }
      }

      // vector x - y (x, y have same layout)
      template<class T2>
      typename SxVecResult<decltype(T(0.) + T2(0.)), Layout>::VecType
      operator- (const SxVecRef<T2, Layout> &x) const&
      {
         return SxVecCompute<algoSelect(std::is_same<T, T2>::value, Layout, Layout)>
                ::subtract (*this, x);
      }

      // vector x - y (x, y have different layout)
      template<class T2, StrideType Layout2>
      SxVector<decltype(T(0.) + T2(0.))>
      operator- (const SxVecRef<T2, Layout2> &x) const
      {
         return SxVecCompute<algoSelect(std::is_same<T, T2>::value, Layout, Layout2)>
                ::subtract (*this, x);
      }

      // vector x - y (x, y have different layout, x is about to die)
      template<StrideType Layout2>
      typename SxVecResult<T, Layout>::VecType
      operator- (const SxVecRef<T, Layout2> &y) &&
      {
         SX_CHECK (allocMem);
         if (isOnlyRef () && isContiguous(Layout))  {
            //cout << "MOVE subtract" << endl;
            SxVecCompute<algoInPlace(Layout, Layout2)>::subtractInPlace (*this, y);
            auxData.setCombined (auxData, y.auxData);
            return std::move(*this);
         } else {
            return SxVecCompute<algoSelect(true, Layout, Layout2)>::subtract (*this, y);
         }
      }

      // vector x * y (x, y have same layout)
      template<class T2>
      typename SxVecResult<decltype(T(0.) + T2(0.)), Layout>::VecType
      operator* (const SxVecRef<T2, Layout> &x) const&
      {
         return SxVecCompute<algoSelect(std::is_same<T, T2>::value, Layout, Layout)>
                ::multiply (*this, x);
      }

      // vector x * y (x, y have different layout)
      template<class T2, StrideType Layout2>
      SxVector<decltype(T(0.) + T2(0.))>
      operator* (const SxVecRef<T2, Layout2> &x) const
      {
         return SxVecCompute<algoSelect(std::is_same<T, T2>::value, Layout, Layout2)>
                ::multiply (*this, x);
      }

      // vector x * y (x, y have different layout, x is about to die)
      template<StrideType Layout2>
      typename SxVecResult<T, Layout>::VecType
      operator* (const SxVecRef<T, Layout2> &y) &&
      {
         SX_CHECK (allocMem);
         if (isOnlyRef () && isContiguous(Layout))  {
            //cout << "MOVE multiply" << endl;
            SxVecCompute<algoInPlace(Layout, Layout2)>::multiplyInPlace (*this, y);
            auxData.setCombined (auxData, y.auxData);
            return std::move(*this);
         } else {
            return SxVecCompute<algoSelect(true, Layout, Layout2)>::multiply (*this, y);
         }
      }

      // vector x / y (x, y have same layout)
      template<class T2>
      typename SxVecResult<decltype(T(0.) + T2(0.)), Layout>::VecType
      operator/ (const SxVecRef<T2, Layout> &x) const&
      {
         return SxVecCompute<algoSelect(std::is_same<T, T2>::value, Layout, Layout)>
                ::divide (*this, x);
      }

      // vector x / y (x, y have different layout)
      template<class T2, StrideType Layout2>
      SxVector<decltype(T(0.) + T2(0.))>
      operator/ (const SxVecRef<T2, Layout2> &x) const
      {
         return SxVecCompute<algoSelect(std::is_same<T, T2>::value, Layout, Layout2)>
                ::divide (*this, x);
      }

      // vector x / y (x, y have different layout, x is about to die)
      template<StrideType Layout2>
      typename SxVecResult<T, Layout>::VecType
      operator/ (const SxVecRef<T, Layout2> &y) &&
      {
         SX_CHECK (allocMem);
         if (isOnlyRef () && isContiguous(Layout))  {
            //cout << "MOVE divide" << endl;
            SxVecCompute<algoInPlace(Layout, Layout2)>::divideInPlace (*this, y);
            auxData.setCombined (auxData, y.auxData);
            return std::move(*this);
         } else {
            return SxVecCompute<algoSelect(true, Layout, Layout2)>::divide (*this, y);
         }
      }
      // --- vec_op_vec

      // --- The next 72 lines were generated from snippets/SxVector.h snippet vec_op_eq
      /// add vectors in place
      template<StrideType Layout2>
      SxVecRef<T, Layout>& operator+= (const SxVecRef<T, Layout2> &y)
      {
         SxVecCompute<algoInPlace(Layout,Layout2)>::addInPlace (*this, y);
         return *this;
      }

      /// add vectors in place with type mixing
      template<class T2, StrideType Layout2>
      SxVecRef<T, Layout>& operator+= (const SxVecRef<T2, Layout2> &y)
      {
         // no in-place with automatic type (down)conversion
         SX_CHECK ((std::is_same<T, decltype(T(0) + T2(0))>::value));
         SxVecCompute<UseIterator>::addInPlace (*this, y);
         return *this;
      }

      /// subtract vectors in place
      template<StrideType Layout2>
      SxVecRef<T, Layout>& operator-= (const SxVecRef<T, Layout2> &y)
      {
         SxVecCompute<algoInPlace(Layout,Layout2)>::subtractInPlace (*this, y);
         return *this;
      }

      /// subtract vectors in place with type mixing
      template<class T2, StrideType Layout2>
      SxVecRef<T, Layout>& operator-= (const SxVecRef<T2, Layout2> &y)
      {
         // no in-place with automatic type (down)conversion
         SX_CHECK ((std::is_same<T, decltype(T(0) - T2(0))>::value));
         SxVecCompute<UseIterator>::subtractInPlace (*this, y);
         return *this;
      }

      /// multiply vectors in place
      template<StrideType Layout2>
      SxVecRef<T, Layout>& operator*= (const SxVecRef<T, Layout2> &y)
      {
         SxVecCompute<algoInPlace(Layout,Layout2)>::multiplyInPlace (*this, y);
         return *this;
      }

      /// multiply vectors in place with type mixing
      template<class T2, StrideType Layout2>
      SxVecRef<T, Layout>& operator*= (const SxVecRef<T2, Layout2> &y)
      {
         // no in-place with automatic type (down)conversion
         SX_CHECK ((std::is_same<T, decltype(T(0) * T2(0))>::value));
         SxVecCompute<UseIterator>::multiplyInPlace (*this, y);
         return *this;
      }

      /// divide vectors in place
      template<StrideType Layout2>
      SxVecRef<T, Layout>& operator/= (const SxVecRef<T, Layout2> &y)
      {
         SxVecCompute<algoInPlace(Layout,Layout2)>::divideInPlace (*this, y);
         return *this;
      }

      /// divide vectors in place with type mixing
      template<class T2, StrideType Layout2>
      SxVecRef<T, Layout>& operator/= (const SxVecRef<T2, Layout2> &y)
      {
         // no in-place with automatic type (down)conversion
         SX_CHECK ((std::is_same<T, decltype(T(0) / T2(0))>::value));
         SxVecCompute<UseIterator>::divideInPlace (*this, y);
         return *this;
      }
      // --- vec_op_eq

      // Y += a*X
      template<StrideType Layout2>
      inline void plus_assign_ax (const T &a,
                                  const SxVecRef<T,Layout2> &x)
      {
         SX_CHECK(this->getSize () == x.getSize (),
                  this->getSize (), x.getSize ());
         SX_CHECK(this->getSize () > 0);
         SxVecCompute<algoInPlace(Layout,Layout2)>::axpy (*this, a, x);
         SX_VALIDATE_VECTOR (*this);
      }

      // Y += a*X
      template<class T2, StrideType Layout2>
      inline void plus_assign_ax (const T &a,
                                  const SxVecRef<T2,Layout2> &x)
      {
         SX_CHECK(this->getSize () == x.getSize (),
                  this->getSize (), x.getSize ());
         SX_CHECK(this->getSize () > 0);
         SxVecCompute<UseIterator>::axpy (*this, a, x);
         SX_VALIDATE_VECTOR (*this);
      }

      // Matrix-matrix multiplication
      template<StrideType Layout2>
      inline SxVector<T> operator^ (const SxVecRef<T,Layout2> &x) const
      {
         SX_CHECK (this->getSize () > 0);
         SX_CHECK (x.getSize () > 0);
         if (isPacked (Layout) || this->getNCols () > 1 || x.getNRows () == 1) {
            // standard case
            SX_CHECK (this->getNCols () == x.getNRows (),
                      this->getNCols (), x.getNRows ());
            return SxVecCompute<algoMat(Layout,Layout2)>::matmult (*this, x);
         }
         // special case: col vector ^ matrix
         // => vector will be interpreted as row vector
         SX_CHECK (this->getNRows () == x.getNRows (),
                   this->getNRows (), x.getNRows ());
         if (this->getRowStride () == 1)  {
            SxVecRef<T, Compact> rowRef
               = this->getRef<Compact> (0, 1, this->getNRows ());
            SxVector<T> res = SxVecCompute<algoMat(Compact,Layout2)>::matmult (rowRef, x);
            res.reshape (x.getNCols (), 1);
            return res;
         } else {
            SxVecRef<T, Strided> rowRef
               = this->getRef<Strided> (0, 1, this->getNRows (),
                                        this->getRowStride ());
            return SxVecCompute<algoMat(Strided,Layout2)>::matmult (rowRef, x);
         }
      }

      // --- The next 172 lines were generated from snippets/SxVector.h snippet scalar_op
      inline typename SxVecResult<decltype(T(0) + typename SxTypeMapper<T>::TReal(0)), Layout>::VecType
      operator+ (typename SxTypeMapper<T>::TReal x) const&
      {
         return SxVecCompute<algoScalar(Layout)>::add (*this, x);
      }

      inline typename SxVecResult<typename SxTypeMapper<T>::TCmplx, Layout>::VecType
      operator+ (const typename SxTypeMapper<T>::TCmplx &x) const&
      {
         return SxVecCompute<algoScalar(Layout)>::add (*this, x);
      }

      inline typename SxVecResult<decltype(T(0) + typename SxTypeMapper<T>::TReal(0)), Layout>::VecType
      operator+ (typename SxTypeMapper<T>::TReal x) &&
      {
         SX_CHECK (allocMem);
         if (isContiguous(Layout) && isOnlyRef ())
         {
            return SxVecRef<T, Layout>(std::move (*this)) += x;
         } else {
            return SxVecCompute<algoScalar(Layout)>::add (*this, x);
         }
      }

      inline typename SxVecResult<typename SxTypeMapper<T>::TCmplx, Layout>::VecType
      operator+ (const typename SxTypeMapper<T>::TCmplx &x) &&
      {
         SX_CHECK (allocMem);
         if (isContiguous(Layout) && isOnlyRef ())
         {
            return SxVecRef<T, Layout>(std::move (*this)) += x;
         } else {
            return SxVecCompute<algoScalar(Layout)>::add (*this, x);
         }
      }

      inline SxVecRef<T, Layout>&
      operator+= (const T &x)
      {
         SxVecCompute<algoScalar(Layout)>::addInPlace (*this, x);
         return *this;
      }

      inline typename SxVecResult<decltype(T(0) - typename SxTypeMapper<T>::TReal(0)), Layout>::VecType
      operator- (typename SxTypeMapper<T>::TReal x) const&
      {
         return SxVecCompute<algoScalar(Layout)>::subtract (*this, x);
      }

      inline typename SxVecResult<typename SxTypeMapper<T>::TCmplx, Layout>::VecType
      operator- (const typename SxTypeMapper<T>::TCmplx &x) const&
      {
         return SxVecCompute<algoScalar(Layout)>::subtract (*this, x);
      }

      inline typename SxVecResult<decltype(T(0) - typename SxTypeMapper<T>::TReal(0)), Layout>::VecType
      operator- (typename SxTypeMapper<T>::TReal x) &&
      {
         SX_CHECK (allocMem);
         if (isContiguous(Layout) && isOnlyRef ())
         {
            return SxVecRef<T, Layout>(std::move (*this)) -= x;
         } else {
            return SxVecCompute<algoScalar(Layout)>::subtract (*this, x);
         }
      }

      inline typename SxVecResult<typename SxTypeMapper<T>::TCmplx, Layout>::VecType
      operator- (const typename SxTypeMapper<T>::TCmplx &x) &&
      {
         SX_CHECK (allocMem);
         if (isContiguous(Layout) && isOnlyRef ())
         {
            return SxVecRef<T, Layout>(std::move (*this)) -= x;
         } else {
            return SxVecCompute<algoScalar(Layout)>::subtract (*this, x);
         }
      }

      inline SxVecRef<T, Layout>&
      operator-= (const T &x)
      {
         SxVecCompute<algoScalar(Layout)>::subtractInPlace (*this, x);
         return *this;
      }

      inline typename SxVecResult<decltype(T(0) * typename SxTypeMapper<T>::TReal(0)), Layout>::VecType
      operator* (typename SxTypeMapper<T>::TReal x) const&
      {
         return SxVecCompute<algoScalar(Layout)>::multiply (*this, x);
      }

      inline typename SxVecResult<typename SxTypeMapper<T>::TCmplx, Layout>::VecType
      operator* (const typename SxTypeMapper<T>::TCmplx &x) const&
      {
         return SxVecCompute<algoScalar(Layout)>::multiply (*this, x);
      }

      inline typename SxVecResult<decltype(T(0) * typename SxTypeMapper<T>::TReal(0)), Layout>::VecType
      operator* (typename SxTypeMapper<T>::TReal x) &&
      {
         SX_CHECK (allocMem);
         if (isContiguous(Layout) && isOnlyRef ())
         {
            return SxVecRef<T, Layout>(std::move (*this)) *= x;
         } else {
            return SxVecCompute<algoScalar(Layout)>::multiply (*this, x);
         }
      }

      inline typename SxVecResult<typename SxTypeMapper<T>::TCmplx, Layout>::VecType
      operator* (const typename SxTypeMapper<T>::TCmplx &x) &&
      {
         SX_CHECK (allocMem);
         if (isContiguous(Layout) && isOnlyRef ())
         {
            return SxVecRef<T, Layout>(std::move (*this)) *= x;
         } else {
            return SxVecCompute<algoScalar(Layout)>::multiply (*this, x);
         }
      }

      inline SxVecRef<T, Layout>&
      operator*= (const T &x)
      {
         SxVecCompute<algoScalar(Layout)>::multiplyInPlace (*this, x);
         return *this;
      }

      inline typename SxVecResult<decltype(T(0) / typename SxTypeMapper<T>::TReal(0)), Layout>::VecType
      operator/ (typename SxTypeMapper<T>::TReal x) const&
      {
         return SxVecCompute<algoScalar(Layout)>::divide (*this, x);
      }

      inline typename SxVecResult<typename SxTypeMapper<T>::TCmplx, Layout>::VecType
      operator/ (const typename SxTypeMapper<T>::TCmplx &x) const&
      {
         return SxVecCompute<algoScalar(Layout)>::divide (*this, x);
      }

      inline typename SxVecResult<decltype(T(0) / typename SxTypeMapper<T>::TReal(0)), Layout>::VecType
      operator/ (typename SxTypeMapper<T>::TReal x) &&
      {
         SX_CHECK (allocMem);
         if (isContiguous(Layout) && isOnlyRef ())
         {
            return SxVecRef<T, Layout>(std::move (*this)) /= x;
         } else {
            return SxVecCompute<algoScalar(Layout)>::divide (*this, x);
         }
      }

      inline typename SxVecResult<typename SxTypeMapper<T>::TCmplx, Layout>::VecType
      operator/ (const typename SxTypeMapper<T>::TCmplx &x) &&
      {
         SX_CHECK (allocMem);
         if (isContiguous(Layout) && isOnlyRef ())
         {
            return SxVecRef<T, Layout>(std::move (*this)) /= x;
         } else {
            return SxVecCompute<algoScalar(Layout)>::divide (*this, x);
         }
      }

      inline SxVecRef<T, Layout>&
      operator/= (const T &x)
      {
         SxVecCompute<algoScalar(Layout)>::divideInPlace (*this, x);
         return *this;
      }
      // --- scalar_op

      /// Unary minus
      inline typename SxVecResult<T,Layout>::VecType operator- () const&
      {
         return SxVecCompute<algoScalar(Layout)>::minus (*this);
      }

      /// Unary minus (reuse memory, if feasible)
      inline typename SxVecResult<T,Layout>::VecType minus () &&
      {
         return SxVecCompute<algoScalar(Layout)>::minusInPlace (std::move(*this));
      }

      // --- The next 48 lines were generated from snippets/SxVector.h snippet vec_func
      /// Square
      inline typename SxVecResult<T,Layout>::VecType sqr () const&
      {
         return SxVecCompute<algoScalar(Layout)>::sqr (*this);
      }

      /// Square (reuse memory, if feasible)
      inline typename SxVecResult<T,Layout>::VecType sqr () &&
      {
         return SxVecCompute<algoScalar(Layout)>::sqrInPlace (std::move(*this));
      }

      /// Cube (x*x*x)
      inline typename SxVecResult<T,Layout>::VecType cub () const&
      {
         return SxVecCompute<algoScalar(Layout)>::cub (*this);
      }

      /// Cube (x*x*x) (reuse memory, if feasible)
      inline typename SxVecResult<T,Layout>::VecType cub () &&
      {
         return SxVecCompute<algoScalar(Layout)>::cubInPlace (std::move(*this));
      }

      /// Absolute value
      inline typename SxVecResult<typename SxTypeMapper<T>::TReal,Layout>::VecType abs () const&
      {
         return SxVecCompute<algoScalar(Layout)>::abs (*this);
      }

      /// Absolute value (reuse memory, if feasible)
      inline typename SxVecResult<typename SxTypeMapper<T>::TReal,Layout>::VecType abs () &&
      {
         return SxVecCompute<algoScalar(Layout)>::absInPlace (std::move(*this));
      }

      /// Absolute value squared
      inline typename SxVecResult<typename SxTypeMapper<T>::TReal,Layout>::VecType absSqr () const&
      {
         return SxVecCompute<algoScalar(Layout)>::absSqr (*this);
      }

      /// Absolute value squared (reuse memory, if feasible)
      inline typename SxVecResult<typename SxTypeMapper<T>::TReal,Layout>::VecType absSqr () &&
      {
         return SxVecCompute<algoScalar(Layout)>::absSqrInPlace (std::move(*this));
      }
      // --- vec_func

      /// --- contractions
      /// Sum of elements
      inline T sum () const
      {
         return SxVecCompute<algoSelect(true,Layout,Layout)>::sum (*this);
      }

      /// Product of elements
      inline T product () const
      {
         return SxVecCompute<algoSelect(true,Layout,Layout)>::product (*this);
      }

      /// Sum absolute square of elements
      inline TReal normSqr () const
      {
         return SxVecCompute<algoSelect(true,Layout,Layout)>::normSqr (*this);
      }

      /// Norm of vector
      inline TReal norm () const
      {
         return SxVecCompute<algoSelect(true,Layout,Layout)>::norm (*this);
      }

      // --- conjugate and transpose

      /// Conjugate (for complex numbers)
      //inline const SxVecRef<Conjugate<T>, Layout> conj () const;
      inline typename SxVecResult<T, Layout>::VecType conj () const&
      {
         return SxVecCompute<algoScalar(Layout)>::conj (*this);
      }

      /// Conjugate (for complex numbers)
      inline typename SxVecResult<T, Layout>::VecType conj () &&
      {
         return SxVecCompute<algoScalar(Layout)>::conjInPlace (std::move(*this));
      }

      /// Matrix transpose
      inline typename SxVecResult<T, Layout>::VecType transpose () const
      {
         SX_CHECK (this->getSize () > 0, this->getSize ());
         return SxVecCompute<algoTrans(Layout)>::transpose (*this);
      }

      /// Matrix transpose and complex conjugate
      inline typename SxVecResult<T, Layout>::VecType adjoint () const
      {
         SX_CHECK (this->getSize () > 0, this->getSize ());
         return SxVecCompute<algoTrans(Layout)>::adjoint (*this);
      }

      // --- other
      /// Normalize
      inline void normalize ()
      {
         SX_CHECK (getNCols () == 1, getNCols ());
         (*this) /= norm ();
      }

      /** \brief Randomize elements and normalize each column
                 (unless it's a packed matrix) */
      void randomize ();

      // --- matrix functions
      /** \brief In-place rotation
          \param rotMat a square rotation matrix, to be applied from the right.
          This is a memory-friendly in-place rotation
          \f[
          M \leftarrow M U
          \f]
          where (this) M is a N x M matrix, and U (rotMat) is a M x M
          square matrix.
          \author C. Freysoldt
        */
      void rotate(const SxVecRef<T,Compact> &rotMat);

      /** \brief Overlap computation

        This computes the overlap matrix x^T y summing over the first nSum
        rows. It is equivalent to
        \code
        SxVecRef<T,SubMatrix> xSub = x.blockRef(0,0,nSum,x.getNCols ());
        SxVecRef<T,SubMatrix> ySub = y.blockRef(0,0,nSum,y.getNCols ());
        return xSub.adjoint () ^ ySub;
        \endcode
        but avoids the explicit submatrix objects and the adjoint.

        \note This is implemented only for Compact and SubMatrix layout.
        For other layouts, use an explicit adjoint.

      */
      template<StrideType Layout2>
      SxVector<T> overlap (const SxVecRef<T, Layout2> &y, ssize_t nSum) const;

      /** \brief Overlap computation

        This computes the overlap matrix x^T y. It is equivalent to
        \code
        x.adjoint () ^ y
        \endcode
        but avoids the explicit adjoint.

        \note This is implemented only for Compact and SubMatrix layout.
        For other layouts, use an explicit adjoint.

      */
      template<StrideType Layout2>
      inline SxVector<T> overlap (const SxVecRef<T, Layout2> &y) const
      {
         SX_CHECK (this->getNRows () == y.getNRows (),
                   this->getNRows (), y.getNRows ());
         return overlap (y, this->getNRows ());
      }

      /// Get reordered data
      SxVector<T> getSorted (const SxArray<ssize_t> &sortIdx) const;
};

#ifndef NDEBUG
template<class T, StrideType Layout>
inline void SX_VALIDATE_VECTOR (const SxVecRef<T,Layout> &in)
{
   SxVecConstIt<T,Layout> it = in.begin ();
   ssize_t n = in.getSize ();
   for (ssize_t i = 0; i < n; ++i, ++it)
      SX_CHECK_NUM (*it);
}
#endif /* NDEBUG */


/** \brief Vector class

    \author Christoph Freysoldt, freysoldt@mpie.de
    @note This vector class is a rewrite of the original Dirac vector class by
    \author Sixten Boeck
    */
template<class T>
class SX_EXPORT_MATH SxVector : public SxVecRef<T, Compact>
{
   public:
      // --- Constructors
      /// empty constructor
      inline SxVector () {}

      /// Copy constructor
      inline SxVector (const SxVector<T> &in);

      /// Copy constructor
      inline SxVector (const SxVecRef<T, Compact> &in)
      {
         if (in.allocMem) copy (in);
      }

      /// Constructor: copy from a reference
      template<class T2, StrideType Layout>
      SxVector (const SxVecRef<T2, Layout> &in)
      {
         SxVecRef<T,Compact>::operator= (in);
      }

      /// Move constructor
      inline SxVector (SxVector<T> &&in)
         : SxVecRef<T,Compact> (static_cast<SxVecRef<T,Compact>&&>(in))
      {
         // empty
      }

      /// Constructor: take over from rvalue reference
      inline SxVector (SxVecRef<T, Compact> &&in)
         : SxVecRef<T,Compact> (std::move (in))
      {
         // SxVector shouldn't be aliased
         SX_CHECK (!allocMem || allocMem->refCounter == 1,
                   allocMem ? int(allocMem->refCounter) : 0);
      }

      /// Constructor from array
      SxVector (const SxArray<T> &in)
         : SxVecRef<T, Compact> (in.getSize ())
      {
         if (size == 0) return;
         alloc (size);
         set (in);
      }

      /// Constructor from array (with typecast)
      template<class T2>
      explicit SxVector (const SxArray<T2> &in);

      /// Constructor from list
      SxVector (const SxList<T> &in);

      /// Constructor from stack
      inline SxVector (const SxStack<T> &in)
         : SxVecRef<T, Compact> (in.getSize ())
      {
         if (size == 0) return;
         alloc (size);
         in.exportStack (this->elements, size);
         SX_VALIDATE_VECTOR (*this);
      }

      /// Constructor from SxVector3
      inline SxVector (const SxVector3<T> &in)
         : SxVecRef<T, Compact> (3)
      {
         alloc (3);
         elements[0] = in(0);
         elements[1] = in(1);
         elements[2] = in(2);
      }

      /// Constructor from SxMatrix3
      inline SxVector (const SxMatrix3<T> &in)
         : SxVecRef<T, Compact> (9)
      {
         alloc (9);
         this->reshape (3,3);
         elements[0] = in(0,0);
         elements[1] = in(1,0);
         elements[2] = in(2,0);
         elements[3] = in(0,1);
         elements[4] = in(1,1);
         elements[5] = in(2,1);
         elements[6] = in(0,2);
         elements[7] = in(1,2);
         elements[8] = in(2,2);
      }

      /// Explicitly sized vector
      inline explicit SxVector (ssize_t n);

      /// Explicitly sized vector
      inline explicit SxVector (ssize_t nRows, ssize_t nCols);

      /// Create from basis
      inline explicit SxVector (const SxBasis &b)
         : SxVector(getNElements (b))
      {
         this->auxData.basisPtr = &b;
      }

      /// Create from basis and initialize (usually zero)
      inline explicit SxVector (const T val, const SxBasis &b)
         : SxVector(getNElements (b))
      {
         this->auxData.basisPtr = &b;
         this->set (val);
      }

   protected:
      template <SxVecAlgorithm> friend class SxVecCompute;
      template <class T2, StrideType Layout> friend class SxVecRef;
      /// Constructor for result from binary vector-vector operations
      template<StrideType Layout2>
      SxVector (SxAllocMem* mem, const SxVecLayout<Layout2> &layout,
                const SxAuxData &aux1, const SxAuxData &aux2)
      : SxVecRef<T,Compact>(mem, layout, aux1, aux2)
      {
         // empty
      }
      /// Constructor for result from matrix operations
      SxVector (SxAllocMem* mem, ssize_t rows, ssize_t cols,
                const SxAuxData &auxIn)
      : SxVecRef<T,Compact>(mem, mem ? (T*)mem->mem : NULL,
                            rows, 1, cols, rows, auxIn)
      {
         SX_VALIDATE_VECTOR (*this);
      }
   public:
      // --- change size
      /// Change size (loses data)
      inline void resize (ssize_t n);

      /** \brief Change size (keep and/or initialize data)
          @param keep    whether current values should be copied
          @param fillVal how to set new values
        */
      inline void resize (ssize_t n, bool keep,
                          const T& fillVal = T(0.));

      /// Change size (loses data)
      inline void reformat (ssize_t nRowsIn, ssize_t nColsIn);

      /// Reformat and set auxData according from existing reference
      template<class T2, StrideType Layout>
      inline void reformat (const SxVecRef<T2, Layout> &in)
      {
         reformat (in.getNRows (), in.getNCols ());
         this->auxData = in.auxData;
      }

      // --- assignment
      /// Copy data
      inline void copy (const SxVecRef<T,Compact> &);
      /// Copy data
      template <class T2, StrideType Layout>
      inline void copy (const SxVecRef<T2,Layout> &);

      /// Assignment (copies)
      inline SxVector<T>& operator= (const SxVecRef<T,Compact> &in)
      {
         if (in.allocMem) {
            copy (in);
         } else {
            this->unref ();
         }
         return *this;
      }

      /// Assignment (copies)
      inline SxVector<T>& operator= (const SxVector<T> &in)
      {
         return operator=(static_cast<const SxVecRef<T,Compact>&>(in));
      }

      /// Set from a vecref (copies if necessary)
      template <class T2, StrideType Layout>
      inline SxVector<T>& operator= (const SxVecRef<T2,Layout> &in)
      {
         if (this->getSize () == in.getSize () ) {
            SX_CHECK (!isPacked (Layout));
            // size is same, overwrite layout and auxData
            this->reshape (in.getNRows (), in.getNCols ());
            this->auxData = in.auxData;
         } else {
            if (this->allocMem) this->unref ();
         }
         SxVecRef<T,Compact>::operator= (in);
         return *this;
      }

      /// Move assignment: take over from rvalue reference
      inline SxVector<T>& operator= (SxVecRef<T, Compact> &&);

      /// Move assignment: take over from rvalue reference
      inline SxVector<T>& operator= (SxVector<T> &&in)
      {
         return operator= (static_cast<SxVecRef<T,Compact>&&>(in));
      }

      /** \brief Set values from stack
          @param stack the stack
          @param stackSize size of the stack
          @param offset offset into the vector

          \note The vector size must be large enough to receive all values.
       */
      inline void set (const SxStack<T> &stack, size_t stackSize, ssize_t offset = 0);

      // make visible other set routines
      using SxVecRef<T,Compact>::set;

      // --- getting references
      using SxVecRef<T,Compact>::operator ();
      /// Get a contiguous reference to part of vector
      inline SxVecRef<T,Compact> operator() (const SxIdx &idx)
      {
         SX_CHECK (idx.start >= 0, idx.start);
         SX_CHECK (idx.end < this->getSize (), idx.end, this->getSize ());
         return SxVecRef<T, Compact> (allocMem, elements + idx.start,
                                      idx.end - idx.start + 1);
      }
      /// Get a contiguous reference to part of vector
      inline const SxVecRef<T,Compact> operator() (const SxIdx &idx) const
      {
         SX_CHECK (idx.start >= 0, idx.start);
         SX_CHECK (idx.end < this->getSize (), idx.end, this->getSize ());
         return SxVecRef<T, Compact> (allocMem, elements + idx.start,
                                      idx.end - idx.start + 1);
      }
      /// Get a submatrix
      inline SxVecRef<T,SubMatrix>
      blockRef (ssize_t rowOffset, ssize_t colOffset,
                ssize_t nRowsIn, ssize_t nColsIn)
      {
         SX_CHECK (nRowsIn > 0, nRowsIn);
         SX_CHECK (nColsIn > 0, nColsIn);
         SX_CHECK(rowOffset + nRowsIn <= nRows, rowOffset, nRowsIn, nRows);
         SX_CHECK(colOffset + nColsIn <= nCols, colOffset, nColsIn, nCols);
         return SxVecRef<T, SubMatrix>
            (allocMem, elements + rowOffset + colOffset * nRows,
             nRowsIn, 1, nColsIn, nRows);
      }

      /// Get a submatrix (const)
      inline const SxVecRef<T,SubMatrix>
      blockRef (ssize_t rowOffset, ssize_t colOffset,
                ssize_t nRowsIn, ssize_t nColsIn) const
      {
         SX_CHECK (nRowsIn > 0, nRowsIn);
         SX_CHECK (nColsIn > 0, nColsIn);
         SX_CHECK(rowOffset + nRowsIn <= nRows, rowOffset, nRowsIn, nRows);
         SX_CHECK(colOffset + nColsIn <= nCols, colOffset, nColsIn, nCols);
         return SxVecRef<T, SubMatrix>
            (allocMem, elements + rowOffset + colOffset * nRows,
             nRowsIn, 1, nColsIn, nRows);
      }

      // --- functions that make sense on 1D vectors
      /// Find minimum value
      T minval (ssize_t *idx=NULL) const;
      /// Find maximum value
      T maxval (ssize_t *idx=NULL) const;
      /// Find minimum value
      inline T minval (int *idx) const
      {
         ssize_t i;
         T res = minval (&i);
         SX_CHECK (int(i) == i);
         *idx = int(i);
         return res;
      }
      /// Find maximum value
      inline T maxval (int *idx) const
      {
         ssize_t i;
         T res = maxval (&i);
         SX_CHECK (int(i) == i);
         *idx = int(i);
         return res;
      }

      /// Get sorting index
      inline SxArray<ssize_t> getSortIdx () const;

      /// 1D Integrate (Simpson rule)
      T integrate (double step) const;

      using SxVecRef<T,Compact>::sum;

      /// Sum of elements
      inline T sum (ssize_t startIdx, ssize_t endIdx=-1) const
      {
         SX_CHECK (elements);
         ssize_t nSum = (endIdx<0) ? (size-startIdx) : (endIdx+1-startIdx);
         return SxVecCompute<Simd>::sum (elements + startIdx, nSum);
      }

      // --- matrix functions
      /// Invert a square matrix
      inline SxVector<T> inverse () const&;
      /// Invert a square matrix
      inline SxVector<T> inverse () &&;

      /// pseudo-inverse of a arbitrary matrix
      inline SxVector<T>
      pseudoInverse (typename SxTypeMapper<T>::TReal threshold = 1e-10) const&;

      /** \brief Solve A ^ x = b (in the linear least square sense)

          This solves |A ^ x - b|^2 = min. via singular value decomposition
        */
      SxVector<T> solve (const SxVector<T>  &b) const&;
      SxVector<T> solve (const SxVector<T>  &b) &&;
      SxVector<T> solve (      SxVector<T> &&b) const&;
      SxVector<T> solve (      SxVector<T> &&b) &&;

      /// Cholesky decomposition
      SxVector<T> choleskyDecomposition (enum UPLO = LowerLeft) const;

      /// Check whether matrix is Hermitian
      inline bool isHermitian () const
      {
         return SxVecCompute<Simd>::isHermitian (*this);
      }

      // --- printing
      /** Print out a vector either in simple output format or in vector
          format (similar to \URL[Mathematica]{http://www.wolfram.com}) */
      inline void print (bool vectorForm=false) const;
      void print (std::ostream &, bool vectorForm=false) const;
      /** Recusive print */
      void  println (std::ostream &, bool) const;

   protected:
      // --- make some SxVecRef stuff visible
      using SxVecRef<T,Compact>::allocMem;
      using SxVecRef<T,Compact>::nRows;
      using SxVecRef<T,Compact>::nCols;
      using SxVecRef<T,Compact>::size;
      using SxVecRef<T,Compact>::alloc;
   public:
      using SxVecRef<T,Compact>::isOnlyRef;
      using SxVecRef<T,Compact>::elements;
      using SxVecRef<T,Compact>::getSize;

};

template<class T, StrideType Layout>
SxVecRef<T,Layout>::operator SxMatrix3<T> ()
{
   SX_CHECK (getNCols () == 3, getNCols ());
   SX_CHECK (getNRows () == 3, getNRows ());
   SX_CHECK (!isPacked(Layout));
   return SxMatrix3<T> ((*this)(0,0), (*this)(0,1), (*this)(0,2),
                        (*this)(1,0), (*this)(1,1), (*this)(1,2),
                        (*this)(2,0), (*this)(2,1), (*this)(2,2));
}

template<class T, StrideType Layout>
SxVecRef<T,Layout>::SxVecRef (SxVecRef<T,Layout> &in)
   : SxVecLayout<Layout> (in),
     allocMem (in.allocMem),
     elements (in.elements),
     auxData (in.auxData)
{
   allocMem->ref ();
}

template<class T, StrideType Layout>
SxVecRef<T,Layout>::SxVecRef (SxVecRef<T,Layout> &&in)
   : SxVecLayout<Layout> (in),
     allocMem (in.allocMem),
     elements (in.elements),
     auxData (in.auxData)
{
   in.allocMem = NULL;
   in.elements = NULL;
   // in.unset ();
   // cout << "MOVE" << endl;
}

/// Copy constructor from different type (copies)
template<class T, StrideType Layout>
template<class T2>
SxVecRef<T,Layout>::SxVecRef (const SxVecRef<T2, Layout> &in)
   : auxData (in.auxData)
{
   if (!in.allocMem)  {
      allocMem = NULL;
      elements = NULL;
      return;
   }
   alloc (in.getSize ());
   /* constexpr */ if (isPacked(Layout))  {
      copy (this, in); // copy layout
   } else {
      SxVecLayout<Layout>::set (in.getNRows (), in.getNCols (),
                                1, in.getNRows ());
   }
   SX_CHECK (getSize () == in.getSize (), getSize (), in.getSize ());
   // ignore shape
   SxVecIt<T,Layout> it = begin ();
   SxVecConstIt<T2,Layout> inIt = in.begin ();
   for ( ; it != end () ; ++it, ++inIt)  {
      SX_CHECK_NUM (*inIt);
      *it = static_cast<T>(*inIt);
   }
}

// Layout cast constructor -> references data
template<class T, StrideType Layout>
template<StrideType Layout2>
SxVecRef<T,Layout>::SxVecRef (SxVecRef<T,Layout2> &in)
   : SxVecLayout<Layout> (in.getNRows (), in.getNCols (),
                          in.getRowStride (), in.getColStride ()),
     allocMem (in.allocMem),
     elements (in.elements),
     auxData (in.auxData)
{
   SX_CHECK (!isPacked(Layout));
   SX_CHECK (!isPacked(Layout2));
   // check that Layout-cast doesn't affect strides
   // if this fails, you are forced to copy the data
   // via newref (SxVector<T> (ref))
   SX_CHECK (getRowStride () == in.getRowStride (),
             getRowStride (), in.getRowStride ());
   SX_CHECK (getColStride () == in.getColStride (),
             getColStride (), in.getColStride ());
   allocMem->ref ();
}

template<class T, StrideType Layout>
SxVecRef<T,Layout>::SxVecRef (SxAllocMem* allocIn,
                              T*          data,
                              ssize_t     nRowsIn,
                              ssize_t     rowStrideIn,
                              ssize_t     nColsIn,
                              ssize_t     colStrideIn)
   : SxVecLayout<Layout> (nRowsIn, nColsIn, rowStrideIn, colStrideIn),
     allocMem(allocIn), elements(data)
{
   SX_CHECK(allocIn);
   allocIn->ref (); // register this reference
#ifndef NDEBUG
   if (isPacked(Layout))  {
      SX_CHECK(allocMem->isInside (elements, this->size),
               *allocMem, elements, this->size);
   } else {
      ssize_t rStride = getRowStride (), cStride = getColStride ();
      SX_CHECK(allocMem->isInside (elements,  (nCols-1) * cStride
                                        + (nRows-1)*rStride + 1),
               *allocMem, elements, (nCols-1) * cStride, (nRows-1)*rStride+1);
   }
#endif
}

template<class T, StrideType Layout>
SxVecRef<T,Layout>::SxVecRef (SxAllocMem* allocIn,
                              T*          data,
                              ssize_t     nRowsIn,
                              ssize_t     rowStrideIn,
                              ssize_t     nColsIn,
                              ssize_t     colStrideIn,
                              const SxAuxData &auxIn)
   : SxVecLayout<Layout> (nRowsIn, nColsIn, rowStrideIn, colStrideIn),
     allocMem(allocIn), elements(data), auxData (auxIn)
{
   SX_CHECK(allocIn);
   allocIn->ref (); // register this reference
#ifndef NDEBUG
   if (isPacked(Layout))  {
      SX_CHECK(allocMem->isInside (elements, this->size),
               *allocMem, elements, this->size);
   } else {
      ssize_t rStride = getRowStride (), cStride = getColStride ();
      SX_CHECK(allocMem->isInside (elements,  (nCols-1) * cStride
                                        + (nRows-1)*rStride + 1),
               *allocMem, elements, (nCols-1) * cStride, (nRows-1)*rStride+1);
   }
#endif
}

/// Copy assignment (associates if unassociated, otherwise sets)
template<class T, StrideType Layout>
inline SxVecRef<T,Layout>& SxVecRef<T,Layout>::operator= (const SxVecRef<T,Layout>& in)
{
   if (&in == this) return *this;
   if (!in.allocMem)  {
      unref ();
      return *this;
   }
   if (!allocMem)  {
      SxVecLayout<Layout>::operator= (in);
      allocMem = in.allocMem;
      allocMem->ref ();
      elements = in.elements;
      auxData = in.auxData;
   } else {
      operator<<= (in);
   }
   return *this;
}

/// Copy assignment (sets)
template<class T, StrideType Layout>
template<StrideType Layout2>
inline SxVecRef<T,Layout>&
SxVecRef<T,Layout>::operator<<= (const SxVecRef<T,Layout2>& in)
{
   SX_CHECK ((void*)this != (const void*)&in);
   SX_CHECK (getSize () == in.getSize (), getSize (), in.getSize ());
   // ignore shape
   SxVecIt<T,Layout> it = begin ();
   SxVecConstIt<T,Layout2> inIt = in.begin ();
   for ( ; it != end () ; ++it, ++inIt)  {
      SX_CHECK_NUM (*inIt);
      *it = *inIt;
   }
   return *this;
}

/// Copy assignment (sets)
template<class T, StrideType Layout>
inline SxVecRef<T,Layout>&
SxVecRef<T,Layout>::operator<<= (const SxArray<T>& in)
{
   SX_CHECK (getSize () == in.getSize (), getSize (), in.getSize ());
   // ignore shape
   SxVecIt<T,Layout> it = begin ();
   typename SxArray<T>::ConstIterator inIt = in.begin ();
   for ( ; it != end () ; ++it, ++inIt)  {
      SX_CHECK_NUM (*inIt);
      *it = *inIt;
   }
   return *this;
}

/// Assignment from different SxVecRef (copies)
template<class T, StrideType Layout>
template<class T2, StrideType Layout2>
SxVecRef<T,Layout>&
SxVecRef<T,Layout>::operator= (const SxVecRef<T2, Layout2> &in)
{
   if (!in.allocMem)  {
      unref ();
      return *this;
   }
   if (!allocMem)  {
      alloc (in.getSize ());
      /* constexpr */ if (isPacked(Layout))  {
         SX_CHECK (Layout == Layout2);
         copy (this, in); // copy layout
      } else {
         SX_CHECK (!isPacked(Layout2));
         SxVecLayout<Layout>::set (in.getNRows (), in.getNCols (),
                                   1, in.getNRows ());
      }
      auxData = in.auxData;
   }
   SX_CHECK (getSize () == in.getSize (), getSize (), in.getSize ());
   // ignore shape
   SxVecIt<T,Layout> it = begin ();
   SxVecConstIt<T2,Layout2> inIt = in.begin ();
   for ( ; it != end () ; ++it, ++inIt)  {
      SX_CHECK_NUM (*inIt);
      *it = static_cast<T>(*inIt);
   }
   return *this;
}

template<class T, StrideType Layout>
template<StrideType Layout2>
SxVecRef<T,Layout>& SxVecRef<T,Layout>::operator= (SxVecRef<T, Layout2> &&in)
{
   SX_CHECK ((void*)this != (void*)&in); // X = std::move(X) is insane
   if (!in.allocMem)  {
      unref ();
      return *this;
   }
   SX_CHECK (!isPacked(Layout)); // to be implemented
   SX_CHECK (!isPacked(Layout2)); // to be implemented
   if (isOnlyRef ())  {
      bool row1  = in.getRowStride () == 1;
      bool colOk = (in.getNCols () == 1) ||
                   (in.getColStride () == in.getNRows () * in.getRowStride ());

      if (   (Layout == Compact && row1 && colOk)
          || (Layout == Strided && colOk)
          || (Layout == SubMatrix && row1)
          || (Layout == GeneralMatrix))
      {
         if (allocMem) unref_ ();
         SxVecLayout<Layout>::set (in.getNRows (), in.getNCols (),
                                   in.getRowStride (), in.getColStride ());
         allocMem = in.allocMem;
         elements = in.elements;
         auxData = in.auxData;
         in.allocMem = NULL;
         in.elements = NULL;
         in.unset ();
         return *this;
      }
   }
   return operator= (static_cast<const SxVecRef<T, Layout2>&>(in));
}

template<class T, StrideType Layout>
void SxVecRef<T,Layout>::set (const T& in)
{
   for (SxVecIt<T,Layout> it = begin (); it != end () ; ++it)
      *it = in;
}

template<class T, StrideType Layout>
void SxVecRef<T,Layout>::set (const SxArray<T>& in)
{
   SX_CHECK (in.getSize () == getSize (), in.getSize (), getSize ());
   typename SxArray<T>::ConstIterator inIt = in.begin ();
   for (SxVecIt<T,Layout> it = begin (); it != end () ; ++it)  {
      SX_CHECK_NUM (*inIt);
      *it = *inIt++;
   }
}

// ---

template<class T, StrideType Layout>
inline SxVecRef<T,colStrideType(Layout)> SxVecRef<T,Layout>::colRef (ssize_t c)
{
   SX_CHECK (!isPacked(Layout));
   SX_CHECK (c >= 0 && c < getNCols (), c, getNCols ());
   return SxVecRef<T, colStrideType(Layout)>
      (allocMem, elements + c * getColStride (),
       nRows, getRowStride (), 1, getColStride (), auxData);
}

template<class T, StrideType Layout>
inline const SxVecRef<T, colStrideType(Layout)>
SxVecRef<T,Layout>::colRef (ssize_t c) const
{
   SX_CHECK (!isPacked(Layout));
   SX_CHECK (c >= 0 && c < getNCols (), c, getNCols ());
   return SxVecRef<T, colStrideType(Layout)>
      (allocMem, elements + c * getColStride (),
       nRows, getRowStride (), 1, getColStride (), auxData);
}

template<class T, StrideType Layout>
inline SxVecRef<T, Strided> SxVecRef<T, Layout>::rowRef (ssize_t r)
{
   SX_CHECK (!isPacked(Layout));
   SX_CHECK (r >= 0 && r < getNRows (), r, getNRows ());
   return SxVecRef<T, Strided> (allocMem, elements + r * getRowStride (), nCols, getColStride ());
}

template<class T, StrideType Layout>
inline const SxVecRef<T, Strided> SxVecRef<T, Layout>::rowRef (ssize_t r) const
{
   SX_CHECK (!isPacked(Layout));
   SX_CHECK (r >= 0 && r < getNRows (), r, getNRows ());
   return SxVecRef<T, Strided> (allocMem, elements + r * getRowStride (), nCols, getColStride ());
}


// Reference to diagonal (of square matrix)
template<class T, StrideType Layout>
inline SxVecRef<T, Strided> SxVecRef<T, Layout>::diag ()
{
   SX_CHECK (nRows == nCols, nRows, nCols);
   return SxVecRef<T, Strided>(allocMem, elements, nRows,
                               getColStride () + getRowStride ());
}

// Reference to diagonal (of square matrix) -- const
template<class T, StrideType Layout>
inline const SxVecRef<T, Strided> SxVecRef<T, Layout>::diag () const
{
   SX_CHECK (nRows == nCols, nRows, nCols);
   return SxVecRef<T, Strided>(allocMem, elements, nRows,
                               getColStride () + getRowStride ());
}

// Reference to single component
template<class T, StrideType Layout>
inline SxVecRef<T,subMatrixLayout(Layout)>
SxVecRef<T,Layout>::compRef (int iComp)
{
   SX_CHECK (!isPacked(Layout));
   SX_CHECK (iComp >= 0 && iComp < auxData.nComp, iComp, (int)auxData.nComp);
   SX_CHECK (nRows % auxData.nComp == 0, nRows, (int)auxData.nComp);
   ssize_t nElem = nRows / ssize_t(auxData.nComp);
   SX_CHECK ((!auxData.basisPtr) || (getNElements (*auxData.basisPtr) == nElem),
             getNElements (*auxData.basisPtr), nElem, (int)auxData.nComp);

   SxVecRef<T, subMatrixLayout(Layout)> res
      (allocMem, elements + iComp * nElem * getRowStride (),
       nElem, getRowStride (), nCols, getColStride (), auxData);
   res.auxData.nComp = 1;
   return res;
}

// Reference to single component -- const
template<class T, StrideType Layout>
inline const SxVecRef<T,subMatrixLayout(Layout)>
SxVecRef<T,Layout>::compRef (int iComp) const
{
   return const_cast<SxVecRef<T,Layout>*>(this)->compRef (iComp);
}

template<>
class SxVectorReIm<true>
{
   public:
      /// --- Get real part (copy)
      template<class T, StrideType Layout>
      static inline SxVecRef<T, Layout> real (const SxVecRef<SxComplex<T>,Layout> &vec)
      {
         return SxVecCompute<Simd>::real (vec);
      }

      /// --- Get imaginary part (copy)
      template<class T, StrideType Layout>
      static inline SxVecRef<T, Layout> imag (const SxVecRef<SxComplex<T>,Layout> &vec)
      {
         return SxVecCompute<Simd>::imag (vec);
      }
} ;

template<>
class SxVectorReIm<false>
{
   public:
      /// --- Get real part
      template<class T, StrideType Layout>
      static inline SxVecRef<T,realPartLayout(Layout)>
      real (const SxVecRef<SxComplex<T>, Layout> &vec)
      {
         return SxVecRef<T,realPartLayout(Layout)>
            (vec.allocMem, (T*)vec.elements, vec.nRows, 2 * vec.getRowStride (),
             vec.nCols, 2 * vec.getColStride (), vec.auxData);
      }

      /// --- Get imaginary part
      template<class T, StrideType Layout>
      static inline SxVecRef<T,realPartLayout(Layout)>
      imag (const SxVecRef<SxComplex<T>,Layout> &vec)
      {
         return SxVecRef<T,realPartLayout(Layout)>
            (vec.allocMem, ((T*)vec.elements)+1, vec.nRows, 2 * vec.getRowStride (),
             vec.nCols, 2 * vec.getColStride (), vec.auxData);
      }
} ;

template <class T, StrideType Layout>
template<StrideType Layout2>
SxVecRef<T, Layout2>
SxVecRef<T,Layout>::getRef (ssize_t offset,
                            ssize_t nRowsIn, ssize_t nColsIn,
                            ssize_t rStride, ssize_t cStride) const
{
   SX_CHECK (allocMem);
   SX_CHECK (offset + nRowsIn * nColsIn <= getSize (),
             offset, nRowsIn*nColsIn, getSize ());
   SX_CHECK (rStride % getRowStride () == 0, rStride, getRowStride ());
   if (Layout2 == SubMatrix || Layout2 == GeneralMatrix)  {
      if (cStride < 0)  {
         cStride = getColStride ();
      }
      SX_CHECK (offset % getNRows () + (nRowsIn - 1) * rStride / getRowStride ()
                < getNRows ());
      SX_CHECK (cStride % getColStride () == 0, cStride, getColStride ());
   } // else: cStride is ignored anyway
   return SxVecRef<T, Layout2> (allocMem, elements + offset * getRowStride (),
                                nRowsIn, rStride, nColsIn, cStride);
}

// --- change shape
template<class T, StrideType Layout>
SxVecRef<T,Layout>&
SxVecRef<T,Layout>::reshape (ssize_t nRowsIn, ssize_t nColsIn)
{
   SX_CHECK (!isPacked(Layout));
   SX_CHECK (nRowsIn * nColsIn == nRows * nCols,
             nRows, nCols, nRowsIn, nColsIn);
   SX_CHECK (nRows * getRowStride () == getColStride (),
             nRows, getRowStride (), getColStride ());
   SxVecLayout<Layout>::set (nRowsIn, nColsIn, getRowStride (),
                             nRowsIn * getRowStride ());
   return *this;
}

namespace {
   class SxGetRandom {
      public:
         template<class T>
         inline operator SxComplex<T> () const
         {
            return SxComplex<T>(SxRandom::get (), SxRandom::get ());
         }
         template<class T>
         inline operator T () const
         {
            return T(SxRandom::get ());
         }
   };
}


template<class T, StrideType Layout>
void SxVecRef<T, Layout>::randomize ()
{
   for (Iterator it = begin (); it != end (); ++it)
      *it = SxGetRandom ();

   /* constexpr */
   if (!isPacked(Layout)) {
      for (ssize_t ic = 0; ic < this->getNCols (); ic++)
         colRef(ic).normalize ();
   }
}

// ---

template<class T, StrideType Layout>
std::ostream& operator<< (std::ostream& out, const SxVecRef<T,Layout> &vec)
{
   if (vec.getNCols () > 1 && vec.getNRows () > 1 && !isPacked(Layout))
   {
      out << '[';
      for (ssize_t ir = 0; ir < vec.getNRows (); ir++)  {
         out << '[' << vec(ir,0);
         for (ssize_t ic = 1; ic < vec.getNCols (); ic++)  {
            out << ',' << ' ' << vec(ir, ic);
         }
         out << ']';
         if (ir != vec.getNRows () - 1) out << ',';
      }
      out << ']';
   } else {
      typename SxVecRef<T,Layout>::ConstIterator it = vec.begin ();
      out << '[' << *it++;
      for ( ; it != vec.end (); ++it)
         out << ',' << ' ' << *it;
      out << ']';
   }
   return out;
}

template<class T>
void SxVector<T>::print (bool vectorForm) const
{
   print (cout, vectorForm);
   cout << endl;
}

template<class T>
void SxVector<T>::print (std::ostream &s, bool vectorForm) const
{
   int prec = static_cast<int>(s.precision());
   ios::fmtflags flags = s.flags ();
   if (vectorForm)  {
      s.setf (ios::fixed, ios::floatfield);
      s.precision (3);
      s.fill(' ');
      s << setw(7);
   }

   println (s, vectorForm);

   if (vectorForm)  { s << endl; }

   // --- restore predefined format
//   s.setf (ios_base::fmtflags(0), ios_base::floatfield);
   s.precision (prec);
   s.flags (flags);
}
template<class T>
void SxVector<T>::println (std::ostream &s, bool vectorForm) const
{
   ssize_t i, r, c;
   int width = static_cast<int>(s.width());

   if ( !vectorForm )  s << "[";
   if (getSize () == 0)  {            // empty

      s << "empty";

   } else if (this->getNCols () <= 1)  {     // vector print out

      for (i=0; i < this->getSize () - 1; i++)
         s << setw(width) << elements[i] << ", ";
      s << setw(width) << elements[this->getSize () - 1];

   } else {

      for (r=0; r < this->getNRows () ; r++)  {  // matrix print out
         if ( !vectorForm )  s << "[";
         for (c=0; c < (this->getNCols ())-1; c++)  {
            s << setw(width) << (*this)(r,c) << ", ";
         }
         s << setw(width) << (*this)(r,this->getNCols ()-1);
         if ( vectorForm) {s << endl;} else {s << "]";}
         if (!vectorForm && r < this->getNRows () -1)  s << ",";
      }
   }
   if ( !vectorForm )  s << "]";
}


// --- The next 116 lines were generated from snippets/SxVector.h snippet vec_op_dyingvec
// vector x + y (x, y are compact, y is about to die)
template<class T>
inline SxVector<T>
operator+ (const SxVecRef<T, Compact> &x, SxVecRef<T,Compact> &&y)
{
   if (y.isOnlyRef ())  {
      // cout << "MOVE(y) add" << endl;
      return SxVecCompute<Simd>::add (x, std::move(y));
   } else {
      return SxVecCompute<Simd>::add
         (x, static_cast<const SxVecRef<T,Compact>&>(y));
   }
}

// vector x + y (x, y are compact and both are about to die)
template<class T>
inline SxVector<T>
operator+ (SxVecRef<T, Compact> &&x, SxVecRef<T,Compact> &&y)
{
   if (y.isOnlyRef ())  {
      // cout << "MOVE(y) add" << endl;
      return SxVecCompute<Simd>::add
         (static_cast<const SxVecRef<T,Compact>&>(x), std::move(y));
   } else {
      return SxVecCompute<Simd>::add
         (std::move(x), static_cast<const SxVecRef<T,Compact>&>(y));
   }
}

// vector x - y (x, y are compact, y is about to die)
template<class T>
inline SxVector<T>
operator- (const SxVecRef<T, Compact> &x, SxVecRef<T,Compact> &&y)
{
   if (y.isOnlyRef ())  {
      // cout << "MOVE(y) subtract" << endl;
      return SxVecCompute<Simd>::subtract (x, std::move(y));
   } else {
      return SxVecCompute<Simd>::subtract
         (x, static_cast<const SxVecRef<T,Compact>&>(y));
   }
}

// vector x - y (x, y are compact and both are about to die)
template<class T>
inline SxVector<T>
operator- (SxVecRef<T, Compact> &&x, SxVecRef<T,Compact> &&y)
{
   if (y.isOnlyRef ())  {
      // cout << "MOVE(y) subtract" << endl;
      return SxVecCompute<Simd>::subtract
         (static_cast<const SxVecRef<T,Compact>&>(x), std::move(y));
   } else {
      return SxVecCompute<Simd>::subtract
         (std::move(x), static_cast<const SxVecRef<T,Compact>&>(y));
   }
}

// vector x * y (x, y are compact, y is about to die)
template<class T>
inline SxVector<T>
operator* (const SxVecRef<T, Compact> &x, SxVecRef<T,Compact> &&y)
{
   if (y.isOnlyRef ())  {
      // cout << "MOVE(y) multiply" << endl;
      return SxVecCompute<Simd>::multiply (x, std::move(y));
   } else {
      return SxVecCompute<Simd>::multiply
         (x, static_cast<const SxVecRef<T,Compact>&>(y));
   }
}

// vector x * y (x, y are compact and both are about to die)
template<class T>
inline SxVector<T>
operator* (SxVecRef<T, Compact> &&x, SxVecRef<T,Compact> &&y)
{
   if (y.isOnlyRef ())  {
      // cout << "MOVE(y) multiply" << endl;
      return SxVecCompute<Simd>::multiply
         (static_cast<const SxVecRef<T,Compact>&>(x), std::move(y));
   } else {
      return SxVecCompute<Simd>::multiply
         (std::move(x), static_cast<const SxVecRef<T,Compact>&>(y));
   }
}

// vector x / y (x, y are compact, y is about to die)
template<class T>
inline SxVector<T>
operator/ (const SxVecRef<T, Compact> &x, SxVecRef<T,Compact> &&y)
{
   if (y.isOnlyRef ())  {
      // cout << "MOVE(y) divide" << endl;
      return SxVecCompute<Simd>::divide (x, std::move(y));
   } else {
      return SxVecCompute<Simd>::divide
         (x, static_cast<const SxVecRef<T,Compact>&>(y));
   }
}

// vector x / y (x, y are compact and both are about to die)
template<class T>
inline SxVector<T>
operator/ (SxVecRef<T, Compact> &&x, SxVecRef<T,Compact> &&y)
{
   if (y.isOnlyRef ())  {
      // cout << "MOVE(y) divide" << endl;
      return SxVecCompute<Simd>::divide
         (static_cast<const SxVecRef<T,Compact>&>(x), std::move(y));
   } else {
      return SxVecCompute<Simd>::divide
         (std::move(x), static_cast<const SxVecRef<T,Compact>&>(y));
   }
}
// --- vec_op_dyingvec

/// --- scalar product (real values)
template<class T, StrideType Layout1, StrideType Layout2>
T dot(const SxVecRef<T,Layout1> &x, const SxVecRef<T,Layout2> &y)
{
   SX_CHECK(x.getSize () == y.getSize (), x.getSize (), y.getSize ());
   SX_CHECK(x.getNRows () == y.getNRows (), x.getNRows (), y.getNRows ());
   SX_CHECK(x.getNCols () == y.getNCols (), x.getNCols (), y.getNCols ());
   T res = 0;
   SxVecConstIt<T, Layout1> xIt = x.begin ();
   SxVecConstIt<T, Layout2> yIt = y.begin ();
   T res2 = 0, res3 = 0, res4 = 0;
   for (ssize_t n = x.getSize () / 4; n>0; --n)  {
      res += *xIt++ * *yIt++;
      res2 += *xIt++ * *yIt++;
      res3 += *xIt++ * *yIt++;
      res4 += *xIt++ * *yIt++;
   }
   res += res2 + res3 + res4;
   for ( ; xIt != x.end (); ++xIt, ++yIt)
      res += *xIt * *yIt;
   SX_CHECK_NUM(res);
   return res;
}

/// --- scalar product (complex values)
template<class T, StrideType Layout1, StrideType Layout2>
SxComplex<T> dot(const SxVecRef<SxComplex<T>,Layout1> &x,
                 const SxVecRef<SxComplex<T>,Layout2> &y)
{
   SX_CHECK(x.getSize () == y.getSize (), x.getSize (), y.getSize ());
   SX_CHECK(x.getNRows () == y.getNRows (), x.getNRows (), y.getNRows ());
   SX_CHECK(x.getNCols () == y.getNCols (), x.getNCols (), y.getNCols ());
   SxComplex<T> res(0.,0.);
   SxVecConstIt<SxComplex<T>, Layout1> xIt = x.begin ();
   SxVecConstIt<SxComplex<T>, Layout2> yIt = y.begin ();
   for ( ; xIt != x.end (); ++xIt, ++yIt)
      res += (*xIt).conj () * *yIt;
   SX_CHECK_NUM(res);
   return res;
}


template<class T, StrideType Layout>
void SxVecRef<T,Layout>::rotate(const SxVecRef<T,Compact> &rotMat)
{
   SX_CHECK (Layout == Compact);
   SX_CHECK (getSize () > 0);
   SX_CHECK (rotMat.getNRows () == rotMat.getNCols (),
             rotMat.getNRows (), rotMat.getNCols ());
   SX_CHECK (this->getNCols () == rotMat.getNCols (),
             this->getNCols (), rotMat.getNCols ());
   ::inPlaceRot (elements, rotMat.elements, (int)getNRows (), (int)getNCols ());
   SX_VALIDATE_VECTOR (*this);
}

template<class T, StrideType Layout>
template<StrideType Layout2>
SxVector<T> SxVecRef<T, Layout>::overlap (const SxVecRef<T, Layout2> &y, ssize_t nSum) const
{
   SX_CHECK (Layout  == Compact || Layout  == SubMatrix);
   SX_CHECK (Layout2 == Compact || Layout2 == SubMatrix);
   SX_CHECK (this->getNRows () >= nSum, this->getNRows (), nSum);
   SX_CHECK (y.getNRows () >= nSum, y.getNRows (), nSum);
   // check that int is enough to hold values
   SX_CHECK ((int)this->getNCols () == this->getNCols (), this->getNCols ());
   SX_CHECK ((int)y.getNCols () == y.getNCols (), y.getNCols ());
   SX_CHECK ((int)nSum == nSum, nSum);
   SX_CHECK ((int)this->getColStride () == this->getColStride (),
             this->getColStride ());
   SX_CHECK ((int)y.getColStride () == y.getColStride (), y.getColStride ());

   int nX = (int)this->getNCols (), nY = (int)y.getNCols ();
   SxVector<T> res(nX, nY);
   res.auxData.setCombined (auxData, y.auxData);
   res.auxData.basisPtr = NULL; // doesn't make sense after row contraction
   ::matovlp (res.elements, this->elements, y.elements,
              (int)this->getColStride (), nX, (int)y.getColStride (), nY, (int)nSum);
   SX_VALIDATE_VECTOR(res);
   return res;
}

/// Get reordered data
template<class T, StrideType Layout>
SxVector<T> SxVecRef<T, Layout>::getSorted (const SxArray<ssize_t> &sortIdx) const
{
   SX_CHECK (getSize () == sortIdx.getSize (),
             getSize (), sortIdx.getSize ());
   SX_CHECK (getSize () > 0, getSize ());
   SxAllocMem* newAlloc = SxAllocMem::create (sizeof(T) * getSize ());
   T* newElem = (T*)newAlloc->mem;
   for (ssize_t i = 0; i < getSize (); ++i)  {
      newElem[i] = elements[getIdx (sortIdx(i))];
      SX_CHECK_NUM(newElem[i]);
   }
   return SxVector<T>(newAlloc, *this, this->auxData, this->auxData);
}

// --- SxVector implementation
// --- Constructors
/// Constructor: copy from a reference
template<class T>
SxVector<T>::SxVector (const SxVector<T> &in)
   : SxVecRef<T,Compact> (static_cast<const SxVecLayout<Compact>&>(in),
                          in.auxData)
{
   SX_CHECK (getSize () > 0);
   memcpy (elements, in.elements, in.getSize () * sizeof(T));
}

// Constructor from array with type cast
template<class T>
template<class T2>
SxVector<T>::SxVector (const SxArray<T2> &in)
   : SxVecRef<T,Compact> (in.getSize ())
{
   if (in.getSize () == 0) return;
   alloc (in.getSize ());
   typename SxArray<T2>::ConstIterator inIt = in.begin ();
   for (SxVecIt<T,Compact> it = this->begin (); it != this->end () ; ++it)  {
      SX_CHECK_NUM (*inIt);
      *it = T(*inIt++);
   }
}


/// Constructor from list
template<class T>
SxVector<T>::SxVector (const SxList<T> &in)
   : SxVecRef<T, Compact> (in.getSize ())
{
   if (size == 0) return;
   alloc (size);
   typename SxList<T>::ConstIterator inIt = in.begin ();
   for (int i = 0; i < size; i++, ++inIt)
      elements[i] = *inIt;
}

/// Explicitly sized vector
template<class T>
SxVector<T>::SxVector (ssize_t n)
   : SxVecRef<T,Compact> (n)
{
   if (n > 0) alloc(n);
}

/// Explicitly sized vector
template<class T>
SxVector<T>::SxVector (ssize_t nRowsIn, ssize_t nColsIn)
   : SxVecRef<T,Compact> (nRowsIn * nColsIn)
{
   SX_CHECK (nRowsIn * nColsIn > 0);
   alloc(nRowsIn * nColsIn);
   this->reshape (nRowsIn, nColsIn);
}

/// Copy data
//inline SxVector<T> getCopy () const;

// --- change size
/// Change size (loses data)
template<class T>
void SxVector<T>::resize (ssize_t n)
{
   if (size == n && isOnlyRef ()) return;
   this->unref_ ();
   SxVecLayout<Compact>::set (n, 1);
   if (n > 0) {
      alloc (n);
   } else {
      allocMem = NULL;
      elements = NULL;
      return;
   }
}

/** \brief Change size (keep and/or initialize data)
    @param keep    whether current values should be copied
    @param fillVal how to set new values
  */
template<class T>
void SxVector<T>::resize (ssize_t n, bool keep, const T& fillVal)
{
   if (n == 0) { this->unref (); return; }
   if (size == n && isOnlyRef ())
   {
      if (!keep) SxVecRef<T,Compact>::set (fillVal);
      return;
   }
   SxAllocMem* newAlloc = SxAllocMem::create (sizeof(T) * n);
   ssize_t i = 0;
   if (keep && elements && size > 0)  {
      // copy current data
      T* newElem = (T*)newAlloc->mem;
      ssize_t copySize = min(n,size);
      for ( ; i < copySize; ++i) newElem[i] = elements[i];
   }
   // deallocate old allocation
   this->unref_ ();
   // set new allocation
   allocMem = newAlloc;
   elements = (T*)newAlloc->mem;
   SxVecLayout<Compact>::set (n, 1);
   // --- fill
   for ( ; i < n; ++i) elements[i] = fillVal;
}

/// Change size (loses data)
template<class T>
void SxVector<T>::reformat (ssize_t nRowsIn, ssize_t nColsIn)
{
   SX_CHECK (nRowsIn > 0, nRowsIn);
   SX_CHECK (nColsIn > 0, nColsIn);
   nRows = nRowsIn;
   nCols = nColsIn;
   ssize_t newSize = nRows * nCols;
   if (size != newSize || !isOnlyRef ())  {
      SX_CHECK (size <= 0 || allocMem, size);
      this->unref_ ();
      alloc (newSize);
   }
   size = newSize;
}

// --- assignment
/// Copy data
template<class T>
inline void SxVector<T>::copy (const SxVecRef<T,Compact> &in)
{
   reformat (in.getNRows (), in.getNCols ());
   memcpy (elements, in.elements, sizeof(T) * getSize ());
   this->auxData = in.auxData;
}

/// Copy data
template<class T>
template <class T2, StrideType Layout>
inline void SxVector<T>::copy (const SxVecRef<T2,Layout> &in)
{
   if (!isPacked(Layout))  {
      reformat (in.getNRows (), in.getNCols ());
      SxVecIt<T2,Layout> it = in.begin ();
      for (T* dest = elements; it != in.end; ++it)
         *dest++ = *it;
   } else {
      // TODO: expand?
      SX_EXIT;
   }
}

/// Move assignment: take over from rvalue reference
template<class T>
inline SxVector<T>& SxVector<T>::operator= (SxVecRef<T, Compact> &&in)
{
   SX_CHECK (this != &in); // X = std::move(X) is insane
   if (!in.allocMem)  {
      this->unref ();
      return *this;
   }
   if (in.allocMem == this->allocMem)  {
      // handle cases like
      //   x = x(SxIdx(...))
      // where the RHS may reference x entirely
      this->unref ();
   }
   this->auxData = in.auxData;
   if (in.isOnlyRef ())  {
      this->unref_ ();
      allocMem = in.allocMem;
      elements = in.elements;
      SxVecLayout<Compact>::set (in.getNRows (), in.getNCols ());
      in.allocMem = NULL;
      in.elements = NULL;
      in.SxVecLayout<Compact>::unset ();
   } else {
      if (getSize () != in.getSize () || !isOnlyRef ())  {
         this->unref_ ();
         alloc (in.getSize ());
      }
      SxVecLayout<Compact>::set (in.getNRows (), in.getNCols ());
      SxVecRef<T,Compact>::operator<<= (in);
   }
   return *this;
}

// set values from stack
template<class T>
void SxVector<T>::set (const SxStack<T> &stack, size_t stackSize,
                       ssize_t offset)
{
   SX_CHECK (offset >= 0, offset);
   SX_CHECK ((size_t)getSize () >= size_t(offset) + stackSize,
             getSize (), offset, stackSize);
   stack.exportStack(elements + offset, stackSize);
}

// --- functions that make sense on 1D vectors
/// Find minimum value
template<class T>
T SxVector<T>::minval (ssize_t *idx) const
{
   T res = elements[0];
   ssize_t i = 0, minIdx = 0;
   for (SxVecConstIt<T,Compact> it = this->begin (); it != this->end (); it++, i++)
   {
      if (*it < res)  {
         minIdx = i;
         res = *it;
      }
   }
   if (idx) *idx = minIdx;
   return res;
}
/// Find maximum value
template<class T>
T SxVector<T>::maxval (ssize_t *idx) const
{
   T res = elements[0];
   ssize_t i = 0, maxIdx = 0;
   for (SxVecConstIt<T,Compact> it = this->begin (); it != this->end (); it++, i++)
   {
      if (*it > res)  {
         maxIdx = i;
         res = *it;
      }
   }
   if (idx) *idx = maxIdx;
   return res;
}

/// Get sorting index
template<class T>
SxArray<ssize_t> SxVector<T>::getSortIdx () const
{
   SxArray<ssize_t> sortIdx(getSize ());
   if (getSize() == 1)  {
      sortIdx(0) = 0;
      return sortIdx;
   }

   ::sort (sortIdx.elements, elements, getSize());
   return sortIdx;
}

/// Reorder data (in-place)
template<class T, StrideType Layout>
void SxVecRef<T, Layout>::sortByIdx (const SxArray<ssize_t> &sortIdx)
{
   SX_CHECK (getSize () == sortIdx.getSize (),
             getSize (), sortIdx.getSize ());
   SX_CHECK (getSize () > 0, getSize ());
   if (getSize () == 1) return;
   if (Layout == Compact && isOnlyRef ()) {
      // shift data to a temporary vector and return sorted in a new allocation
      *this = SxVector<T> (std::move (*this)).getSorted (sortIdx);
   } else {
      // get memory from allocation cache
      T* tmp = (T*)SxAllocation::get (sizeof(T) * getSize ());

      // --- copy elements in order
      SxVecIt<T, Layout> it = begin ();
      for (ssize_t i = 0; i < getSize (); ++i)
         tmp[i] = *it++;
      SX_CHECK (it == end ());

      // --- fill in copied elements in new order
      it = begin ();
      for (ssize_t i = 0; i < getSize (); ++i)
         elements[i] = tmp[sortIdx(i)];
      // return memory
      SxAllocation::retain (tmp, sizeof(T) * getSize ());
   }
}

/// 1D Integrate (Simpson rule)
template<class T>
T SxVector<T>::integrate (double step) const
{
   ssize_t n = getSize ();
   SX_CHECK (n > 2, n);
   SX_CHECK (step > 0., step);
   T res = 0;
   const double third  = 1./3., third2 = 2./3., third4 = 4./3.,
                eight3 = 3./8., eight9 = 9./8.;
   if ((n & 1) == 0)  {
      // Simpson 3/8 rule for last 4 points
      res = eight3 * (elements[n-1] + elements[n-4])
          + eight9 * (elements[n-2] + elements[n-3]);
      n -= 3;
   }
   res += third * elements[0];
   for (ssize_t i = 1; i < n; i += 2)  {
      res += third4 * elements[i]
           + third2 * elements[i+1];
   }
   res -= third * elements[n-1];
   res *= step;
   SX_CHECK_NUM(res);
   return res;
}

// Invert a square matrix
template<class T>
SxVector<T> SxVector<T>::inverse () const&
{
   SX_CHECK (getSize () > 0);

   // only square matrices here; to get pseudo-inverse, use "solve"
   SX_CHECK (nRows == nCols, nRows, nCols);
   SX_CHECK (int(nRows) == nRows, nRows);
   SX_CHECK (int(nCols) == nCols, nCols);
   SxAllocMem *resMem = SxAllocMem::create (sizeof(T) * getSize ());
   memcpy (resMem->mem, elements, sizeof(T) * getSize ());
   matInverse ((T*)resMem->mem, (int)nRows, (int)nCols);
   resMem->refCounter--; // down to zero, but don't deallocate...
   // ... because it will be incremented again here:
   return SxVector<T> (resMem, nRows, nCols, this->auxData);
}

// Invert a about-to-die square matrix
template<class T>
SxVector<T> SxVector<T>::inverse () &&
{
   SX_CHECK (getSize () > 0);
   // only square matrices here; to get pseudo-inverse, use "solve"
   SX_CHECK (nRows == nCols, nRows, nCols);

   if (!isOnlyRef ()) return static_cast<const SxVector<T>*>(this)->inverse ();
   // cout << "MOVE inverse" << endl;
   SX_CHECK (int(nRows) == nRows, nRows);
   SX_CHECK (int(nCols) == nCols, nCols);
   matInverse (elements, (int)nRows, (int)nCols);
   return std::move(*this);
}

template<class T>
SxVector<T>
SxVector<T>::pseudoInverse (typename SxTypeMapper<T>::TReal threshold) const&
{
   SX_CHECK(getSize() > 0);
   SX_CHECK(threshold > 0.);
   SX_CHECK(nRows == int(nRows), nRows);
   SX_CHECK(nCols == int(nCols), nCols);
   typedef typename SxTypeMapper<T>::TReal TRe;

   SxVector<T> res(nCols, nRows);
   res.auxData = this->auxData;
   res.auxData.basisPtr = NULL;

   // get a copy of current matrix
   memcpy (res.elements, elements, getSize () * sizeof(T));

   ssize_t minMN = nRows < nCols ? nRows : nCols;

   SxAllocMem *left   = SxAllocMem::create (nRows * minMN * sizeof(T));
   SxAllocMem *rightH = SxAllocMem::create (nCols * minMN * sizeof(T));
   SxAllocMem *valsM  = SxAllocMem::create (minMN * sizeof(TRe));
   TRe *vals = (TRe*)valsM->mem;

   singularValueDecomp (res.elements, (int)nRows, (int)nCols,
                        vals, (T*)left->mem, (T*)rightH->mem,
                        false);

   if (vals[0] < threshold)  {
      cout << "Warning: pseudoInverse with small singular value " << vals[0]
           << " < " << threshold << endl;
   }
   // find number of singular values larger than threshold * max value
   ssize_t nSing = 0;
   for ( ; nSing < minMN ; nSing++) if (vals[nSing] < vals[0] * threshold) break;

   // scale columns of U -> U Sigma^-1
   T* U = (T*)left->mem;
   for (ssize_t iSing = 0; iSing < nSing; iSing++)  {
      TRe inv = 1. / vals[iSing];
      for (ssize_t ir = 0; ir < nRows; ++ir)
         *U++ *= inv;
   }
   // adjoint -> Sigma^-1 U^H with res.elements as workspace
   // TODO: in-place adjoint
   SxVecCompute<Simd>::adjoint (nRows, nSing, (T*)left->mem, nRows,
                                res.elements, nSing);
   // copy back to left->mem
   memcpy (left->mem, res.elements, nSing * nRows * sizeof(T));
   // multiply Sigma^-1 U^H from left with (V^H)^H=V
   matovlp (res.elements, (T*)rightH->mem, (T*)left->mem,
            (int)minMN, (int)nCols, (int)nSing, (int)nRows, (int)nSing);

   valsM->unref ();
   SxAllocMem::destroy (valsM);
   left->unref ();
   SxAllocMem::destroy (left);
   rightH->unref ();
   SxAllocMem::destroy (rightH);
   return res;
}

// solve linear least square via singular value decomposition
// version: keep original A and b
template<class T>
SxVector<T> SxVector<T>::solve (const SxVector<T> &b) const&
{
   SX_CHECK(getSize() > 0);
   SX_CHECK(b.getSize() > 0);
   SX_CHECK(nRows == b.nRows, nRows, b.nRows);
   SX_CHECK(nRows == int(nRows), nRows);
   SX_CHECK(nCols == int(nCols), nCols);


   // get a copy of b
   SxAllocMem *res = SxAllocMem::create (b.getSize () * sizeof(T));
   memcpy (res->mem, b.elements, b.getSize () * sizeof(T));

   // get a copy of current matrix
   SxAllocMem *work = SxAllocMem::create (getSize () * sizeof(T));
   memcpy (work->mem, elements, getSize () * sizeof(T));

   solveLinEq((T*)work->mem, int(nRows), int(nCols),
              (T*)res->mem, int(b.nCols));
   work->unref ();
   SxAllocMem::destroy (work);

   if (nCols == nRows)  {
      // --- result work space is fully used: return it as a vector
      res->refCounter--;
      return SxVector<T> (res, nCols, b.nCols, SxAuxData ());
   }

   // --- copy result from result work space
   SxVector<T> result(nCols * b.nCols);
   result.reshape (nCols, b.nCols);
   size_t colBytes = nCols * sizeof(T);
   for(ssize_t c = 0; c < b.nCols; c++)   {
      memcpy (result.elements + c * nCols, (T*)res->mem + c * nRows, colBytes);
   }
   res->unref ();
   SxAllocMem::destroy (res);

   SX_VALIDATE_VECTOR(result);
   return result;
}

// solve linear least square via singular value decomposition
// version: keep original b, A is no longer needed
template<class T>
SxVector<T> SxVector<T>::solve (const SxVector<T> &b) &&
{
   SX_CHECK(getSize() > 0);
   SX_CHECK(b.getSize() > 0);
   SX_CHECK(nRows == b.nRows, nRows, b.nRows);
   SX_CHECK(nRows == int(nRows), nRows);
   SX_CHECK(nCols == int(nCols), nCols);

   // cout << "solve: temp A" << endl;

   // get a copy of b
   SxAllocMem *res = SxAllocMem::create (b.getSize () * sizeof(T));
   memcpy (res->mem, b.elements, b.getSize () * sizeof(T));

   solveLinEq(elements, int(nRows), int(nCols), (T*)res->mem, int(b.nCols));
#ifndef NDEBUG
   ssize_t nC = nCols, nR = nRows;
   // invalidate this
   this->unref ();
   nCols = nC; // needed below
   nRows = nR;
#endif

   if (nCols == nRows)  {
      // --- result work space is fully used: return it as a vector
      res->refCounter--;
      return SxVector<T> (res, b.nRows, b.nCols, SxAuxData ());
   }

   // --- copy result from result work space
   SxVector<T> result(nCols * b.nCols);
   result.reshape (nCols, b.nCols);
   size_t colBytes = nCols * sizeof(T);
   for(ssize_t c = 0; c < b.nCols; c++)   {
      memcpy (result.elements + c * nCols, (T*)res->mem + c * nRows, colBytes);
   }
   res->unref ();
   SxAllocMem::destroy (res);
#ifndef NDEBUG
   nRows = nCols = 0;
#endif

   SX_VALIDATE_VECTOR(result);
   return result;
}


// solve linear least square via singular value decomposition
// version: keep original A, b is no longer needed
template<class T>
inline SxVector<T> SxVector<T>::solve (SxVector<T> &&b) const&
{
   SX_CHECK(getSize() > 0);
   SX_CHECK(b.getSize() > 0);
   SX_CHECK(nRows == b.nRows, nRows, b.nRows);
   SX_CHECK(nRows == int(nRows), nRows);
   SX_CHECK(nCols == int(nCols), nCols);

   // cout << "solve: temp b" << endl;

   // get a copy of current matrix
   SxAllocMem *work = SxAllocMem::create (getSize () * sizeof(T));
   memcpy (work->mem, elements, getSize () * sizeof(T));

   solveLinEq((T*)work->mem, int(nRows), int(nCols), b.elements, int(b.nCols));
   work->unref ();
   SxAllocMem::destroy (work);

   if (nCols == nRows)  {
      // --- result work space is fully used: return it as a vector
      b.auxData = SxAuxData ();
      return SxVector<T> (std::move (b));
   }

   // --- copy result from result work space
   SxVector<T> result(nCols * b.nCols);
   result.reshape (nCols, b.nCols);
   size_t colBytes = nCols * sizeof(T);
   for(ssize_t c = 0; c < b.nCols; c++)   {
      memcpy (result.elements + c * nCols, b.elements + c * nRows, colBytes);
   }
#ifndef NDEBUG
   // invalidate b
   b.unref ();
#endif

   SX_VALIDATE_VECTOR(result);
   return result;
}

// solve linear least square via singular value decomposition
// version: A and b are no longer needed
template<class T>
inline SxVector<T> SxVector<T>::solve (SxVector<T> &&b) &&
{
   SX_CHECK(getSize() > 0);
   SX_CHECK(b.getSize() > 0);
   SX_CHECK(nRows == b.nRows, nRows, b.nRows);
   SX_CHECK(nRows == int(nRows), nRows);
   SX_CHECK(nCols == int(nCols), nCols);

   // cout << "solve: temp A & b" << endl;

   solveLinEq(elements, int(nRows), int(nCols), b.elements, int(b.nCols));
#ifndef NDEBUG
   ssize_t nC = nCols, nR = nRows;
   // invalidate this
   this->unref ();
   nCols = nC; // needed below
   nRows = nR;
#endif

   if (nCols == nRows)  {
      // --- result work space is fully used: return it as a vector
      b.auxData = SxAuxData ();
      return SxVector<T> (std::move (b));
   }

   // --- copy result from result work space
   SxVector<T> result(nCols * b.nCols);
   result.reshape (nCols, b.nCols);
   size_t colBytes = nCols * sizeof(T);
   for(ssize_t c = 0; c < b.nCols; c++)   {
      memcpy (result.elements + c * nCols, b.elements + c * nRows, colBytes);
   }
#ifndef NDEBUG
   // invalidate b
   b.unref ();
   nRows = nCols = 0;
#endif

   SX_VALIDATE_VECTOR(result);
   return result;
}

// Cholesky decomposition
template<class T>
SxVector<T> SxVector<T>::choleskyDecomposition (enum UPLO uplo) const
{
   SX_CHECK (getSize () > 0);
   SX_CHECK (nRows == nCols, nRows, nCols);
   SX_CHECK (elements);
   SX_CHECK (int(nRows) == nRows, nRows);
   int n = int(nRows);
   SxVector<T> res(getSize ());
   res.reshape (n,n);
   // no auxData

   if (uplo == LowerLeft)  {
      // fill upper right with 0, copy lower left
      for (int j = 0; j < n; ++j)  {
         T* dest = res.elements + j * n;
         for (int i = 0; i < j; i++)
            *dest++ = 0.;
         T* src = elements + j * (n + 1);
         for (int i = j; i < n; i++)
            *dest++ = *src++;
      }
   } else {
      // copy upper right, fill lower left with 0
      SX_CHECK (uplo == UpperRight);
      for (int j = 0; j < n; ++j)  {
         T* src = elements + j * n;
         T* dest = res.elements + j * n;
         for (int i = 0; i <= j; i++)
            *dest++ = *src++;
         for (int i = j + 1; i < n; i++)
            *dest++ = 0.;
      }
   }
   // call LAPACK wrapper from SxBlasLib
   ::cholesky (res.elements, uplo, NULL, (int)n);
   SX_VALIDATE_VECTOR (res);
   return res;
}


// --- printing
/** Print out a vector either in simple output format or in vector
    format (similar to \URL[Mathematica]{http://www.wolfram.com}) */
//inline void print (bool vectorForm=false) const;
//inline void print (std::ostream &, bool vectorForm=false) const;
/** Recusive print */
//void  println (std::ostream &, bool) const;

// --- scalar - vector operations (scalar first)
template<class T, StrideType Layout>
inline
typename SxVecResult<T, Layout>::VecType
operator+ (T x, const SxVecRef<T, Layout> &y)
{
   return y + x;
}

template<class T, StrideType Layout>
inline
typename SxVecResult<T, Layout>::VecType
operator* (T x, const SxVecRef<T, Layout> &y)
{
   return y * x;
}

template<class T, StrideType Layout>
inline
typename SxVecResult<T, Layout>::VecType
operator- (T x, const SxVecRef<T, Layout> &y)
{
   return SxVecCompute<algoScalar(Layout)>::subtract (x, y);
}

template<class T, StrideType Layout>
inline
typename SxVecResult<T, Layout>::VecType
operator/ (T x, const SxVecRef<T, Layout> &y)
{
   return SxVecCompute<algoScalar(Layout)>::divide (x, y);
}

// --- binary operations: real scalar and complex vector
template<class T, StrideType Layout>
inline
typename SxVecResult<SxComplex<T>, Layout>::VecType
operator+ (T x, const SxVecRef<SxComplex<T>, Layout> &y)
{
   return SxVecCompute<algoScalar(Layout)>::add (y, x); // swap order!
}

template<class T, StrideType Layout>
inline
typename SxVecResult<SxComplex<T>, Layout>::VecType
operator- (T x, const SxVecRef<SxComplex<T>, Layout> &y)
{
   return SxVecCompute<algoScalar(Layout)>::subtract (x, y);
}

template<class T, StrideType Layout>
inline
typename SxVecResult<SxComplex<T>, Layout>::VecType
operator* (T x, const SxVecRef<SxComplex<T>, Layout> &y)
{
   return SxVecCompute<algoScalar(Layout)>::multiply (y, x); // swap order!
}

template<class T, StrideType Layout>
inline
typename SxVecResult<SxComplex<T>, Layout>::VecType
operator/ (T x, const SxVecRef<SxComplex<T>, Layout> &y)
{
   return SxVecCompute<algoScalar(Layout)>::divide (x, y);
}

// --- binary operations: complex scalar and real vector
template<class T, StrideType Layout>
inline
typename SxVecResult<SxComplex<T>, Layout>::VecType
operator+ (const SxComplex<T> &x, const SxVecRef<T, Layout> &y)
{
   return SxVecCompute<algoScalar(Layout)>::add (y, x); // swap order!
}

template<class T, StrideType Layout>
inline
typename SxVecResult<SxComplex<T>, Layout>::VecType
operator- (const SxComplex<T> &x, const SxVecRef<T, Layout> &y)
{
   return SxVecCompute<algoScalar(Layout)>::subtract (x, y);
}

template<class T, StrideType Layout>
inline
typename SxVecResult<SxComplex<T>, Layout>::VecType
operator* (const SxComplex<T> &x, const SxVecRef<T, Layout> &y)
{
   return SxVecCompute<algoScalar(Layout)>::multiply (y, x); // swap order!
}

template<class T, StrideType Layout>
inline
typename SxVecResult<SxComplex<T>, Layout>::VecType
operator/ (const SxComplex<T> &x, const SxVecRef<T, Layout> &y)
{
   return SxVecCompute<algoScalar(Layout)>::divide (x, y);
}


// --- functions
/// Exponential function of vector
template<class T, StrideType Layout>
inline SxVector<T> exp (const SxVecRef<T, Layout> &x)
{
   return SxVecCompute<algoScalar(Layout)>::exp (x);
}

/// Exponential function of vector (reuse memory, if possible)
template<class T, StrideType Layout>
inline typename SxVecResult<T, Layout>::VecType exp (SxVecRef<T, Layout> &&x)
{
   return SxVecCompute<algoScalar(Layout)>::expInPlace (std::move(x));
}

// --- The next 56 lines were generated from snippets/SxVector.h snippet func_of_double
/// Square root function of vector
template<StrideType Layout>
inline SxVector<double> sqrt (const SxVecRef<double, Layout> &x)
{
   return SxVecCompute<algoScalar(Layout)>::sqrt (x);
}

/// Square root function of vector (reuse memory, if possible)
template<StrideType Layout>
inline typename SxVecResult<double, Layout>::VecType sqrt (SxVecRef<double, Layout> &&x)
{
   return SxVecCompute<algoScalar(Layout)>::sqrtInPlace (std::move(x));
}

/// Logarithm function of vector
template<StrideType Layout>
inline SxVector<double> log (const SxVecRef<double, Layout> &x)
{
   return SxVecCompute<algoScalar(Layout)>::log (x);
}

/// Logarithm function of vector (reuse memory, if possible)
template<StrideType Layout>
inline typename SxVecResult<double, Layout>::VecType log (SxVecRef<double, Layout> &&x)
{
   return SxVecCompute<algoScalar(Layout)>::logInPlace (std::move(x));
}

/// Error function of vector
template<StrideType Layout>
inline SxVector<double> erf (const SxVecRef<double, Layout> &x)
{
   return SxVecCompute<algoScalar(Layout)>::erf (x);
}

/// Error function of vector (reuse memory, if possible)
template<StrideType Layout>
inline typename SxVecResult<double, Layout>::VecType erf (SxVecRef<double, Layout> &&x)
{
   return SxVecCompute<algoScalar(Layout)>::erfInPlace (std::move(x));
}

/// Complementary error function of vector
template<StrideType Layout>
inline SxVector<double> erfc (const SxVecRef<double, Layout> &x)
{
   return SxVecCompute<algoScalar(Layout)>::erfc (x);
}

/// Complementary error function of vector (reuse memory, if possible)
template<StrideType Layout>
inline typename SxVecResult<double, Layout>::VecType erfc (SxVecRef<double, Layout> &&x)
{
   return SxVecCompute<algoScalar(Layout)>::erfcInPlace (std::move(x));
}
// --- func_of_double
/// Power function of vector x
template<StrideType Layout>
inline SxVector<double> pow (const SxVecRef<double, Layout> &x, double a)
{
   // --- treat a few special cases
   if (a == double(0))  {
      SxVector<double> res(x.getSize ());
      res.reshape (x.getNRows (), x.getNCols ());
      res.set (1.);
      res.auxData = x.auxData;
      return res;
   } else if (a == double(1)) {
      return x;
   } else if (a == double(2)) {
      return x.sqr ();
   } else if (a == double(3)) {
      return x.cub ();
   }
   return SxVecCompute<UseIterator>::pow (x,a);
}


#include <SxVecCompute.hpp>

template <class T>
SxVector<T> mergeCols (const SxVecRef<T,Compact> &x, const SxVecRef<T,Compact> &y)
{
   SX_CHECK (x.getSize () > 0, x.getSize ());
   SX_CHECK (y.getSize () > 0, y.getSize ());
   SX_CHECK (x.getNRows () == y.getNRows (), x.getNRows (), y.getNRows ());
   int n = (int)x.getNRows (), nx = (int)x.getNCols (), ny = (int)y.getNCols ();
   int nxy = nx + ny;
   SxVector<T> res(n, nxy);
   res(SxIdx(0, n * nx - 1))       <<= x;
   res(SxIdx(n * nx, n * nxy - 1)) <<= y;
   res.auxData = x.auxData;
   return res;
}

#include <SxLoopMPI.h>
#ifdef USE_LOOPMPI
template<>
inline void SxLoopMPI::sum (SxVector<int> &inout)
{
   // actual implementation is in SxLoopMPI
   SxLoopMPI::sum (inout.elements, inout.elements, inout.getSize ());
}

template<>
inline void SxLoopMPI::sum (SxVector<double> &inout)
{
   // actual implementation is in SxLoopMPI
   SxLoopMPI::sum (inout.elements, inout.elements, inout.getSize ());
}

template<>
inline void SxLoopMPI::sum (SxVector<SxComplex16> &inout)
{
   // actual implementation is in SxLoopMPI
   SxLoopMPI::sum ((double*)inout.elements, (double*)inout.elements,
                   inout.getSize () * 2);
}

template<>
inline void SxLoopMPI::sum (SxVecRef<double,PackedSymMatrix> &inout)
{
   // actual implementation is in SxLoopMPI
   SxLoopMPI::sum (inout.elements, inout.elements, inout.getSize ());
}


template <>
inline void SxLoopMPI::bcast (SxVector<double> &inout, int source)
{
   SxLoopMPI::bcast(inout.elements, inout.getSize(), source);
}

template <>
inline void SxLoopMPI::bcast (SxVector<SxComplex16> &inout, int source)
{
   SxLoopMPI::bcast((double*)inout.elements, 2 * inout.getSize(), source);
}

template <>
inline void SxLoopMPI::bcast (SxVecRef<double> &inout, int source)
{
   SxLoopMPI::bcast(inout.elements, inout.getSize(), source);
}

template <>
inline void SxLoopMPI::bcast (SxVecRef<SxComplex16> &inout, int source)
{
   SxLoopMPI::bcast((double*)inout.elements, 2 * inout.getSize(), source);
}
#else
template<> inline void SxLoopMPI::sum (SxVector<int> &) { }
template<> inline void SxLoopMPI::sum (SxVector<double> &) { }
template<> inline void SxLoopMPI::sum (SxVector<SxComplex16> &) { }
template<> inline void SxLoopMPI::sum (SxVecRef<double,PackedSymMatrix> &) { }
template<> inline void SxLoopMPI::bcast (SxVector<double> &, int) { }
template<> inline void SxLoopMPI::bcast (SxVector<SxComplex16> &, int) { }
template<> inline void SxLoopMPI::bcast (SxVecRef<double> &, int) { }
template<> inline void SxLoopMPI::bcast (SxVecRef<SxComplex16> &, int) { }
#endif


#endif /* _SX_VECTOR_H_ */
