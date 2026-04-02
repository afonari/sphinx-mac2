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

#ifndef _SX_VEC_LAYOUT_H_
#define _SX_VEC_LAYOUT_H_

#include <SxMath.h>

#ifdef __GNUC__
#define forced_inline __attribute__((always_inline))
#else
#define forced_inline
#endif

/// Types of memory layout
enum StrideType { Compact, Strided, SubMatrix, GeneralMatrix, PackedSymMatrix, BandedMatrix };

constexpr bool isPacked (StrideType layout)
{
   return (layout == PackedSymMatrix || layout == BandedMatrix);
}

constexpr bool isContiguous (StrideType layout)
{
   return (   layout == Compact
           || layout == PackedSymMatrix
           || layout == BandedMatrix);
}

constexpr StrideType realPartLayout(StrideType layout)
{
   return  (layout == Compact   || layout == Strided)       ? Strided
        : ((layout == SubMatrix || layout == GeneralMatrix) ? GeneralMatrix
             /* layout is packed  */                        : layout);

}

constexpr StrideType colStrideType(StrideType layout)
{
   return (layout == Strided || layout == GeneralMatrix) ? Strided : Compact;
}

constexpr StrideType subMatrixLayout (StrideType layout)
{
   if (layout == Compact || layout == SubMatrix) return SubMatrix;
   if (layout == Strided || layout == GeneralMatrix) return GeneralMatrix;
   return layout;
}

// --- Layout class (independent of scalar type)
template<StrideType Layout> class SxVecLayout;
template <> class SxVecLayout<Compact>
{
   public:
      /// Number of rows
      ssize_t nRows;
      /// Number of columns
      ssize_t nCols;
      /// Number of elements
      ssize_t size;

      /// Get row stride
      constexpr ssize_t getRowStride () const { return 1; }
      /// Get column stride
      inline ssize_t getColStride () const { return nRows; }

      /// Constructor with sizes
      inline SxVecLayout (ssize_t nRowsIn, ssize_t nColsIn, ssize_t = 0, ssize_t = 0)
         : nRows(nRowsIn), nCols(nColsIn),
           size (nRowsIn * nColsIn)
      { /* empty */ }

      /// Copy constructor
      inline SxVecLayout (const SxVecLayout<Compact> &) = default;

      /// Constructor from a different Layout
      template<StrideType Layout2>
      inline SxVecLayout(const SxVecLayout<Layout2> &in)
        : nRows(in.nRows), nCols(in.nCols)
      {
         size = nRows * nCols;
      }

      /// Set sizes
      inline void set (ssize_t nRowsIn, ssize_t nColsIn, ssize_t = 1, ssize_t = 1)
      {
         nRows = nRowsIn;
         nCols = nColsIn;
         size  = nRowsIn * nColsIn;
      }

      /// Invalidate
      inline void unset ()
      {
         nRows = nCols = size = 0;
      }

      /// Empty constructor
      inline SxVecLayout () : nRows(0), nCols (0), size (0) {}

      /// Set strides (not used)
      inline void setStrides (ssize_t, ssize_t) {}

   protected:
      /// Linear access (index -> memory offset)
      ssize_t getIdx (ssize_t i) const forced_inline { return i; }
      /// 2-index access ((row,col) -> memory offset)
      ssize_t getIdx (ssize_t i, ssize_t j) const forced_inline
      {
         return i + nRows * j;
      }
};

template <> class SxVecLayout<Strided>
{
   public:
      /// Number of rows
      ssize_t nRows;
      /// Number of columns
      ssize_t nCols;
      /// Number of elements
      ssize_t size;

      /// Stride
      ssize_t rowStride;

      /// Get row stride
      inline ssize_t getRowStride () const { return rowStride; }
      /// Get column stride
      inline ssize_t getColStride () const { return rowStride * nRows; }

      /// Constructor with sizes
      inline SxVecLayout (ssize_t nRowsIn, ssize_t nColsIn, ssize_t stride, ssize_t)
         : nRows(nRowsIn), nCols(nColsIn),
           size (nRowsIn * nColsIn), rowStride (stride)
      { /* empty */ }

      /// Set sizes
      inline void set (ssize_t nRowsIn, ssize_t nColsIn, ssize_t stride, ssize_t)
      {
         nRows = nRowsIn;
         nCols = nColsIn;
         size  = nRowsIn * nColsIn;
         rowStride = stride;
      }

      /// Invalidate
      inline void unset ()
      {
         nRows = nCols = size = 0;
      }

      /// Empty constructor
      inline SxVecLayout () : nRows(0), nCols (0), size (0) {}

      /// Set row strides (2nd stride is not used)
      inline void setStrides (ssize_t stride, ssize_t) { rowStride = stride; }
   protected:
      /// Linear access (index -> memory offset)
      ssize_t getIdx (ssize_t i) const forced_inline { return i * rowStride; }
      /// 2-index access ((row,col) -> memory offset)
      ssize_t getIdx (ssize_t i, ssize_t j) const forced_inline
      {
         return (i + nRows * j) * rowStride;
      }
};

template <> class SxVecLayout<SubMatrix>
{
   public:
      /// Number of rows
      ssize_t nRows;
      /// Number of columns
      ssize_t nCols;
      /// Number of elements
      ssize_t size;

      /// Allocated length of columns (=BLAS3 "lda")
      ssize_t colStride;

      /// Get row stride
      constexpr ssize_t getRowStride () const { return 1; }
      /// Get column stride
      inline ssize_t getColStride () const { return colStride; }

      /// Constructor with sizes
      inline SxVecLayout (ssize_t nRowsIn, ssize_t nColsIn, ssize_t, ssize_t stride)
         : nRows(nRowsIn), nCols(nColsIn),
           size (nRowsIn * nColsIn), colStride (stride)
      { /* empty */ }

      /// Set sizes
      inline void set (ssize_t nRowsIn, ssize_t nColsIn, ssize_t, ssize_t stride)
      {
         nRows = nRowsIn;
         nCols = nColsIn;
         size  = nRowsIn * nColsIn;
         colStride = stride;
      }

      /// Invalidate
      inline void unset ()
      {
         nRows = nCols = size = 0;
      }

      /// Empty constructor
      inline SxVecLayout () : nRows(0), nCols (0), size (0) {}

      /// Set col stride (1st stride is not used)
      inline void setStrides (ssize_t, ssize_t stride) { colStride = stride; }
   protected:
      /// Linear access (index -> memory offset)
      ssize_t getIdx (ssize_t i) const forced_inline
      {
         return i + (colStride - nRows) * (i / nRows);
      }
      /// 2-index access ((row,col) -> memory offset)
      ssize_t getIdx (ssize_t i, ssize_t j) const forced_inline
      {
         return i + colStride * j;
      }
};

template <> class SxVecLayout<GeneralMatrix>
{
   public:
      /// Number of rows
      ssize_t nRows;
      /// Number of columns
      ssize_t nCols;
      /// Number of elements
      ssize_t size;

      /// Row stride
      ssize_t rowStride;
      /// Column stride
      ssize_t colStride;

      /// Get row stride
      inline ssize_t getRowStride () const { return rowStride; }
      /// Get column stride
      inline ssize_t getColStride () const { return colStride; }

      /// Constructor with sizes
      inline SxVecLayout (ssize_t nRowsIn, ssize_t nColsIn,
                          ssize_t rowStrideIn, ssize_t colStrideIn)
         : nRows(nRowsIn), nCols(nColsIn), size (nRowsIn * nColsIn),
           rowStride (rowStrideIn), colStride (colStrideIn)
      { /* empty */ }

      /// Set sizes
      inline void set (ssize_t nRowsIn, ssize_t nColsIn,
                       ssize_t rowStrideIn, ssize_t colStrideIn)
      {
         nRows = nRowsIn;
         nCols = nColsIn;
         size  = nRowsIn * nColsIn;
         rowStride = rowStrideIn;
         colStride = colStrideIn;
      }

      /// Invalidate
      inline void unset ()
      {
         nRows = nCols = size = 0;
      }

      /// Empty constructor
      inline SxVecLayout () : nRows(0), nCols (0), size (0) {}

      /// Set strides
      inline void setStrides (ssize_t rowStrideIn, ssize_t colStrideIn)
      {
         rowStride = rowStrideIn;
         colStride = colStrideIn;
      }
   protected:
      /// Linear access (index -> memory offset)
      ssize_t getIdx (ssize_t i) const forced_inline
      {
         return i * rowStride + (colStride - nRows * rowStride) * (i / nRows);
      }
      /// 2-index access ((row,col) -> memory offset)
      ssize_t getIdx (ssize_t i, ssize_t j) const forced_inline
      {
         return i * rowStride + colStride * j;
      }
};

template <> class SxVecLayout<PackedSymMatrix>
{
   public:
      /// Number of rows
      ssize_t nRows;
      /// Number of columns
      ssize_t nCols;
      /// Number of elements
      ssize_t size;

      /// Get row stride (should never be called)
      inline ssize_t getRowStride () const { SX_EXIT; }
      /// Get column stride (should never be called)
      inline ssize_t getColStride () const { SX_EXIT; }

      /// Constructor with sizes
      inline SxVecLayout (ssize_t nRowsIn, ssize_t nColsIn, ssize_t, ssize_t)
         : nRows(nRowsIn), nCols(nColsIn)
      {
         SX_CHECK (nRows == nCols, nRows, nCols);
         size = nRows * (nRows + 1)/2;
      }

      /// Set sizes
      inline void set (ssize_t nRowsIn, ssize_t nColsIn, ssize_t, ssize_t)
      {
         SX_CHECK (nRowsIn == nColsIn, nRowsIn, nColsIn);
         nCols = nRows = nRowsIn;
         size = nRows * (nRows + 1)/2;
      }

      /// Invalidate
      inline void unset ()
      {
         nRows = nCols = size = 0;
      }

      /// Empty constructor
      inline SxVecLayout () : nRows(0), nCols (0), size (0) {}

      /// Set strides (not used)
      inline void setStrides (ssize_t, ssize_t) { }
   protected:
      /// Linear access (index -> memory offset)
      ssize_t getIdx (ssize_t i) const forced_inline { return i; }
      /// 2-index access ((row,col) -> memory offset)
      ssize_t getIdx (ssize_t i, ssize_t j) const forced_inline
      {
         SX_CHECK (i <= j, i, j); // only upper
         return i + (j * (j + 1)) / 2;
      }
};

template <> class SxVecLayout<BandedMatrix>
{
   public:
      /// Number of rows
      ssize_t nRows;
      /// Number of columns
      ssize_t nCols;
      /// Number of elements
      ssize_t size;

      /// Number of subdiagonals
      int nSubDiag;
      /// Number of superdiagonals
      int nSupDiag;

      /// Get row stride (should never be called)
      inline ssize_t getRowStride () const { SX_EXIT; }
      /// Get column stride (should never be called)
      inline ssize_t getColStride () const { SX_EXIT; }

      /// Constructor with sizes
      inline SxVecLayout (ssize_t nRowsIn, ssize_t nColsIn, ssize_t nSub, ssize_t nSup)
         : nRows(nRowsIn), nCols(nColsIn),
           size ((nSub+nSup+1) * nColsIn),
           nSubDiag(int(nSub)), nSupDiag(int(nSup))
      { /* empty */ }

      /// Set sizes
      inline void set (ssize_t nRowsIn, ssize_t nColsIn, ssize_t nSub, ssize_t nSup)
      {
         nRows = nRowsIn;
         nCols = nColsIn;
         size = (nSub+nSup+1) * nColsIn;
         nSubDiag = int(nSub);
         nSupDiag = int(nSup);
      }

      /// Invalidate
      inline void unset ()
      {
         nRows = nCols = size = 0;
         nSubDiag = nSupDiag = -1;
      }

      /// Empty constructor
      inline SxVecLayout ()
         : nRows(0), nCols (0), size (0), nSubDiag(-1), nSupDiag(-1)
      {}

      /// Set strides (not used)
      inline void setStrides (ssize_t, ssize_t) { }
   protected:
      /// Linear access (index -> memory offset)
      ssize_t getIdx (ssize_t i) const forced_inline { return i; }
      /// 2-index access ((row,col) -> memory offset)
      ssize_t getIdx (ssize_t i, ssize_t j) const forced_inline
      {
         SX_CHECK (i >= j - nSupDiag, i, j, nSupDiag);
         SX_CHECK (i <= j + nSubDiag, i, j, nSubDiag);
         //return (i - j - nSupDiag + j * (nSubDiag + nSupDiag + 1));
         return (i - nSupDiag + j * (nSubDiag + nSupDiag));
      }
};

template<StrideType Layout>
void copy (SxVecLayout<Layout> *dest, const SxVecLayout<Layout> &src)
{
   *dest = src;
}

template<StrideType Layout, StrideType Layout2>
void copy (SxVecLayout<Layout> *dest, const SxVecLayout<Layout2> &src)
{
   dest->set (src.nRows, src.nCols, src.getRowStride (), src.getColStride ());
}


#endif /* _SX_VEC_LAYOUT_H_ */
