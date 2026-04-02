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

#ifndef _SX_SYM_MATRIX_H_
#define _SX_SYM_MATRIX_H_

#include <SxMath.h>
#include <SxVector.h>

namespace {
   /// conjugate of real number = real number
   template <class T>
   T inline _sxconj (T in) { return in; }

   /// conjugate of complex number
   template <class T>
   SxComplex<T> inline _sxconj (const SxComplex<T> &in) { return in.conj (); }
}

/** \brief Packed symmetric/Hermitean matrix

    \author C. Freysoldt, freysoldt@mpie.de */
template <class T>
class SX_EXPORT_MATH SxSymMatrix : public SxVecRef<T, PackedSymMatrix>
{
   public:
      /// Empty constructor
      SxSymMatrix () { }

      /// Sized constructor
      explicit SxSymMatrix (ssize_t squareSize)
         : SxVecRef<T, PackedSymMatrix> (squareSize, 1, squareSize)
      {
         this->alloc (this->size);
      }

      /// Copy from reference
      SxSymMatrix (SxVecRef<T, PackedSymMatrix> &&in)
         : SxVecRef<T, PackedSymMatrix> (in)
      { /* empty */ }

      /// Copy from reference
      inline SxSymMatrix (const SxVecRef<T, PackedSymMatrix> &in);

      /// Copy from reference
      SxSymMatrix (const SxSymMatrix<T> &in)
         : SxSymMatrix (static_cast<const SxVecRef<T,PackedSymMatrix>&>(in))
      { }

      /// Copy from reference
      inline void operator= (const SxVecRef<T,PackedSymMatrix> &in);

      /// Steal data from dying object
      inline void operator= (SxVecRef<T,PackedSymMatrix> &&in)
      {
         if (this == &in) return;
         this->unref_ ();
         SxVecRef<T,PackedSymMatrix>::operator= (std::move (in));
      }

      /// Steal data from dying object
      inline void operator= (SxSymMatrix<T> &&in)
      {
         operator= (static_cast<SxVecRef<T,PackedSymMatrix>&&>(in));
      }

      /// Constructor from list (input: row-wise ordering!)
      explicit SxSymMatrix (const SxList<T> &in)
      {
         // in.size = N*(N+1)/2 ; N^2 < 2 * in.size < (N+1)^2
         ssize_t squareSize = (ssize_t)sqrt(2. * (double)in.getSize ());
         this->resize (squareSize);
         SX_CHECK (this->size == in.getSize (), this->size, in.getSize ());
         typename SxList<T>::ConstIterator it = in.begin ();
         for (ssize_t i = 0; i < squareSize; ++i)
            for (ssize_t j = i; j < squareSize; ++j)
               (*this)(i,j) = *it++;
      }

      /// Resize
      void resize (ssize_t N)
      {
         if (N == this->getSize ()) return;
         SX_CHECK (N >=0, N);
         this->unref_ ();
         this->SxVecLayout<PackedSymMatrix>::set (N, N, -1, -1);
         this->alloc (this->getSize ());
      }

      /// Inverse of matrix
      SxSymMatrix<T> inverse () &&;

      /// Inverse of matrix
      SxSymMatrix<T> inverse () const&
      {
         return SxSymMatrix<T>(*this).inverse ();
      }

      /// Expand to full matrix
      SxVector<T> expand () const;

      void print () const;

   protected:
      /// copy data
      inline void copyData (const SxVecRef<T, PackedSymMatrix> &in)
      {
         ssize_t n = this->getSize ();
         SX_CHECK (in.getSize () == n, in.getSize (), n);
         this->alloc (n);
         memcpy (this->elements, in.elements, n * sizeof(T));
         this->auxData = in.auxData;
      }

};

template<class T>
SxSymMatrix<T>::SxSymMatrix (const SxVecRef<T, PackedSymMatrix> &in)
   : SxVecRef<T, PackedSymMatrix> (in.getNRows (), 1, in.getNCols ())
{
   copyData (in);
}

template<class T>
void SxSymMatrix<T>::operator= (const SxVecRef<T,PackedSymMatrix> &in)
{
   if (this == &in) return;
   this->unref_ ();
   // copy layout
   this->SxVecLayout<PackedSymMatrix>::operator= (in);

   copyData (in);
}

template<class T>
SxSymMatrix<T> SxSymMatrix<T>::inverse () &&
{
   SX_CHECK ((int)this->getNRows () == this->getNRows (), this->getNRows ());
   matInverseTri (this->elements, (int)this->getNRows (), UpperRight);
   return SxSymMatrix<T> (std::move (*this));
}

template<class T> inline T sxconjugate(const T x) { return x; }
template<class T> inline SxComplex<T> sxconjugate(const SxComplex<T>& x) { return x.conj (); }

template<class T>
SxVector<T> SxSymMatrix<T>::expand () const
{
   ssize_t N = this->getNRows ();
   SxVector<T> res(N, N);
   // split into blocks of max 2048 x 2048 (= 4 million elements)
   for (ssize_t jj = 0; jj < N; jj+= 2048)  {
      ssize_t mm = (N-jj) > 2048 ? 2048 : (N-jj);
      for (ssize_t ii = 0; ii <= jj; ii+=2048)  {
         ssize_t nn = (N-ii) > 2048 ? 2048 : (N-ii);
         // --- use cache-oblivious transpose on these blocks
         SxCacheOblivious task(uint16_t(nn), uint16_t(mm), 128);
         while (task.splitProblem ())  {
            ssize_t offX = ii + task.offX;
            ssize_t offY = jj + task.offY;
            if (offX + task.nX <= offY)  {
               for (uint16_t j = 0; j < task.nY; j++)  {
                  for (uint16_t i = 0; i < task.nX; i++)  {
                     const T& x = this->operator()(offX + i,offY + j);
                     res(offX + i,offY + j) = x;
                     res(offY + j,offX + i) = sxconjugate(x);
                  }
               }
            } else if (offX <= offY + task.nY) {
               for (ssize_t j = offY; j < offY + task.nY; j++)  {
                  for (ssize_t i = offX; i < min(offX + task.nX, j+1); i++)  {
                     const T& x = this->operator()(i,j);
                     res(i,j) = x;
                     res(j,i) = sxconjugate(x);
                  }
               }
            }
         }
      }
   }
   res.auxData = this->auxData;
   return res;
}

template<class T>
void SxSymMatrix<T>::print () const
{
   ssize_t n = this->getNCols ();
   for (ssize_t i = 0; i < n; ++i)  {
      // lower part (transpose and conj of upper part)
      cout << (*this)(0,i);
      for (ssize_t j = 1; j <= i; ++j)
         cout << ' ' << _sxconj((*this)(j,i));
      // upper part
      for (ssize_t j = i + 1; j < n; ++j)
         cout << ' ' << (*this)(i,j);
      cout << endl;
   }
}
#endif /* _SX_SYM_MATRIX_H_ */
