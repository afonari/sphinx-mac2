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

#ifndef _SX_VECTOR3_H_
#define _SX_VECTOR3_H_

#include <stdio.h>
#include <math.h>
#include <SxComplex.h>
#include <SxList.h>
#include <SxMathLib.h>
#include <SxArray.h>
#include <iostream>
#include <SxTypeMapper.h>
#include <SxPrecision.h>
#include <SxConstants.h>
#include <SxString.h>


// resolve ambiguity (should never be used)
template<class T>
inline SxComplex<T> operator+(const SxComplex<T> &c, int i)
{
   return c + T(i);
}

/**
  @ingroup Numerics
  \note Different to the vector classes, SxVector3 does not support reference
  counting. Instead, asVec3Ref may be used to make any 3 numbers behave as
  a SxVector3. For this, it is absolutely crucial that this SxVector3 remains
  memory-compatible to (T::Type)[3].
  Under no circumstances, virtual functions or additional data may be added.
  */
template<class T>
class SxVector3
{
   public:
      /** \brief The vector elements */
      T v[3];

      inline SxVector3 ();
      inline explicit SxVector3 (const T &xyz)
      {
         v[0] = xyz, v[1] = xyz, v[2] = xyz;
      }
      inline SxVector3 (const T &xVal, const T &yVal, const T &zVal)
      {
         v[0] = xVal, v[1] = yVal, v[2] = zVal;
      }
      inline SxVector3 (const SxVector3<T> &in) = default;
      template<class B>
         inline SxVector3 (const SxVector3<B> &in);
      inline explicit SxVector3 (const SxList<T> &);
      inline explicit SxVector3 (const SxArray<T> &);
      
      /** \brief Reinterpret 3 numbers in memory as a SxVector3 
        \example
        \code
SxVector<double> m(3,nAtoms);
SxVector3<double> t;
for (int ia = 0; ia < nAtoms; ia++)  {
   // m.colRef(ia) += t; // would not work
   SxVector3<double>::toVec3Ref(m.colRef(ia).elements) += t; // works
}
       \endcode
       @return a SxVector3 reference to the memory pointed to by in.
       */
      static SxVector3<T> & toVec3Ref(T in[3])
      {
         return reinterpret_cast<SxVector3<T> &> (*in); 
      }
      static const SxVector3<T> & toVec3Ref(const T in[3])
      { 
         return reinterpret_cast<const SxVector3<T> &> (*in); 
      }

      inline       T &operator() (ssize_t i);
      inline const T &operator() (ssize_t i) const;
      inline       T &operator() (SxAutoLoop &i);
      inline const T &operator() (SxAutoLoop &i) const;
      /** Assign scalar value to vector */
      inline SxVector3<T>     operator= (const T &s);
      /** Assign vector to vector */
      inline SxVector3<T> &   operator= (const SxVector3<T> &in) = default;
      template<class B>
      inline SxVector3<T> &   operator= (const SxVector3<B> &in);
      inline void             operator+= (const SxVector3<T> &r);
      inline void             operator+= (const T &s);
      inline void             operator-= (const T &s);
      inline void             operator/= (const SxVector3<T> &r);
      inline void             operator*= (const T &s);
      inline void             operator/= (const T &s);
      /** the unary (-) operator */
      inline SxVector3<T>     operator-  () const;
      inline void             operator-= (const SxVector3<T> &r);
      inline bool             operator== (const SxVector3<T> &in) const;
      inline bool operator!= (const SxVector3<T> &in) const {
         return ! operator== (in);
      }

      inline SxVector3<T>  x (const SxVector3<T> &in) const;
      inline SxVector3<typename SxTypeMapper<T>::TReal> absSqr () const;
      inline void set (const T &);
      inline T sum () const;
      inline T product () const;
      /** \brief Normalize the vector
        @return Reference to the normalized vector
        \note The return value allows for on-the-fly normalization
              of expressions:
        \code
n12 = (a1.x(a2)).normalize ();
        \endcode
        \note For getting a normalized vector (without changing the original)
              use
        \code
normalizedCopy = SxVector3<?> (original).normalize ();
        \endcode
        */
      inline SxVector3<T> &normalize ();
      /// Square of length of vector
      inline typename SxTypeMapper<T>::TReal normSqr () const;
      /// Length of vector
      inline typename SxTypeMapper<T>::TReal norm () const;
      /** minimum element */
      inline T minval () const;
      /** maximum element */
      inline T maxval () const;

      /// SxVector3<T> + int
      SxVector3<decltype(T(0) + int(0))>
      operator+ (int);

      void print () const;

};


template<class T>
SxVector3<T>::SxVector3 ()
{
   v[0] = (T)0;
   v[1] = (T)0;
   v[2] = (T)0;
}



template<class T>
   template<class B>
      SxVector3<T>::SxVector3 (const SxVector3<B> &in)
{
   v[0] = (T)in.v[0];
   v[1] = (T)in.v[1];
   v[2] = (T)in.v[2];
}


template<class T>
SxVector3<T>::SxVector3 (const SxList<T> &list)
{
   SX_CHECK (list.getSize() == 3, list.getSize());
   v[0] = list(0); v[1] = list(1); v[2] = list(2);
}

template<class T>
SxVector3<T>::SxVector3 (const SxArray<T> &array)
{
   SX_CHECK (array.getSize() == 3, array.getSize());
   v[0] = array(0); v[1] = array(1); v[2] = array(2);
}

template<class T>
T &SxVector3<T>::operator() (ssize_t i)
{
   SX_CHECK (i >= 0 && i < 3, i); 
   return v[i];
}


template<class T>
const T &SxVector3<T>::operator() (ssize_t i) const
{
   SX_CHECK (i >= 0 && i < 3, i); 
   return v[i];
}

template<class T>
T &SxVector3<T>::operator() (SxAutoLoop &i)
{
   SX_CHECK (i >= 0 && i < 3, i); 
   i.setLimit (3);
   return v[i];
}


template<class T>
const T &SxVector3<T>::operator() (SxAutoLoop &i) const
{
   SX_CHECK (i >= 0 && i < 3, i); 
   i.setLimit (3);
   return v[i];
}


template<class T>
SxVector3<T> SxVector3<T>::operator= (const T &s)
{
   v[0] = v[1] = v[2] = s;
   return *this;
}


template<class T>
template<class B>
SxVector3<T> &SxVector3<T>::operator= (const SxVector3<B> &in)
{
   v[0] = static_cast<T>(in.v[0]);
   v[1] = static_cast<T>(in.v[1]);
   v[2] = static_cast<T>(in.v[2]);
   return *this;
}



template<class T>
void SxVector3<T>::operator+= (const SxVector3<T> &r)
{
   v[0] += r.v[0]; v[1] += r.v[1]; v[2] += r.v[2];
}

template<class T>
void SxVector3<T>::operator+= (const T &s)
{
   v[0] += s; v[1] += s; v[2] += s;
}

template<class T>
void SxVector3<T>::operator/= (const SxVector3<T> &r)
{
   SX_CHECK_DIV (r.v[0]);
   SX_CHECK_DIV (r.v[1]);
   SX_CHECK_DIV (r.v[2]);
   v[0] /= r.v[0]; v[1] /= r.v[1]; v[2] /= r.v[2];
}

template<class T>
void SxVector3<T>::operator*= (const T &s)
{
   v[0] *= s; v[1] *= s; v[2] *= s;
}


template<class T>
void SxVector3<T>::operator/= (const T &s)
{
   SX_CHECK_DIV (s);
   v[0] /= s; v[1] /= s; v[2] /= s;
}


template<class T>
SxVector3<T> SxVector3<T>::operator- () const
{
   return SxVector3<T> (-v[0], -v[1], -v[2]);
}


template<class T>
void SxVector3<T>::operator-= (const SxVector3<T> &r)
{
   v[0] -= r.v[0]; v[1] -= r.v[1]; v[2] -= r.v[2];
}

template<class T>
void SxVector3<T>::operator-= (const T &s)
{
   v[0] -= s; v[1] -= s; v[2] -= s;
}

template<class T>
bool SxVector3<T>::operator== (const SxVector3<T> &in) const
{
   for (int i=0; i<3; i++)
      if (v[i] != in.v[i] )   return false;
   return true;
}

template<>
inline bool SxVector3<float>::operator== (const SxVector3<float> &in) const
{
   for (int i=0; i<3; i++)
      if ( fabs(v[i] - in.v[i]) > 1e-10 )   return false;
   return true;
}

template<>
inline bool SxVector3<double>::operator== (const SxVector3<double> &in) const
{
   for (int i=0; i<3; i++)
      if ( fabs(v[i] - in.v[i]) > 1e-10 )   return false;
   return true;
}


template<class T>
SxVector3<T> SxVector3<T>::x (const SxVector3<T> &in) const
{
   return SxVector3<T> ( v[1]*in.v[2] - v[2]*in.v[1],
		                   v[2]*in.v[0] - v[0]*in.v[2],
		                   v[0]*in.v[1] - v[1]*in.v[0] );
}

// --- absSqr for different types
inline float absSqr (float x)       { return x*x; }
inline double absSqr (double x)     { return x*x; }
inline int absSqr (int i)           { return i*i; }
inline long int absSqr (long int i) { return i*i; }

template<class T>
inline T absSqr (const SxComplex<T> &z)
{ 
   return z.absSqr ();
}
// ------------------------------

template<class T>
SxVector3<typename SxTypeMapper<T>::TReal> SxVector3<T>::absSqr () const
{
   return SxVector3<typename SxTypeMapper<T>::TReal>
      (::absSqr (v[0]), ::absSqr (v[1]), ::absSqr (v[2]));
}

template<class T>
void SxVector3<T>::set (const T &in)
{
   v[0] = in; v[1] = in; v[2] = in;
}



template<class T>
T SxVector3<T>::sum () const
{
   return v[0] + v[1] + v[2];
}


template<class T>
T SxVector3<T>::product () const
{
   return v[0] * v[1] * v[2];
}

template<class T>
SxVector3<T> &SxVector3<T>::normalize ()
{
   typename SxTypeMapper<T>::TReal nrm = norm (), c;
   SX_CHECK_DIV (nrm);
   c = typename SxTypeMapper<T>::TReal(1) / nrm;
   v[0] *= c; v[1] *= c; v[2] *= c;

   SX_CHECK_NUM (v[0]);
   SX_CHECK_NUM (v[1]);
   SX_CHECK_NUM (v[2]);
   
   return *this;
}


template<class T>
void SxVector3<T>::print () const
{
   sxprintf ("[%f, %f, %f]\n", (float)v[0], (float)v[1], (float)v[2]);
}

template<class T>
typename SxTypeMapper<T>::TReal SxVector3<T>::normSqr () const
{
   return ::absSqr(v[0]) + ::absSqr(v[1]) + ::absSqr(v[2]);
}

template<class T>
typename SxTypeMapper<T>::TReal SxVector3<T>::norm () const
{
   return typename SxTypeMapper<T>::TReal(sqrt (normSqr ()));
}

template <>
inline int SxVector3<int>::norm () const
{
   // There is no sqrt for integers
   SX_EXIT;
   return -1;
}

template <class T>
T SxVector3<T>::minval () const
{
   if (v[0] < v[1])  {  // v[0] < v[1]
      if (v[0] < v[2])  return v[0];
      else              return v[2];
   }  else  {           // v[1] <= v[0]
      if (v[1] < v[2])  return v[1];
      else              return v[2];
   }
}

// --- template specializations for complex vectors
template<>
inline SxComplex8 SxVector3<SxComplex8>::minval () const
{
   SX_EXIT; // complex numbers do not have minval
   return SxComplex8 ();
}

template<>
inline SxComplex16 SxVector3<SxComplex16>::minval () const
{
   SX_EXIT; // complex numbers do not have minval
   return SxComplex16 ();
}

template <class T>
inline T SxVector3<T>::maxval () const
{
   if (v[0] > v[1])  {  // v[0] > v[1]
      if (v[0] > v[2])  return v[0];
      else              return v[2];
   }  else  {           // v[1] >= v[0]
      if (v[1] > v[2])  return v[1];
      else              return v[2];
   }
}

// --- template specializations for complex vectors
template<>
inline SxComplex8 SxVector3<SxComplex8>::maxval () const
{
   SX_EXIT; // complex numbers do not have maxval
   return SxComplex8 ();
}

template<>
inline SxComplex16 SxVector3<SxComplex16>::maxval () const
{
   SX_EXIT; // complex numbers do not have maxval
   return SxComplex16 ();
}

template <class T>
SxVector3<decltype(T(0) + int(0))>
SxVector3<T>::operator+ (int i)
{
   return SxVector3<decltype(T(0) + int(0))>
          (v[0] + i, v[1] + i, v[2] + i);
}

// --- The next 36 lines were generated from ./snippets/SxVector3.h snippet VecOpVec
// --- SxVector3<?> = SxVector3<A> + SxVector3<B>
template<class A,class B>
inline SxVector3<decltype(A(0)+B(0))>
operator+ (const SxVector3<A> &a, const SxVector3<B> &b)
{
   return SxVector3<decltype(A(0)+B(0))>
      (a.v[0] + b.v[0], a.v[1] + b.v[1], a.v[2] + b.v[2]);
}

// --- SxVector3<?> = SxVector3<A> - SxVector3<B>
template<class A,class B>
inline SxVector3<decltype(A(0)-B(0))>
operator- (const SxVector3<A> &a, const SxVector3<B> &b)
{
   return SxVector3<decltype(A(0)-B(0))>
      (a.v[0] - b.v[0], a.v[1] - b.v[1], a.v[2] - b.v[2]);
}

// --- SxVector3<?> = SxVector3<A> * SxVector3<B>
template<class A,class B>
inline SxVector3<decltype(A(0)*B(0))>
operator* (const SxVector3<A> &a, const SxVector3<B> &b)
{
   return SxVector3<decltype(A(0)*B(0))>
      (a.v[0] * b.v[0], a.v[1] * b.v[1], a.v[2] * b.v[2]);
}

// --- SxVector3<?> = SxVector3<A> / SxVector3<B>
template<class A,class B>
inline SxVector3<decltype(A(0)/B(0))>
operator/ (const SxVector3<A> &a, const SxVector3<B> &b)
{
   return SxVector3<decltype(A(0)/B(0))>
      (a.v[0] / b.v[0], a.v[1] / b.v[1], a.v[2] / b.v[2]);
}
// --- VecOpVec


// --- The next 180 lines were generated from snippets/SxVector3.h snippet ScalarOpVec
// --- SxVector3<?> = int * SxVector3<B>
template<class B>
inline SxVector3<decltype(int(0)*B(0))>
operator* (int a, const SxVector3<B> &b)
{
   return SxVector3<decltype(int(0)*B(0))>
          (a * b.v[0], a * b.v[1], a * b.v[2]);
}

// --- SxVector3<?> = SxVector3<A> * int
template<class A>
inline SxVector3<decltype(A(0)*int(0))>
operator* (const SxVector3<A> &a, int b)
{
   return SxVector3<decltype(A(0)*int(0))>
          (a.v[0] * b, a.v[1] * b, a.v[2] * b);
}

// --- SxVector3<?> = float * SxVector3<B>
template<class B>
inline SxVector3<decltype(float(0)*B(0))>
operator* (float a, const SxVector3<B> &b)
{
   return SxVector3<decltype(float(0)*B(0))>
          (a * b.v[0], a * b.v[1], a * b.v[2]);
}

// --- SxVector3<?> = SxVector3<A> * float
template<class A>
inline SxVector3<decltype(A(0)*float(0))>
operator* (const SxVector3<A> &a, float b)
{
   return SxVector3<decltype(A(0)*float(0))>
          (a.v[0] * b, a.v[1] * b, a.v[2] * b);
}

// --- SxVector3<?> = double * SxVector3<B>
template<class B>
inline SxVector3<decltype(double(0)*B(0))>
operator* (double a, const SxVector3<B> &b)
{
   return SxVector3<decltype(double(0)*B(0))>
          (a * b.v[0], a * b.v[1], a * b.v[2]);
}

// --- SxVector3<?> = SxVector3<A> * double
template<class A>
inline SxVector3<decltype(A(0)*double(0))>
operator* (const SxVector3<A> &a, double b)
{
   return SxVector3<decltype(A(0)*double(0))>
          (a.v[0] * b, a.v[1] * b, a.v[2] * b);
}

// --- SxVector3<?> = SxComplex8 * SxVector3<B>
template<class B>
inline SxVector3<decltype(SxComplex8(0)*B(0))>
operator* (const SxComplex8 &a, const SxVector3<B> &b)
{
   return SxVector3<decltype(SxComplex8(0)*B(0))>
          (a * b.v[0], a * b.v[1], a * b.v[2]);
}

// --- SxVector3<?> = SxVector3<A> * SxComplex8
template<class A>
inline SxVector3<decltype(A(0)*SxComplex8(0))>
operator* (const SxVector3<A> &a, const SxComplex8 &b)
{
   return SxVector3<decltype(A(0)*SxComplex8(0))>
          (a.v[0] * b, a.v[1] * b, a.v[2] * b);
}

// --- SxVector3<?> = SxComplex16 * SxVector3<B>
template<class B>
inline SxVector3<decltype(SxComplex16(0)*B(0))>
operator* (const SxComplex16 &a, const SxVector3<B> &b)
{
   return SxVector3<decltype(SxComplex16(0)*B(0))>
          (a * b.v[0], a * b.v[1], a * b.v[2]);
}

// --- SxVector3<?> = SxVector3<A> * SxComplex16
template<class A>
inline SxVector3<decltype(A(0)*SxComplex16(0))>
operator* (const SxVector3<A> &a, const SxComplex16 &b)
{
   return SxVector3<decltype(A(0)*SxComplex16(0))>
          (a.v[0] * b, a.v[1] * b, a.v[2] * b);
}

// --- SxVector3<?> = int / SxVector3<B>
template<class B>
inline SxVector3<decltype(int(0)/B(0))>
operator/ (int a, const SxVector3<B> &b)
{
   return SxVector3<decltype(int(0)/B(0))>
          (a / b.v[0], a / b.v[1], a / b.v[2]);
}

// --- SxVector3<?> = SxVector3<A> / int
template<class A>
inline SxVector3<decltype(A(0)/int(0))>
operator/ (const SxVector3<A> &a, int b)
{
   return SxVector3<decltype(A(0)/int(0))>
          (a.v[0] / b, a.v[1] / b, a.v[2] / b);
}

// --- SxVector3<?> = float / SxVector3<B>
template<class B>
inline SxVector3<decltype(float(0)/B(0))>
operator/ (float a, const SxVector3<B> &b)
{
   return SxVector3<decltype(float(0)/B(0))>
          (a / b.v[0], a / b.v[1], a / b.v[2]);
}

// --- SxVector3<?> = SxVector3<A> / float
template<class A>
inline SxVector3<decltype(A(0)/float(0))>
operator/ (const SxVector3<A> &a, float b)
{
   return SxVector3<decltype(A(0)/float(0))>
          (a.v[0] / b, a.v[1] / b, a.v[2] / b);
}

// --- SxVector3<?> = double / SxVector3<B>
template<class B>
inline SxVector3<decltype(double(0)/B(0))>
operator/ (double a, const SxVector3<B> &b)
{
   return SxVector3<decltype(double(0)/B(0))>
          (a / b.v[0], a / b.v[1], a / b.v[2]);
}

// --- SxVector3<?> = SxVector3<A> / double
template<class A>
inline SxVector3<decltype(A(0)/double(0))>
operator/ (const SxVector3<A> &a, double b)
{
   return SxVector3<decltype(A(0)/double(0))>
          (a.v[0] / b, a.v[1] / b, a.v[2] / b);
}

// --- SxVector3<?> = SxComplex8 / SxVector3<B>
template<class B>
inline SxVector3<decltype(SxComplex8(0)/B(0))>
operator/ (const SxComplex8 &a, const SxVector3<B> &b)
{
   return SxVector3<decltype(SxComplex8(0)/B(0))>
          (a / b.v[0], a / b.v[1], a / b.v[2]);
}

// --- SxVector3<?> = SxVector3<A> / SxComplex8
template<class A>
inline SxVector3<decltype(A(0)/SxComplex8(0))>
operator/ (const SxVector3<A> &a, const SxComplex8 &b)
{
   return SxVector3<decltype(A(0)/SxComplex8(0))>
          (a.v[0] / b, a.v[1] / b, a.v[2] / b);
}

// --- SxVector3<?> = SxComplex16 / SxVector3<B>
template<class B>
inline SxVector3<decltype(SxComplex16(0)/B(0))>
operator/ (const SxComplex16 &a, const SxVector3<B> &b)
{
   return SxVector3<decltype(SxComplex16(0)/B(0))>
          (a / b.v[0], a / b.v[1], a / b.v[2]);
}

// --- SxVector3<?> = SxVector3<A> / SxComplex16
template<class A>
inline SxVector3<decltype(A(0)/SxComplex16(0))>
operator/ (const SxVector3<A> &a, const SxComplex16 &b)
{
   return SxVector3<decltype(A(0)/SxComplex16(0))>
          (a.v[0] / b, a.v[1] / b, a.v[2] / b);
}
// --- ScalarOpVec
// --- <?> = SxVector3<A> ^ SxVector3<B>
template<class A,class B>
decltype(A(0)*B(0))
operator^ (const SxVector3<A> &a, const SxVector3<B> &b)
{
   return a.v[0] * b.v[0] + a.v[1] * b.v[1]  + a.v[2] * b.v[2];
}


template<class A,class B>
inline decltype(A(0)*B(0))
dot (const SxVector3<A> &a, const SxVector3<B> &b)
{
   return a ^ b;
}


template<class A,class B>
inline SxVector3<decltype(A(0)*B(0))>
cross (const SxVector3<A> &a, const SxVector3<B> &b)
{
   return SxVector3<decltype(A(0)*B(0))> (
           a.v[1] * b.v[2] - a.v[2] * b.v[1],
           a.v[2] * b.v[0] - a.v[0] * b.v[2],
           a.v[0] * b.v[1] - a.v[1] * b.v[0]);
}

/// int + SxVector3<T>
template <class T>
SxVector3<decltype(int(0)+T(0))>
operator+ (int i, const SxVector3<T> &vec)
{
   return SxVector3<decltype(int(0)+T(0))>
          (vec.v[0] + i, vec.v[1] + i, vec.v[2] + i);
}

inline SxVector3<int> operator%(const SxVector3<int> &a,
                                const SxVector3<int> &b)
{
   SX_CHECK (b(0) != 0 && b(1) != 0 && b(2) != 0, b(0), b(1), b(2));
   return SxVector3<int> (a(0) % b(0), a(1) % b(1), a(2) % b(2));
}

inline SxVector3<int> operator%(const SxVector3<int> &a, int b)
{
   SX_CHECK (b != 0);
   return SxVector3<int> (a(0) % b, a(1) % b, a(2) % b);
}

inline SxVector3<int> operator%(int a, const SxVector3<int> &b)
{
   return SxVector3<int> (a % b(0), a % b(1), a % b(2));
}

inline SxVector3<int> & operator%=(SxVector3<int> &a, const SxVector3<int> &b)
{
   a(0) %= b(0); a(1) %= b(1); a(2) %= b(2);
   return a;
}

inline SxVector3<int> & operator%=(SxVector3<int> &a, int b)
{
   a(0) %= b; a(1) %= b; a(2) %= b;
   return a;
}

//------------------------------------------------------------------------------
inline float getAngle (const SxVector3<float> &a, const SxVector3<float> &b)
{
   float angle = acosf((a ^ b) / sqrtf((a^a) * (b^b)));
   if (angle < 0.f) angle += static_cast<float>(PI);
   return angle;
}

inline double getAngle (const SxVector3<double> &a, const SxVector3<double> &b)
{
   double angle = acos((a ^ b) / sqrt((a^a) * (b^b)));
   if (angle < 0.) angle += PI;
   return angle;
}

//------------------------------------------------------------------------------

inline SxVector3<double> round(const SxVector3<double> &x)
{
   return SxVector3<double> (round(x(0)), round(x(1)), round(x(2)));
}

inline SxVector3<double> floor(const SxVector3<double> &x)
{
   return SxVector3<double> (floor(x(0)), floor(x(1)), floor(x(2)));
}

inline SxVector3<double> ceil(const SxVector3<double> &x)
{
   return SxVector3<double> (ceil(x(0)), ceil(x(1)), ceil(x(2)));
}


//------------------------------------------------------------------------------


template<class T>
inline std::ostream& operator<< (std::ostream &s, const SxVector3<T> &in)
{
   return s << "{" << in.v[0] << "," << in.v[1] << "," << in.v[2] << "}";
}


#endif // _SX_VECTOR3_H_

