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

#ifndef _SX_MATRIX3_H_
#define _SX_MATRIX3_H_

#include <SxVector3.h>
#include <SxComplex.h>
#include <SxError.h>
#include <SxList.h>
#include <SxMathLib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <SxString.h>
#include <SxTypeMapper.h>

/** This class defines operators for simple 3x3 matrices. 
    Instead of the common C order of arrays here the Fortran style, 
    i.e. first dimension varies fastest, is used in order to stay consistent
    with the vector class.
    @ingroup Numerics
    @author Sixten Boeck
 */
template<class T>
class SxMatrix3
{
   public:
      T m[3][3];

      inline SxMatrix3 ();
      /** \brief Constructor as multiple of unit matrix

           val  0   0,
            0  val  0,
            0   0  val
      */
      inline explicit SxMatrix3 (const T &val);
      inline SxMatrix3 (const SxMatrix3<int> &mat);
      inline SxMatrix3 (const SxMatrix3<float> &mat);
      inline SxMatrix3 (const SxMatrix3<double> &mat);
      inline SxMatrix3 (const SxMatrix3<SxComplex8> &mat);
      inline SxMatrix3 (const SxMatrix3<SxComplex16> &mat);
      inline SxMatrix3 (const SxVector3<T> &v1,
                        const SxVector3<T> &v2,
                        const SxVector3<T> &v3);
      inline SxMatrix3 (const T &e00, const T &e01, const T &e02,
                        const T &e10, const T &e11, const T &e12,
                        const T &e20, const T &e21, const T &e22);
      inline explicit SxMatrix3 (const SxList<T> &);

      inline void set (const T &);
      inline void setCol (ssize_t i, const SxVector3<T> &colIn);

      /** Use matrix(x,y) in order to dereference an element */
      inline SxVector3<T>        operator() (ssize_t i);
      inline const SxVector3<T>  operator() (ssize_t i) const;
      /** i and j are the row and column, resp. */
      inline T &operator() (ssize_t i, ssize_t j);
      /** i and j are the row and column, resp. */
      inline const T &operator() (ssize_t i, ssize_t j) const;
      /** i and j are the row and column, resp. */
      inline T &operator() (SxAutoLoop &i, SxAutoLoop &j);
      /** i and j are the row and column, resp. */
      inline const T &operator() (SxAutoLoop &i, SxAutoLoop &j) const;

      /** Assign a scalar value to all matrix elements */
      inline SxMatrix3<T> &operator=  (const T &s);
      /** Assign matrix to matrix */
      inline SxMatrix3<T> &operator=  (const SxMatrix3<T> &in);
      /** Multiply matrix by a scalar value */
      //inline SxMatrix3<T> operator*  (const T &s) const;
      inline SxMatrix3<T> operator*= (const T &s);
      /** \brief Add a matrix elementwise */
      SxMatrix3<T> operator+= (const SxMatrix3<T> &in);
      /** \brief Subtract a matrix elementwise */
      SxMatrix3<T> operator-= (const SxMatrix3<T> &in);
      /** Divide matrix by a scalar value */
      inline SxMatrix3<T> operator/= (const T &s);
      inline bool operator== (const SxMatrix3<T> &) const;
      inline bool operator!= (const SxMatrix3<T> &in) const;
      inline SxVector3<T> row (ssize_t i) const;
      inline SxVector3<T> col (ssize_t i) const;
      inline T sum () const;
      inline T tr  () const;
      inline SxMatrix3<typename SxTypeMapper<T>::TReal> absSqr () const;

      /** Evaluate inverse of matrix */
      inline SxMatrix3<T> inverse () const;
      /** Compute determinant of matrix */
      inline T determinant () const;
      /** Return transpose of matrix */
      inline SxMatrix3<T> transpose () const;
      /** Return trace of matrix */
      inline T trace () const;

      void print () const;
};


template<class T>
SxMatrix3<T>::SxMatrix3 ()
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = (T)0.0;
}


template<class T>
SxMatrix3<T>::SxMatrix3 (const T &val)
{
   for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
         m[i][j] = (i == j) ? val : (T)0;
}



template<class T>
SxMatrix3<T>::SxMatrix3 (const SxMatrix3<int> &mat)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = (T)mat.m[i][j];
}
template<class T>
SxMatrix3<T>::SxMatrix3 (const SxMatrix3<float> &mat)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = (T)mat.m[i][j];
}
template<class T>
SxMatrix3<T>::SxMatrix3 (const SxMatrix3<double> &mat)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = (T)mat.m[i][j];
}
template<class T>
SxMatrix3<T>::SxMatrix3 (const SxMatrix3<SxComplex8> &mat)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = (T)mat.m[i][j];
}
template<class T>
SxMatrix3<T>::SxMatrix3 (const SxMatrix3<SxComplex16> &mat)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = (T)mat.m[i][j];
}
template<>
inline SxMatrix3<int>::SxMatrix3 (const SxMatrix3<float> &mat)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = toInt(mat.m[i][j]);
}
template<>
inline SxMatrix3<int>::SxMatrix3 (const SxMatrix3<double> &mat)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = toInt(mat.m[i][j]);
}



template<class T>
SxMatrix3<T>::SxMatrix3 (const SxVector3<T> &v1,
                         const SxVector3<T> &v2,
                         const SxVector3<T> &v3)
{
    m[0][0] = v1.v[0]; m[1][0] = v2.v[0]; m[2][0] = v3.v[0];
    m[0][1] = v1.v[1]; m[1][1] = v2.v[1]; m[2][1] = v3.v[1];
    m[0][2] = v1.v[2]; m[1][2] = v2.v[2]; m[2][2] = v3.v[2];
}


template<class T>
SxMatrix3<T>::SxMatrix3 (const T &e00, const T &e01, const T &e02,
                         const T &e10, const T &e11, const T &e12,
                         const T &e20, const T &e21, const T &e22)
{
   m[0][0] = e00; m[0][1] = e01; m[0][2] = e02;
   m[1][0] = e10; m[1][1] = e11; m[1][2] = e12;
   m[2][0] = e20; m[2][1] = e21; m[2][2] = e22;
}

template<class T>
SxMatrix3<T>::SxMatrix3 (const SxList<T> &in)
{
   SX_CHECK (in.getSize () == 9, in.getSize ());
   typename SxList<T>::ConstIterator it = in.begin();
   m[0][0] = *it++; m[0][1] = *it++; m[0][2] = *it++;
   m[1][0] = *it++; m[1][1] = *it++; m[1][2] = *it++;
   m[2][0] = *it++; m[2][1] = *it++; m[2][2] = *it;
}


template<class T>
void SxMatrix3<T>::set (const T &v)
{
   m[0][0] = m[0][1] = m[0][2] =
   m[1][0] = m[1][1] = m[1][2] =
   m[2][0] = m[2][1] = m[2][2] = v;
}

template<class T>
void SxMatrix3<T>::setCol (ssize_t i, const SxVector3<T> &colIn)
{
   m[0][i] = colIn(0);
   m[1][i] = colIn(1);
   m[2][i] = colIn(2);
}


template<class T>
SxVector3<T> SxMatrix3<T>::operator() (ssize_t i)
{
   SX_CHECK (i >= 0 && i < 3, i);
   return SxVector3<T> (m[0][i], m[1][i], m[2][i]);
}

template<class T>
const SxVector3<T> SxMatrix3<T>::operator() (ssize_t i) const
{
   SX_CHECK (i >= 0 && i < 3, i);
   return SxVector3<T> (m[0][i], m[1][i], m[2][i]);
}


template<class T>
T &SxMatrix3<T>::operator() (ssize_t i, ssize_t j)
{
   SX_CHECK (i >= 0 && i < 3 && j >= 0 && j < 3, i, j);
   return m[i][j];
}


template<class T>
const T& SxMatrix3<T>::operator() (ssize_t i, ssize_t j) const
{
   SX_CHECK (i >= 0 && i < 3 && j >= 0 && j < 3, i, j);
   return m[i][j];
}

template<class T>
T &SxMatrix3<T>::operator() (SxAutoLoop& i, SxAutoLoop& j)
{
   i.setLimit (3);
   j.setLimit (3);
   return m[i.i][j.i];
}


template<class T>
const T& SxMatrix3<T>::operator() (SxAutoLoop& i, SxAutoLoop& j) const
{
   i.setLimit (3);
   j.setLimit (3);
   return m[i.i][j.i];
}



template<class T>
SxMatrix3<T> &SxMatrix3<T>::operator= (const T &s)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = s;
   return *this;
}

template<class T>
SxMatrix3<T> &SxMatrix3<T>::operator= (const SxMatrix3<T> &in)
{
   if (&in == this)  return *this;

   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = in.m[i][j];
   return *this;
}


template<class T>
SxMatrix3<T> SxMatrix3<T>::operator*= (const T &s)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] *= s;
   return *this;
}

template<class T>
SxMatrix3<T> SxMatrix3<T>::operator/= (const T &s)
{
   // doesn't work properly for SxComplex<?>:
   // SX_CHECK ( fabs ((double)s) > 1e-10 , (double)s);
   int i, j;
   T d = 1. / s;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] *= d;
   return *this;
}

template <class T>
SxMatrix3<T> SxMatrix3<T>::operator+= (const SxMatrix3<T> &in)
{
   int i,j;
   for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
         this->m[i][j] += in.m[i][j];
   return *this;
}

template <class T>
SxMatrix3<T> SxMatrix3<T>::operator-= (const SxMatrix3<T> &in)
{
   int i,j;
   for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
         this->m[i][j] -= in.m[i][j];
   return *this;
}

template<class T>
bool SxMatrix3<T>::operator!= (const SxMatrix3<T> &in) const
{
   return !operator==(in);
}

template<class T>
bool SxMatrix3<T>::operator== (const SxMatrix3<T> &in) const
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         if (m[i][j] != in.m[i][j] )   return false;
   return true;
}


template<>
inline bool SxMatrix3<float>::operator== (const SxMatrix3<float> &in) const
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         if ( fabsf(m[i][j] - in.m[i][j]) > 1e-5 )  return false;
   return true;
}


template<>
inline bool SxMatrix3<double>::operator== (const SxMatrix3<double> &in) const
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         if ( fabs(m[i][j] - in.m[i][j]) > 1e-10 )  return false;
   return true;
}



template<class T>
SxVector3<T> SxMatrix3<T>::row (ssize_t i) const
{
   SX_CHECK (i >= 0 && i < 3, i);
   return   SxVector3<T> (m[i][0], m[i][1], m[i][2]);
}


template<class T>
SxVector3<T> SxMatrix3<T>::col (ssize_t i) const
{
   SX_CHECK (i >= 0 && i < 3, i);
   return   SxVector3<T> (m[0][i], m[1][i], m[2][i]);
}




template<class T>
T SxMatrix3<T>::sum () const
{
   return   m[0][0] + m[0][1] + m[0][2]
          + m[1][0] + m[1][1] + m[1][2]
          + m[2][0] + m[2][1] + m[2][2];
}

template<class T>
T SxMatrix3<T>::tr () const
{
   return   m[0][0] + m[1][1] + m[2][2];
}


template<>
inline SxMatrix3<int> SxMatrix3<int>::absSqr () const
{
   SxMatrix3<int> res;
   int i, j;
   for (i=0; i < 3; i++)
      for (j=0; j < 3; j++)
         res.m[i][j] = m[i][j] * m[i][j];
   return res;
}

template<>
inline SxMatrix3<float> SxMatrix3<float>::absSqr () const
{
   SxMatrix3<float> res;
   int i, j;
   for (i=0; i < 3; i++)
      for (j=0; j < 3; j++)
         res.m[i][j] = m[i][j] * m[i][j];
   return res;
}

template<>
inline SxMatrix3<double> SxMatrix3<double>::absSqr () const
{
   SxMatrix3<double> res;
   int i, j;
   for (i=0; i < 3; i++)
      for (j=0; j < 3; j++)
         res.m[i][j] = m[i][j] * m[i][j];
   return res;
}

template<>
inline SxMatrix3<float> SxMatrix3<SxComplex8>::absSqr () const
{
   SxMatrix3<float> res;
   int i, j;
   for (i=0; i < 3; i++)
      for (j=0; j < 3; j++)
         res.m[i][j] = m[i][j].absSqr ();
   return res;
}

template<>
inline SxMatrix3<double> SxMatrix3<SxComplex16>::absSqr () const
{
   SxMatrix3<double> res;
   int i, j;
   for (i=0; i < 3; i++)
      for (j=0; j < 3; j++)
         res.m[i][j] = m[i][j].absSqr ();
   return res;
}


template<class T>
T SxMatrix3<T>::determinant() const
{
  return (m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2])+
          m[1][0] * (m[2][1] * m[0][2] - m[0][1] * m[2][2])+
	       m[2][0] * (m[0][1] * m[1][2] - m[1][1] * m[0][2]));
}

template<class T>
T SxMatrix3<T>::trace () const
{
   return (m[0][0] + m[1][1] + m[2][2]);
}


template<class T>
SxMatrix3<T> SxMatrix3<T>::transpose () const
{
   SxMatrix3<T> trans;
   int i,j;
   for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
         trans.m[i][j] = m[j][i];
   
   return trans;
}


template<class T>
SxMatrix3<T> SxMatrix3<T>::inverse () const
{
  const SxMatrix3<T> &M = *this;

  SxMatrix3<T> inv (M(1).x (M(2)), 
                    M(2).x (M(0)), 
                    M(0).x (M(1)));

  SX_CHECK (fabs((double)determinant()) > 1e-10, (double)determinant());
  inv *= ((double)1.) / determinant ();
  return inv;
}

template<class T>
void SxMatrix3<T>::print () const
{
   sxprintf ("[[%f, %f, %f]\n",  (float)m[0][0], (float)m[0][1], (float)m[0][2]);
   sxprintf (" [%f, %f, %f]\n",  (float)m[1][0], (float)m[1][1], (float)m[1][2]);
   sxprintf (" [%f, %f, %f]]\n", (float)m[2][0], (float)m[2][1], (float)m[2][2]);
}

// --- The next 220 lines were generated from ./snippets/SxMatrix3.h snippet ScalarOpMat
// SxMatrix3<?> = int * SxMatrix3<A>
template<class B>
inline SxMatrix3<decltype(int(0)*B(0))>
operator* (int a, const SxMatrix3<B> &b)
{
   return SxMatrix3<decltype(int(0)*B(0))>
      (a * b.m[0][0], a * b.m[0][1], a * b.m[0][2],
       a * b.m[1][0], a * b.m[1][1], a * b.m[1][2],
       a * b.m[2][0], a * b.m[2][1], a * b.m[2][2]);
}

// SxMatrix3<?> = SxMatrix3<A> * int
template<class A>
inline SxMatrix3<decltype(A(0)*int(0))>
operator* (const SxMatrix3<A> &a, int b)
{
   return SxMatrix3<decltype(A(0)*int(0))>
      (a.m[0][0] * b, a.m[0][1] * b, a.m[0][2] * b,
       a.m[1][0] * b, a.m[1][1] * b, a.m[1][2] * b,
       a.m[2][0] * b, a.m[2][1] * b, a.m[2][2] * b);
}

// SxMatrix3<?> = float * SxMatrix3<A>
template<class B>
inline SxMatrix3<decltype(float(0)*B(0))>
operator* (float a, const SxMatrix3<B> &b)
{
   return SxMatrix3<decltype(float(0)*B(0))>
      (a * b.m[0][0], a * b.m[0][1], a * b.m[0][2],
       a * b.m[1][0], a * b.m[1][1], a * b.m[1][2],
       a * b.m[2][0], a * b.m[2][1], a * b.m[2][2]);
}

// SxMatrix3<?> = SxMatrix3<A> * float
template<class A>
inline SxMatrix3<decltype(A(0)*float(0))>
operator* (const SxMatrix3<A> &a, float b)
{
   return SxMatrix3<decltype(A(0)*float(0))>
      (a.m[0][0] * b, a.m[0][1] * b, a.m[0][2] * b,
       a.m[1][0] * b, a.m[1][1] * b, a.m[1][2] * b,
       a.m[2][0] * b, a.m[2][1] * b, a.m[2][2] * b);
}

// SxMatrix3<?> = double * SxMatrix3<A>
template<class B>
inline SxMatrix3<decltype(double(0)*B(0))>
operator* (double a, const SxMatrix3<B> &b)
{
   return SxMatrix3<decltype(double(0)*B(0))>
      (a * b.m[0][0], a * b.m[0][1], a * b.m[0][2],
       a * b.m[1][0], a * b.m[1][1], a * b.m[1][2],
       a * b.m[2][0], a * b.m[2][1], a * b.m[2][2]);
}

// SxMatrix3<?> = SxMatrix3<A> * double
template<class A>
inline SxMatrix3<decltype(A(0)*double(0))>
operator* (const SxMatrix3<A> &a, double b)
{
   return SxMatrix3<decltype(A(0)*double(0))>
      (a.m[0][0] * b, a.m[0][1] * b, a.m[0][2] * b,
       a.m[1][0] * b, a.m[1][1] * b, a.m[1][2] * b,
       a.m[2][0] * b, a.m[2][1] * b, a.m[2][2] * b);
}

// SxMatrix3<?> = SxComplex8 * SxMatrix3<A>
template<class B>
inline SxMatrix3<decltype(SxComplex8(0)*B(0))>
operator* (const SxComplex8 &a, const SxMatrix3<B> &b)
{
   return SxMatrix3<decltype(SxComplex8(0)*B(0))>
      (a * b.m[0][0], a * b.m[0][1], a * b.m[0][2],
       a * b.m[1][0], a * b.m[1][1], a * b.m[1][2],
       a * b.m[2][0], a * b.m[2][1], a * b.m[2][2]);
}

// SxMatrix3<?> = SxMatrix3<A> * SxComplex8
template<class A>
inline SxMatrix3<decltype(A(0)*SxComplex8(0))>
operator* (const SxMatrix3<A> &a, const SxComplex8 &b)
{
   return SxMatrix3<decltype(A(0)*SxComplex8(0))>
      (a.m[0][0] * b, a.m[0][1] * b, a.m[0][2] * b,
       a.m[1][0] * b, a.m[1][1] * b, a.m[1][2] * b,
       a.m[2][0] * b, a.m[2][1] * b, a.m[2][2] * b);
}

// SxMatrix3<?> = SxComplex16 * SxMatrix3<A>
template<class B>
inline SxMatrix3<decltype(SxComplex16(0)*B(0))>
operator* (const SxComplex16 &a, const SxMatrix3<B> &b)
{
   return SxMatrix3<decltype(SxComplex16(0)*B(0))>
      (a * b.m[0][0], a * b.m[0][1], a * b.m[0][2],
       a * b.m[1][0], a * b.m[1][1], a * b.m[1][2],
       a * b.m[2][0], a * b.m[2][1], a * b.m[2][2]);
}

// SxMatrix3<?> = SxMatrix3<A> * SxComplex16
template<class A>
inline SxMatrix3<decltype(A(0)*SxComplex16(0))>
operator* (const SxMatrix3<A> &a, const SxComplex16 &b)
{
   return SxMatrix3<decltype(A(0)*SxComplex16(0))>
      (a.m[0][0] * b, a.m[0][1] * b, a.m[0][2] * b,
       a.m[1][0] * b, a.m[1][1] * b, a.m[1][2] * b,
       a.m[2][0] * b, a.m[2][1] * b, a.m[2][2] * b);
}

// SxMatrix3<?> = int / SxMatrix3<A>
template<class B>
inline SxMatrix3<decltype(int(0)/B(0))>
operator/ (int a, const SxMatrix3<B> &b)
{
   return SxMatrix3<decltype(int(0)/B(0))>
      (a / b.m[0][0], a / b.m[0][1], a / b.m[0][2],
       a / b.m[1][0], a / b.m[1][1], a / b.m[1][2],
       a / b.m[2][0], a / b.m[2][1], a / b.m[2][2]);
}

// SxMatrix3<?> = SxMatrix3<A> / int
template<class A>
inline SxMatrix3<decltype(A(0)/int(0))>
operator/ (const SxMatrix3<A> &a, int b)
{
   return SxMatrix3<decltype(A(0)/int(0))>
      (a.m[0][0] / b, a.m[0][1] / b, a.m[0][2] / b,
       a.m[1][0] / b, a.m[1][1] / b, a.m[1][2] / b,
       a.m[2][0] / b, a.m[2][1] / b, a.m[2][2] / b);
}

// SxMatrix3<?> = float / SxMatrix3<A>
template<class B>
inline SxMatrix3<decltype(float(0)/B(0))>
operator/ (float a, const SxMatrix3<B> &b)
{
   return SxMatrix3<decltype(float(0)/B(0))>
      (a / b.m[0][0], a / b.m[0][1], a / b.m[0][2],
       a / b.m[1][0], a / b.m[1][1], a / b.m[1][2],
       a / b.m[2][0], a / b.m[2][1], a / b.m[2][2]);
}

// SxMatrix3<?> = SxMatrix3<A> / float
template<class A>
inline SxMatrix3<decltype(A(0)/float(0))>
operator/ (const SxMatrix3<A> &a, float b)
{
   return SxMatrix3<decltype(A(0)/float(0))>
      (a.m[0][0] / b, a.m[0][1] / b, a.m[0][2] / b,
       a.m[1][0] / b, a.m[1][1] / b, a.m[1][2] / b,
       a.m[2][0] / b, a.m[2][1] / b, a.m[2][2] / b);
}

// SxMatrix3<?> = double / SxMatrix3<A>
template<class B>
inline SxMatrix3<decltype(double(0)/B(0))>
operator/ (double a, const SxMatrix3<B> &b)
{
   return SxMatrix3<decltype(double(0)/B(0))>
      (a / b.m[0][0], a / b.m[0][1], a / b.m[0][2],
       a / b.m[1][0], a / b.m[1][1], a / b.m[1][2],
       a / b.m[2][0], a / b.m[2][1], a / b.m[2][2]);
}

// SxMatrix3<?> = SxMatrix3<A> / double
template<class A>
inline SxMatrix3<decltype(A(0)/double(0))>
operator/ (const SxMatrix3<A> &a, double b)
{
   return SxMatrix3<decltype(A(0)/double(0))>
      (a.m[0][0] / b, a.m[0][1] / b, a.m[0][2] / b,
       a.m[1][0] / b, a.m[1][1] / b, a.m[1][2] / b,
       a.m[2][0] / b, a.m[2][1] / b, a.m[2][2] / b);
}

// SxMatrix3<?> = SxComplex8 / SxMatrix3<A>
template<class B>
inline SxMatrix3<decltype(SxComplex8(0)/B(0))>
operator/ (const SxComplex8 &a, const SxMatrix3<B> &b)
{
   return SxMatrix3<decltype(SxComplex8(0)/B(0))>
      (a / b.m[0][0], a / b.m[0][1], a / b.m[0][2],
       a / b.m[1][0], a / b.m[1][1], a / b.m[1][2],
       a / b.m[2][0], a / b.m[2][1], a / b.m[2][2]);
}

// SxMatrix3<?> = SxMatrix3<A> / SxComplex8
template<class A>
inline SxMatrix3<decltype(A(0)/SxComplex8(0))>
operator/ (const SxMatrix3<A> &a, const SxComplex8 &b)
{
   return SxMatrix3<decltype(A(0)/SxComplex8(0))>
      (a.m[0][0] / b, a.m[0][1] / b, a.m[0][2] / b,
       a.m[1][0] / b, a.m[1][1] / b, a.m[1][2] / b,
       a.m[2][0] / b, a.m[2][1] / b, a.m[2][2] / b);
}

// SxMatrix3<?> = SxComplex16 / SxMatrix3<A>
template<class B>
inline SxMatrix3<decltype(SxComplex16(0)/B(0))>
operator/ (const SxComplex16 &a, const SxMatrix3<B> &b)
{
   return SxMatrix3<decltype(SxComplex16(0)/B(0))>
      (a / b.m[0][0], a / b.m[0][1], a / b.m[0][2],
       a / b.m[1][0], a / b.m[1][1], a / b.m[1][2],
       a / b.m[2][0], a / b.m[2][1], a / b.m[2][2]);
}

// SxMatrix3<?> = SxMatrix3<A> / SxComplex16
template<class A>
inline SxMatrix3<decltype(A(0)/SxComplex16(0))>
operator/ (const SxMatrix3<A> &a, const SxComplex16 &b)
{
   return SxMatrix3<decltype(A(0)/SxComplex16(0))>
      (a.m[0][0] / b, a.m[0][1] / b, a.m[0][2] / b,
       a.m[1][0] / b, a.m[1][1] / b, a.m[1][2] / b,
       a.m[2][0] / b, a.m[2][1] / b, a.m[2][2] / b);
}
// --- ScalarOpMat

// SxVector3<?> = SxMatrix3<A> ^ SxVector3<B>
template<class A,class B>
inline SxVector3<decltype(A(0)*B(0))>
operator^ (const SxMatrix3<A> &a, const SxVector3<B> &b)
{
   return SxVector3<decltype(A(0)*B(0))>
      (a.m[0][0] * b.v[0] + a.m[0][1] * b.v[1] + a.m[0][2] * b.v[2],
       a.m[1][0] * b.v[0] + a.m[1][1] * b.v[1] + a.m[1][2] * b.v[2],
       a.m[2][0] * b.v[0] + a.m[2][1] * b.v[1] + a.m[2][2] * b.v[2]);
}

// SxVector3<?> = SxVector3<A> ^ SxMatrix3<B>
template<class A,class B>
inline SxVector3<decltype(A(0)*B(0))>
operator^ (const SxVector3<A> &a, const SxMatrix3<B> &b)
{
   return SxVector3<decltype(A(0)*B(0))>
         (b.m[0][0] * a.v[0] + b.m[1][0] * a.v[1] + b.m[2][0] * a.v[2],
          b.m[0][1] * a.v[0] + b.m[1][1] * a.v[1] + b.m[2][1] * a.v[2],
          b.m[0][2] * a.v[0] + b.m[1][2] * a.v[1] + b.m[2][2] * a.v[2]);
}

// SxMatrix3<?> = SxMatrix3<A> ^ SxMatrix3<B>
template<class A,class B>
inline SxMatrix3<decltype(A(0)*B(0))>
operator^ (const SxMatrix3<A> &a, const SxMatrix3<B> &b)
{
   SxMatrix3<decltype(A(0)*B(0))> res;
   int r, c;
   for (c=0; c < 3; c++)
      for (r=0; r < 3; r++)  
         res.m[r][c] = a.m[r][0] * b.m[0][c]
                     + a.m[r][1] * b.m[1][c]
                     + a.m[r][2] * b.m[2][c];
   return res;
}

// --- The next 44 lines were generated from ./snippets/SxMatrix3.h snippet MatOpMat
// SxMatrix3<?> = SxMatrix3<A> + SxMatrix3<B>
template<class A,class B>
inline SxMatrix3<decltype(A(0)+B(0))>
operator+ (const SxMatrix3<A> &a, const SxMatrix3<B> &b)
{
   return SxMatrix3<decltype(A(0)+B(0))>
      (a.m[0][0] + b.m[0][0], a.m[0][1] + b.m[0][1], a.m[0][2] + b.m[0][2],
       a.m[1][0] + b.m[1][0], a.m[1][1] + b.m[1][1], a.m[1][2] + b.m[1][2],
       a.m[2][0] + b.m[2][0], a.m[2][1] + b.m[2][1], a.m[2][2] + b.m[2][2]);
}

// SxMatrix3<?> = SxMatrix3<A> - SxMatrix3<B>
template<class A,class B>
inline SxMatrix3<decltype(A(0)-B(0))>
operator- (const SxMatrix3<A> &a, const SxMatrix3<B> &b)
{
   return SxMatrix3<decltype(A(0)-B(0))>
      (a.m[0][0] - b.m[0][0], a.m[0][1] - b.m[0][1], a.m[0][2] - b.m[0][2],
       a.m[1][0] - b.m[1][0], a.m[1][1] - b.m[1][1], a.m[1][2] - b.m[1][2],
       a.m[2][0] - b.m[2][0], a.m[2][1] - b.m[2][1], a.m[2][2] - b.m[2][2]);
}

// SxMatrix3<?> = SxMatrix3<A> * SxMatrix3<B>
template<class A,class B>
inline SxMatrix3<decltype(A(0)*B(0))>
operator* (const SxMatrix3<A> &a, const SxMatrix3<B> &b)
{
   return SxMatrix3<decltype(A(0)*B(0))>
      (a.m[0][0] * b.m[0][0], a.m[0][1] * b.m[0][1], a.m[0][2] * b.m[0][2],
       a.m[1][0] * b.m[1][0], a.m[1][1] * b.m[1][1], a.m[1][2] * b.m[1][2],
       a.m[2][0] * b.m[2][0], a.m[2][1] * b.m[2][1], a.m[2][2] * b.m[2][2]);
}

// SxMatrix3<?> = SxMatrix3<A> / SxMatrix3<B>
template<class A,class B>
inline SxMatrix3<decltype(A(0)/B(0))>
operator/ (const SxMatrix3<A> &a, const SxMatrix3<B> &b)
{
   return SxMatrix3<decltype(A(0)/B(0))>
      (a.m[0][0] / b.m[0][0], a.m[0][1] / b.m[0][1], a.m[0][2] / b.m[0][2],
       a.m[1][0] / b.m[1][0], a.m[1][1] / b.m[1][1], a.m[1][2] / b.m[1][2],
       a.m[2][0] / b.m[2][0], a.m[2][1] / b.m[2][1], a.m[2][2] / b.m[2][2]);
}
// --- MatOpMat

//------------------------------------------------------------------------------
template<class T>
std::ostream& operator<< (std::ostream &s, const SxMatrix3<T> &in)
{
   return s << "{{" << in.m[0][0] << "," << in.m[0][1] << "," << in.m[0][2] << "},"
            << " {" << in.m[1][0] << "," << in.m[1][1] << "," << in.m[1][2] << "},"
            << " {" << in.m[2][0] << "," << in.m[2][1] << "," << in.m[2][2] << "}}";
}



#endif // _SX_MATRIX3_H_
