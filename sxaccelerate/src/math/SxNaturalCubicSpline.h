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

// Headerguard
#ifndef _NATURALCUBICSPLINE_H_
#define _NATURALCUBICSPLINE_H_

// Including the header-files
// -- SPHInX-Headers
#include <SxVector.h>
#include <SxArray.h>
#include <SxPtr.h>
#include <SxString.h>
#include <SxMath.h>

// -- C++-Standardheaders
#include <iostream>
#include <cmath>
using namespace std;


/** \brief Natural cubic spline interpolation

    \b SxClass = SFHIngX natural cubic spline interpolation

    \author T. Uchdorf, uchdorf@mpie.de
    \author Christoph Freysoldt, freysoldt@mpie.de
 */

class SX_EXPORT_MATH SxNaturalCubicSpline
{
   public:
      // Constructors
      // -- Constructor
      SxNaturalCubicSpline ();
      
      /// Constructor for x-y interpolation
      SxNaturalCubicSpline (const SxVector<double> &x, const SxVector<double> &y);
      /// Constructor for y interpolation (index serves as x)
      SxNaturalCubicSpline (const SxVector<double> &y, bool estimateD2 = false);

      void compute (const SxVector<double> &y, bool estimateD2 = false);

      // -- Destructor
      ~SxNaturalCubicSpline ();

      inline double getVal (double x) const
      {
         return (xVals.getSize () == 0) ? getValY (x) : getValXY (x);
      }

      SxVector<double> getVals (const SxVecRef<double> &x) const;

      /// Get value for i->y(i) even outside original data range
      double getValYExtra (double x) const;
      /// Get derivatives for i->y(i) even outside original data range
      double getDerivYExtra (double x) const;

  protected:

      // Methods
      void print (
            const SxVector<double> &diag1,
            const SxVector<double> &diag2,
            const SxVector<double> &diag3,
            const SxVector<double> &b);

      void compute (const SxVector<double> &x,const SxVector<double> &y);
      
      SxVector<double> symmetricTridiagonalGauss (
            const SxVector<double> & upper,
            const SxVector<double> &mid,
            const SxVector<double> & lower,
            const SxVector<double> &rhs);

      // Members
      SxArray<SxVector<double> > polyCoeff;

      SxVector<double> xVals;

      /// Get value for i->y(i) interpolation 
      double getValY (double x) const;
      /// Get value for x(i)->y(i) interpolation 
      double getValXY (double x) const;
  public:
      const SxArray<SxVector<double> > & getCoeff () const
      {
         return polyCoeff;
      }
};

inline double SxNaturalCubicSpline::getValY (double x) const
{
   int i = int(floor(x));
   if (i == -1 && fabs(x) < 1e-12) return polyCoeff(0)(0);
   const SxVector<double> &data = polyCoeff(0);
   if (i == data.getSize () - 1 && fabs(x - i) < 1e-12) return data(i);
   SX_CHECK (i >= 0 && i < data.getSize () - 1, i, data.getSize ());
   double t = x -i, u = 1. - t;
   SX_CHECK (t >= -1e-12 && t <= 1.0000000000001, t);
   SX_CHECK (u >= -1e-12 && u <= 1.0000000000001, u);
   const SxVector<double> &z = polyCoeff(1);
   // S= (z_{i+1}(x-x_i)^3 + z_i(x_{i+1}-x)^3) / (6 h_i)
   //    + (y_{i+1}/h_i - h_i*z_{i+1}/6)(x-x_i)
   //    + (y_{i}/h_i - h_i*z_{i}/6)(x_{i+1}-x)
   double s = (data(i+1) + z(i+1)*(t*t - 1.)) * t
            + (data(i) + z(i)*(u*u - 1.)) * u;
   return s;
}

inline double SxNaturalCubicSpline::getValYExtra (double x) const
{
   SX_CHECK (xVals.getSize () == 0);
   SX_CHECK (polyCoeff.getSize () == 2, polyCoeff.getSize ());
   int i = int(floor(x));
   const SxVector<double> &data = polyCoeff(0);
   if (i < 0) i = 0;
   else if (i >= data.getSize () - 1) i = static_cast<int>(data.getSize ()) - 2;
   double t = x -i, u = 1. - t;
   const SxVector<double> &z = polyCoeff(1);
   // S= (z_{i+1}(x-x_i)^3 + z_i(x_{i+1}-x)^3) / (6 h_i)
   //    + (y_{i+1}/h_i - h_i*z_{i+1}/6)(x-x_i)
   //    + (y_{i}/h_i - h_i*z_{i}/6)(x_{i+1}-x)
   double s = (data(i+1) + z(i+1)*(t*t - 1.)) * t
            + (data(i) + z(i)*(u*u - 1.)) * u;
   return s;
}

inline double SxNaturalCubicSpline::getDerivYExtra (double x) const
{
   SX_CHECK (xVals.getSize () == 0);
   SX_CHECK (polyCoeff.getSize () == 2, polyCoeff.getSize ());
   int i = int(floor(x));
   const SxVector<double> &data = polyCoeff(0);
   if (i < 0) i = 0;
   else if (i >= data.getSize () - 1) i = static_cast<int>(data.getSize ()) - 2;
   double t = x -i, u = 1. - t;
   const SxVector<double> &z = polyCoeff(1);
   // S'= (3z_{i+1}(x-x_i)^2 - 3z_i(x_{i+1}-x)^2) / (6 h_i)
   //    + (y_{i+1}/h_i - h_i*z_{i+1}/6)
   //    - (y_{i}/h_i - h_i*z_{i}/6)
   double s = (data(i+1) + z(i+1)*(3. * t*t - 1.))
            - (data(i)   + z(i)  *(3. * u*u - 1.));
   return s;
}

#endif /* _NATURALCUBICSPLINE_H_ */
