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

#ifndef _CUBICSPLINE_H_
#define _CUBICSPLINE_H_

// Including the header-files
// -- SPHInX-Headers
#include <SxArray.h>
#include <SxMath.h>
#include <SxError.h>
#include <SxVector.h>

class SX_EXPORT_MATH SxCubicSpline
{

   public:

      //enum WorkingMode {Spline, LeastSquareFit};
      enum SplineType {Hermite, Natural, NaturalHermite, HermiteNatural};
      enum FitType {None, Normal, Extend, MirrorPlane, MirrorPoint};
      // Constructors
      // -- standard constructor
      SxCubicSpline () { /* empty */ }

      /// Constructor for read spline
      SxCubicSpline (const SxVecRef<double> &x,
                     const SxVecRef<double> &spline);

      /// Constructor for x-y interpolation
      SxCubicSpline (const SxVecRef<double> &x,
                     const SxVecRef<double> &y,
                     const enum SplineType splineType,
                     const double slope1 = 0.0,
                     const double slope2 = 0.0);

      // Mixed Constructure due to same variable mask
      // When using Spline mode:
      //   vec1 = x, vec2 = dx, and vec3 = y
      // When using LeastSquareFit mode:
      //   vec1 = xData, vec2 = yData, and vec3 = basis
      SxCubicSpline (const SxVecRef<double> &x,
                     const SxVecRef<double> &dx,
                     const SxVecRef<double> &y,
                     const enum SplineType splineType,
                     const double slope1 = 0.0,
                     const double slope2 = 0.0);

      SxCubicSpline (const SxVecRef<double> &xData,
                     const SxVecRef<double> &yData,
                     const SxVecRef<double> &basis,
                     const enum SplineType splineType,
                     const enum FitType fitType,
                     const double slope1 = 0.0,
                     const double slope2 = 0.0);

      SxVector<double> getYFit () const;

      inline double getY (const double x) const;

      SxVector<double> getY (const SxVecRef<double> &x) const;

      inline double getdYdX (const double x) const;

      SxVector<double> getdYdX (const SxVecRef<double> &x) const;

      void setSpline (const SxVecRef<double> &spline);

      SxVector<double> getSpline () const;

      // calculate Splinepoint dependencys
      SxArray<int> getSplineDep (const enum FitType fitType);

      // set of meta data for dirac Vectors ! have to be set by hand
      int iSpecies;
      int iAtom;
      int n;
      int l;
      int m;

      inline void setXFit (const SxVecRef<double> &in) {xFit = in;};
      inline void setYFit (const SxVecRef<double> &in) {yFit = in;};
      inline void setXVals (const SxVecRef<double> &in) {xVals = in;};
      inline void setYVals (const SxVecRef<double> &in) {yVals = in;};
      inline void setSplineType (const enum SplineType in) {splineType = in;};

      SxVector<double> calcDRDataDRGrid ();


  protected:

      // Methods
      // Cubic Spline Interpolation
      void computeSpline ();
      // Least Square Fit
      void computeLeastSquareFit ();
      // Setup difference in xVals
      void setH ();
      // Setup coeff Matrix for Splinecoeff calculation
      void setT (const enum SplineType sType);
      // Setup right hand side for spline coeff calculation
      void setBeta (const enum SplineType sType);
      // get Spline Index
      inline int getSplineIdx(const double x,
                              const SxVecRef<double>& basis) const;
      /// Non-inline crash from getSplineIdx
      static void errorIdx (double x, const SxVecRef<double> &basis, int idx);
      // mirror
      void mirror ();
      // demirror
      void deMirror ();

      SxVector<double> symmetricTridiagonalGauss (const SxVecRef<double> &Mat, const SxVecRef<double> &rhs);

      void print (const SxVecRef<double> &diag1, const SxVecRef<double> &diag2, const SxVecRef<double> &diag3, const SxVecRef<double> &b);

      SxArray<SxList<int> > getSplineExtend (SxArray<int> &pointsPerSpline);
      SxArray<int> getPointsPerSpline (const SxVecRef<double>& basis, const SxVecRef<double>& xData);
      SxVector<double> cleanZeroDepths(const SxVecRef<double>& xData, const SxArray<int> &pointsPerSpline);
      SxVector<double> dady ();
      SxVector<double> dbdy (const SxVecRef<double> &dA, const SxVecRef<double> &dC);
      SxVector<double> dcdy (const enum SplineType sType);
      SxVector<double> dbetady (const enum SplineType sType);
      SxVector<double> dddy (const SxVecRef<double> &dC);

      // Members
      SxVector<double> polyCoeff;
      SxVector<double> xVals;
      SxVector<double> yVals;
      SxVector<double> hVals;
      SxVector<double> TMat;
      SxVector<double> betaVals;
      SxVector<double> xFit;
      SxVector<double> yFit;
      SxArray<SxList<int> > depArray;
      enum SplineType splineType;
      enum FitType fitType;
      double leftHermite;
      double rightHermite;
};

inline double SxCubicSpline::getY (double x) const
{
  SX_CHECK(polyCoeff.getSize () > 0);
  SX_CHECK(polyCoeff.getSize () % 4 == 0);
  int iSpline = getSplineIdx(x, xVals);
  double diff = x - xVals(iSpline);
  SX_CHECK(polyCoeff.getSize () > 4*iSpline + 3, iSpline, polyCoeff.getSize ());
  double result = polyCoeff(4*iSpline + 0) + diff
              * ( polyCoeff(4*iSpline + 1) + diff
              * ( polyCoeff(4*iSpline + 2) + diff
              *   polyCoeff(4*iSpline + 3)));
  SX_CHECK_NUM(result);
  return result;
}

inline double SxCubicSpline::getdYdX (double x) const
{
  SX_CHECK(polyCoeff.getSize () > 0);
  SX_CHECK(polyCoeff.getSize () % 4 == 0);
  int iSpline = getSplineIdx(x, xVals);
  double diff = x - xVals(iSpline);
  SX_CHECK(polyCoeff.getSize () > 4*iSpline + 3, iSpline, polyCoeff.getSize ());
  double result = polyCoeff(4*iSpline + 1) + diff
              * ( 2.0 * polyCoeff(4*iSpline + 2) + diff
              *   3.0 * polyCoeff(4*iSpline + 3));
  SX_CHECK_NUM(result);
  return result;
}

inline int
SxCubicSpline::getSplineIdx (double x, const SxVecRef<double>& basis) const
{
   SX_CHECK (basis.getSize () > 0);
   int nSpline = (int)basis.getSize ();
   // Check lower bound
   if (fabs(x-basis(0)) < 1e-10 || (x > basis(0) && x < basis(1)) )
      return 0;
   else if (x < basis(0))   {
      /*
      cout << SX_SEPARATOR;
      cout << "Warning!" << endl;
      cout << "Here is SxCubicSpline!" << endl;
      cout << "Value lies beyound left grid point!" << endl;
      cout << "Val: " << x << ", Netpoint: " << basis(0) << endl;
      cout << SX_SEPARATOR;
      */
      return 0;
   }
   // Check upper bound
   if (fabs(x-basis(nSpline-1)) < 1e-10 || (x < basis(nSpline-1) && x > basis(nSpline-2)))
      return nSpline - 2;
   else if (x > basis(nSpline - 1))   {
      /*
      cout << SX_SEPARATOR;
      cout << "Warning !" << endl;
      cout << "Here is SxCubicSpline!" << endl;
      cout << "Value lies beyound right grid point!" << endl;
      cout << "Val: " << x << ", Netpoint: " << basis(nSpline - 1) << endl;
      cout << SX_SEPARATOR;
      */
      return nSpline - 2;
   }
   // intervall is half of nSpline when nSpline is even or [nSpline/2] when nSpline is odd
   int interval = nSpline / 2;
   int result = interval;
   int counter = nSpline;
   // find Idx via interval intersection
   while (counter--) {
      if (result + 1 >= nSpline) break;
      if (x >= basis(result) && x <= basis(result+1)) return result;
      interval = interval / 2;
      if (interval == 0) interval = 1;
      if (x < basis(result)) result -= interval;
      else result += interval;
   }
   errorIdx (x, basis, result);
   return -1;
}



#endif /* _CUBICSPLINE_H_ */
