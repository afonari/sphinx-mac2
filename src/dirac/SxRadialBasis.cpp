// ---------------------------------------------------------------------------
//
//      The ab-initio based multiscale library
//
//                  S / P H I / n X
//
//      Copyright:  Max-Planck-Institut fuer Eisenforschung GmbH
//                  40237 Duesseldorf, Germany
//
//      Contact:    https://sxrepo.mpie.de
//      Authors:    see sphinx/AUTHORS
//      License:    see sphinx/LICENSE
//
// ---------------------------------------------------------------------------

#include <SxRadialBasis.h>
#include <SxRadBasis.h>
#include <SxGBasis.h>
#include <SxYlm.h>
#include <SxCubicSpline.h>


SxRadialBasis::SxRadialBasis ()
{
   //empty
}

SxRadialBasis::~SxRadialBasis ()
{
   //empty/
}

SxRadialBasis::SxRadialBasis (double radMin, double radMax, int nPoints, bool realSpace_, enum gridType mode)
   : realSpace(realSpace_)
{
   set (radMin, radMax, nPoints, mode);
}

SxRadialBasis::SxRadialBasis (const SxVector<double> &basis, bool realSpace_)
   : realSpace(realSpace_)
{
   set (basis);
}

SxRadialBasis::SxRadialBasis (const SxVecRef<CoeffType> &dRadIn, double factor, double rad0, bool realSpace_)
   : realSpace(realSpace_)
{
   set (dRadIn, factor, rad0);
}

void SxRadialBasis::set (double radMin, double radMax, int nPoints, enum gridType mode)
{
   rad.resize(nPoints);
   dRad.resize(nPoints);
   if (mode == Linear)   {
      // linear grid
      double deltaR = (radMax - radMin) / double (nPoints - 1);
      dRad.set(deltaR);
      for(int i = 0; i < nPoints; i++)   {
         rad(i) = radMin + double(i) * dRad(i);
      }
   } else if (mode == Quadratic)   {
      // quadratic grid
      double deltaR = (radMax - radMin) / sqrt(1.0*(nPoints - 1.));
      for(int i = 0; i < nPoints; i++)   {
         rad(i) = radMin + sqrt(1.0*i) * deltaR;
         if(i > 0) dRad(i-1) = rad(i) - rad(i-1);
      }
   } else if (mode == Mixed)   {
      // quadratic grid + linear grid    1:2
      double lenght = (radMax - radMin);
      double quadLenght = lenght / 3.0;
      double deltaRQuad = quadLenght / sqrt(nPoints/3.0 - 1.);
      double deltaRLin = lenght / (nPoints - 1);
      for(int i = 0; i < nPoints; i++)   {
         if (i < nPoints/3.) rad(i) = radMin + sqrt(1.0*i) * deltaRQuad;
         else rad(i) = radMin + i * deltaRLin;
         if(i > 0) dRad(i-1) = rad(i) - rad(i-1);
      }
   } else if (mode == Hyperbel)   {
      // hyperbel grid
      double diff = radMax - radMin;
      double b = (nPoints - 1.) / (2.*diff);
      double a = 1. / (2.*diff);
      for(int i = 0; i < nPoints; i++)   {
         rad(i) = radMin + 1.*i / (b+a*i);
         if(i > 0) dRad(i-1) = rad(i) - rad(i-1);
      }
   } else if (mode == Logarithmic)  {
      SX_CHECK (fabs(radMin) > 1e-12, radMin);
      double logDr = (log(radMax) - log(radMin)) / (nPoints - 1.0);
      for(int i = 0; i < nPoints; i++)   {
         rad(i) = radMin * exp(i * logDr);
         if(i > 0) dRad(i-1) = rad(i) - rad(i-1);
      }
   } else {
      cout << "UNKNOWN MODE IN SXRADGBASIS!" << endl;
      SX_EXIT;
   }

   // last dRad is equal zero conventionally
   dRad(nPoints-1) = 0.0;

   SX_VALIDATE_VECTOR(rad);
   SX_VALIDATE_VECTOR(dRad);
}

void SxRadialBasis::set (const SxVector<CoeffType> &basis)
{
   rad = basis;
   int nPoints = (int)rad.getSize();

   dRad.resize(nPoints);
   dRad(nPoints-1) = 0.0;
   for(int i = 1; i < nPoints; i++)   {
      dRad(i-1) = rad(i) - rad(i-1);
   }

   SX_VALIDATE_VECTOR(dRad);
}

void SxRadialBasis::set (const SxVecRef<CoeffType> &dRadIn, double factor, double rad0)
{
   int nPoints = (int)dRadIn.getSize();
   rad.resize(nPoints);
   rad(0) = rad0;
   dRad.resize(nPoints);
   dRad(nPoints - 1) = 0.0;
   for(int i = 1; i < nPoints; i++)   {
      if (dRadIn(i-1) > 1e-8)
         dRad(i-1) = factor / dRadIn (i-1) /(1.0 * nPoints);
      else {
       cout << "i = " << i << " is " << dRad(i-1) << endl;
         SX_EXIT;
      }
      rad(i) = rad(i-1) + dRad(i-1);
   }
   dRad.print();
   SX_VALIDATE_VECTOR(rad);
}

SxVecRef<SxRadialBasis::CoeffType> SxRadialBasis::radToRad(const SxRadialBasis *basisPtr,
      const SxVecRef<double> &vec) const
{
   if (basisPtr == this)  {
      // --- verify that this is really the identity projector, i.e.
      //     that the provided basis is the same as the vector's basis
      SX_CHECK (vec.getBasisPtr() == this);
      return const_cast<SxVecRef<SxRadialBasis::CoeffType>&>(vec);
   }
   if (realSpace && !basisPtr->realSpace)
      return fourier (basisPtr, vec);
   if (!realSpace && basisPtr->realSpace)
      return fourier (basisPtr, vec);
   // r->r or g->g interpolation is not implemented
   SX_EXIT;
}

SxVector<double> SxRadialBasis::toRadBasis (const SxRadBasis *radBasisPtr,
                                            const SxVecRef<double> &vec) const
{
   SX_CHECK(radBasisPtr);
   if (!realSpace)  {
      /// support transformation of pre-splined vectors
      if (vec.getSize () == 4 * getNElements ())  {
         SxVector<double> y = vec.getRef<Strided> (0, getNElements (), 1, 4);
         y.auxData = vec.auxData;
         return toRadBasis (radBasisPtr, y);
      }
      SX_CLOCK (Timer::RadG2Rad);

      int iSpecies = vec.auxData.is;
      int iAtom    = vec.auxData.ia;
      int n        = vec.auxData.n;
      int l        = vec.auxData.l;
      int m        = vec.auxData.m;
      const SxVector<double> &rAbs = radBasisPtr->radFunc(iSpecies);
      int dimRad = (int)rAbs.getSize();
      SxVector<double> result (dimRad);
      for (int ir = 0; ir < dimRad; ir++)   {
         SxVector<double> jl = jsb(l,rad * rAbs(ir));
         // Psi(r) = int(jl(r*G) * Psi(G) * G^2dG)
         result(ir) = integrate(jl * vec, true);
      }
      result *= sqrt(2./PI);
      result.setBasis(*radBasisPtr);
      result.auxData.is = iSpecies;
      result.auxData.ia = iAtom;
      result.auxData.n = char(n);
      result.auxData.l = char(l);
      result.auxData.m = char(m);

      return result;
   }
   // --- radial to radial interpolation
   SX_CLOCK(Timer::RadR2Rad);

   // Get aux data
   int is = vec.auxData.is;
   int ia = vec.auxData.ia;
   int n  = vec.auxData.n;
   int l  = vec.auxData.l;
   int m  = vec.auxData.m;

   const SxVector<double> &radFunc = radBasisPtr->radFunc(is);
   ssize_t newDim = radFunc.getSize ();
   ssize_t dim = rad.getSize ();
   SxVector<CoeffType> x = rad,
                       y = vec;
   SxCubicSpline spline;

   // Check for boundary
   if (rad(dim-1) < radFunc(newDim-1))  {
      double xLast = rad(dim-1);
      double yLast = vec(dim-1);
      double xPreLast = rad(dim-2);
      double yPreLast = vec(dim-2);
      double slope = (yLast - yPreLast)/(xLast - xPreLast);

      double xLastNew = radFunc(newDim-1);
      double slopeRight = 0.0;

     if (fabs(yLast) < 1e-6 || (slope * yLast) > 0)  {
        // force to zero
        x.resize(dim+1,true);
        y.resize(dim+1,true);
        x(dim) = xLastNew;
        y(dim) = 0.0;
     } else {
        // exponential decay
        double alpha = slope / yLast;
        int nPoints = 100;
        x.resize(dim+nPoints,true);
        y.resize(dim+nPoints,true);
        double delta = (xLastNew - xLast) / nPoints;
        for (int i = 0; i < nPoints; i++)  {
           x(dim + i) = xLast + (i+1) * delta;
           y(dim + i) = yLast * exp(alpha * (i+1) * delta);
        }
        slopeRight = y(dim+nPoints-1) * slope;
     }
     if (l < 2) {
        spline = SxCubicSpline (x,y, SxCubicSpline::NaturalHermite, slopeRight);
     } else  {
        spline = SxCubicSpline (x,y, SxCubicSpline::Hermite, 0.0, slopeRight);
     }
   } else  {
      if (l < 2)
         spline = SxCubicSpline (
               x,
               y,
               SxCubicSpline::NaturalHermite);
      else
         spline = SxCubicSpline (
               x,
               y,
               SxCubicSpline::Hermite);
   }

   TPsi result = spline.getY(radFunc);

   // set aux data
   result.auxData.is = is;
   result.auxData.ia = ia;
   result.auxData.n  = char(n);
   result.auxData.l  = char(l);
   result.auxData.m  = char(m);
   result.setBasis (radBasisPtr);

   return result;
}

SxVector<double> SxRadialBasis::fourier(const SxRadialBasis *radGBasisPtr,
                                        const SxVecRef<double> &vec) const
{
   SX_CLOCK (Timer::RadFourier);
   // Get aux data
   int is = vec.auxData.is;
   int ia = vec.auxData.ia;
   int n  = vec.auxData.n;
   int l  = vec.auxData.l;
   int m  = vec.auxData.m;

   // --- same algorithm for r->g or g->r, but naming is for r->g
   const SxVector<double> &radGFunc = radGBasisPtr->getRadFunc();
   int ng = (int)radGFunc.getSize ();
   SxVector<double> result(ng);


   SxCubicSpline spline;
   if (vec.getSize() == 4 * getNElements ())  {
      spline = SxCubicSpline (rad, vec);
   } else if (vec.getSize() == getNElements ()) {
      if (l < 2)  {
         spline = SxCubicSpline (rad, vec,
               SxCubicSpline::NaturalHermite);
      } else {
         spline = SxCubicSpline (rad, vec,
               SxCubicSpline::Hermite);
      }
   } else  {
      cout << "Incompatible sizes between basis and vector!" << endl;
      SX_EXIT;
   }


   double dRadDense = 0.001 * TWO_PI / radGFunc(ng-1);
   double radMax = rad(rad.getSize()-1);
   int nRadPoints = max(int(radMax/dRadDense) + 1, (int)rad.getSize ());
   nRadPoints = min(nRadPoints, 10 * (int)rad.getSize ());
   SxRadialBasis denseRadR (rad(0), radMax, nRadPoints, realSpace);

   SxVector<double> denseVec = spline.getY(denseRadR.rad);


   for (int ig = 0; ig < ng; ig++)  {
      SxVector<double> jl = jsb (l, radGFunc(ig) * denseRadR.rad);
      // Psi(G) = int(jl(r*G) * Psi(r) * r^2dr)
      result(ig) = denseRadR.integrate(jl * denseVec, false);
   }
   result *= sqrt(2./PI);

   // set aux data
   result.setBasis (radGBasisPtr);
   result.auxData.is = is;
   result.auxData.ia = ia;
   result.auxData.n  = char(n);
   result.auxData.l  = char(l);
   result.auxData.m  = char(m);

   return result;
}

// nsfuerivnkjsnbvowuphyr4 kh0ascn kja vg;
// askdf hkjfhsakdjlf hskadhfjklsadhf adfhjklsadhfkjl


SxVector<PrecCoeffG> SxRadialBasis::toGBasis (const SxGBasis *gBasisPtr,
                                              const SxVecRef<double> &vec) const
{
   SX_CHECK (!realSpace);
   SX_CHECK (vec.getSize() > 0);
   SX_CLOCK (Timer::RadG2Gk);
   SX_CHECK (vec.getBasisPtr () == this);
   SX_CHECK (gBasisPtr->structPtr);
   SX_CHECK (gBasisPtr->structPtr->cell.volume > 0.,
             gBasisPtr->structPtr->cell.volume);

   // Get auxData
   int iSpecies = vec.auxData.is;
   int iAtom    = vec.auxData.ia;
   int n        = vec.auxData.n;
   int l        = vec.auxData.l;
   int m        = vec.auxData.m;


   // Interpolate to the gBasis points
   int lastRad = (int)rad.getSize()-1;
   SxVector<double> gVec = sqrt(gBasisPtr->g2);
   int lastG = (int)gVec.getSize() - 1;
   if (rad(0) > gVec(0)) {
      cout << "Cannot interpolate to outer left region." << endl;
      SX_EXIT;
   }
   if (rad(lastRad) + 1e-10 < gVec(lastG) ) {
      cout << "Cannot interpolate to outer right region." << endl;
      sxprintf("max. g of radial G basis: %.15f\n", rad(lastRad));
      sxprintf("max. g of Gk basis:       %.15f\n", gVec(lastG));
      SX_EXIT;
   }

   SxCubicSpline spline;
   if (vec.getSize() == 4 * getNElements ())  {
      spline = SxCubicSpline (rad, vec);
   } else if (vec.getSize() == getNElements ()) {
      if (l < 2)  {
         spline = SxCubicSpline (rad, vec,
               SxCubicSpline::NaturalHermite);
      } else  {
         spline = SxCubicSpline (rad, vec,
               SxCubicSpline::Hermite);
      }
   } else  {
      cout << "Incompatible sizes between basis and vector!" << endl;
      SX_EXIT;
   }

   SX_START_TIMER(Timer::SplineGetY);
   PsiG result = spline.getY(gVec);
   SX_STOP_TIMER(Timer::SplineGetY);

   // Multiply spherical harmonics if (l,m) is physical
   if (m != NONE_M)   {
      result *= SxYlm::getYlmNormFactor(l,m) * gBasisPtr->getYlm(l,m);
   }

   // Shift to atom iAtom if specified
   if (iAtom != -1)  {
      SX_CLOCK (Timer::PhaseFactors);
      result *= gBasisPtr->getPhaseFactors(iSpecies,iAtom);
   }
   // G+k normalization factor
   result *= sqrt(TWO_PI*TWO_PI*TWO_PI/gBasisPtr->structPtr->cell.volume);

   // set aux data
   result.setBasis (gBasisPtr);
   result.auxData.is = iSpecies;
   result.auxData.ia = iAtom;
   result.auxData.n  = char(n);
   result.auxData.l  = char(l);
   result.auxData.m  = char(m);

   return result;
}

double SxRadialBasis::integrate (const SxVecRef<double> &integrand, bool useSpline) const
{
   SX_CLOCK(Timer::IntegrateRad);
   SX_CHECK(integrand.getSize() > 0);
   SX_VALIDATE_VECTOR(integrand);

   // general Simpson integration sheme 3rd Order
   int dim = (int)rad.getSize ();
   double result = 0.0;
   TPsi f = integrand * rad * rad;
   if (!useSpline)  {
      int i = 0;
      // general Simpson integration sheme 3rd Order
      while (i < dim - 1)  {
         if (i < dim - 3)  {
            result += simpsonTerm(rad(i),rad(i+1),rad(i+2),rad(i+3),f(i),f(i+1),f(i+2),f(i+3));
            i += 3;
         } else if (i == dim - 3)  {
            result += simpsonTerm(rad(i),rad(i+1),rad(i+2),f(i),f(i+1),f(i+2));
            i += 2;
         } else if (i == dim - 2)  {
            result += simpsonTerm(rad(i),rad(i+1),f(i),f(i+1));
            i += 1;
         }
      }
   } else {
      TPsi coeffs = toSpline(f);
      // Integration scheme for cubic Splines S = a x^3 + b x^2 + c x + d
      // Setup Coefficients (integrated analytically S dx) a-> 1/4 a, b -> 1/3 b, c -> 1/2 c, d-> d
      for (int i = 0; i < dim - 1; i++)  {
         double h = rad(i+1)-rad(i);
         result += coeffs(4*i+0) * h
                 + coeffs(4*i+1) * h*h / 2.0
                 + coeffs(4*i+2) * h*h*h / 3.0
                 + coeffs(4*i+3) * h*h*h*h / 4.0;
      }
   }

   return result;
}

double SxRadialBasis::simpsonTerm (double x1, double x2, double x3, double x4,
      double f1, double f2, double f3, double f4) const
{
   double T1,T2,T3,T4, result;

   T1 = f1*(6.*x2*x3 + (2.*(x1-x2-x3)+x4)*x4-x1*(4.*x2+4.*x3-3.*x1)) / ((x2-x1)*(x1-x3));
   T2 = f2*((x1-x4)*(x1-x4)*(x1-2.*x3+x4))/((x2-x1)*(x2-x3)*(x2-x4));
   T3 = f3*((x1-x4)*(x1-x4)*(x1-2.*x2+x4))/((x1-x3)*(x2-x3)*(x3-x4));
   T4 = f4*(6.*x2*x3 + (2.*(x4-x2-x3)+x1)*x1-x4*(4.*x2+4.*x3-3.*x4)) / ((x4-x2)*(x3-x4));;

   result = (x1-x4)/12. * (T1+T2+T3+T4);

   return result;

}

double SxRadialBasis::simpsonTerm (double x1, double x2, double x3, double f1, double f2, double f3) const
{
   double T1,T2,T3, result;

   T1 = f1*(x2-x3)*(2.*x1-3.*x2+x3);
   T2 = f2*(x1-x3)*(x1-x3);
   T3 = f3*(x2-x1)*(x1-3.*x2+2.*x3);

   result = (x3-x1)/(6.*(x1-x2)*(x2-x3)) * (T1+T2+T3);

   return result;
}

double SxRadialBasis::simpsonTerm (double x1, double x2, double f1, double f2) const
{
   double result = 0.5 * (f1+f2) * (x2-x1);

   return result;
}


SxVector<double> SxRadialBasis::jsb (int l, const SxVecRef<double> &z)
{

   SX_CHECK (z.getSize () > 0, z.getSize());

   int zSize = (int)z.getSize ();
   SxVector<double> vec(zSize);
   // --- use generic jsb
#ifdef USE_OPENMP
#   pragma omp parallel for
#endif
   for (int i = 0; i < zSize; ++i)
      vec(i) = SxYlm::jsb(l, z(i));
   SX_VALIDATE_VECTOR (vec);
   return vec;
}

SxVector<double> SxRadialBasis::toSpline (const SxVecRef<double> &vec) const
{
   SX_CLOCK(Timer::VecToSpline);
   // get auxdata
   int is = vec.auxData.is;
   int ia = vec.auxData.ia;
   int n  = vec.auxData.n;
   int l  = vec.auxData.l;
   int m  = vec.auxData.m;

   SxCubicSpline spline;
   if (l < 2)
      spline = SxCubicSpline (rad, vec,
            SxCubicSpline::NaturalHermite);
   else
      spline = SxCubicSpline (rad, vec,
            SxCubicSpline::Hermite);

   SxVector<double> result = spline.getSpline ();

   // set Basis and auxData
   result.setBasis(this);
   result.auxData.is = is;
   result.auxData.ia = ia;
   result.auxData.n  = char(n);
   result.auxData.l  = char(l);
   result.auxData.m  = char(m);

   return result;
}

SxVector<double> SxRadialBasis::toVec (const SxVecRef<double> &vec) const
{
   SX_CLOCK(Timer::SplineToVec);
   // get auxdata
   int is = vec.auxData.is;
   int ia = vec.auxData.ia;
   int n  = vec.auxData.n;
   int l  = vec.auxData.l;
   int m  = vec.auxData.m;

   SxCubicSpline spline (rad, vec);

   SxVector<double> result = spline.getY (rad);

   // set Basis and auxData
   result.setBasis(this);
   result.auxData.is = is;
   result.auxData.ia = ia;
   result.auxData.n  = char(n);
   result.auxData.l  = char(l);
   result.auxData.m  = char(m);

   return result;
}
