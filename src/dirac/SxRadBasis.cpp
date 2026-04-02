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

#include <SxRadBasis.h>
#include <SxRadialBasis.h>
#include <SxGBasis.h>
#include <SxConstants.h>
#include <SxMathLib.h>
#include <stdio.h>
#include <math.h>
#include <SxYlm.h>
#include <SxTextIO.h>
#include <SxCubicSpline.h>


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------


SxRadBasis::SxRadBasis ()
{
   registerMemoryObservers ();
}


SxRadBasis::SxRadBasis (const SxArray<SxVector<double> > &radFunc_,
                        const SxArray<Real8>               &logDr_)
{
   set(radFunc_, logDr_);
}

void SxRadBasis::set (const SxArray<SxVector<double> > &radFunc_,
                      const SxArray<Real8>               &logDr_)
{
   this->radFunc = radFunc_;
   this->logDr   = logDr_;

   registerMemoryObservers ();
}

SxRadBasis::SxRadBasis (double rMin, double rMax, int nPoints)
{
   set(rMin, rMax, nPoints);
}

SxRadBasis::SxRadBasis (const SxVector<double> &rMin, 
                        const SxVector<double> &rMax,
                        const SxVector<int> &nPoints)
{
   set(rMin, rMax, nPoints);
}

void SxRadBasis::set (double rMin, double rMax, int nPoints)
{
   SX_CHECK (nPoints > 1, nPoints);
   SX_CHECK (rMin > 0., rMin);
   SX_CHECK (rMax > rMin, rMax, rMin);

   radFunc.resize (1); // one species
   logDr.resize (1);
   SxVector<double> r(nPoints);
   double ldr = log(rMax/rMin) / double(nPoints - 1);
   for (int i = 0; i < nPoints; i++)
      r(i) = rMin * exp(ldr * i);
   SX_CHECK (fabs(r(nPoints-1)-rMax) < 1e-12 * rMax, r(nPoints-1), rMax);
   radFunc(0) = r;
   this->logDr(0) = ldr;
   registerMemoryObservers ();
}

void SxRadBasis::set (const SxVector<double> &rMin, 
                      const SxVector<double> &rMax,
                      const SxVector<int> &nPoints)
{
   SX_CHECK(rMin.getSize () == rMax.getSize ());
   SX_CHECK(rMin.getSize () == nPoints.getSize ());

   ssize_t nSpecies = rMin.getSize ();
   for (ssize_t iSpecies = 0; iSpecies < nSpecies; iSpecies++)  {
      SX_CHECK (nPoints(iSpecies) > 1, nPoints(iSpecies));
      SX_CHECK (rMin(iSpecies) > 0., rMin(iSpecies));
      SX_CHECK (rMax(iSpecies) > rMin(iSpecies),
                rMax(iSpecies), rMin(iSpecies));

      radFunc.resize (nSpecies); 
      logDr.resize (nSpecies);
      int dim = nPoints(iSpecies);
      SxVector<double> r(dim);
      double ldr = log(rMax(iSpecies)/rMin(iSpecies)) / double(dim - 1);
      for (int i = 0; i < dim; i++)
         r(i) = rMin(iSpecies) * exp(ldr * i);
      SX_CHECK (fabs(r(dim-1)-rMax(iSpecies)) < 1e-12 * rMax(iSpecies),
                r(dim-1), rMax(iSpecies));
      radFunc(iSpecies) = r;
      this->logDr(iSpecies) = ldr;
   }
   registerMemoryObservers ();
}

SxRadBasis::SxRadBasis (const SxString &file)
{
   this->read(file);
}

SxRadBasis::SxRadBasis (const SxBinIO &io)
{
   this->read(io);
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
SxRadBasis::~SxRadBasis ()
{
   deregisterAll ();
}

void SxRadBasis::read(const SxString &file)
{
try  {
      SxBinIO io (file, SxBinIO::BINARY_READ_ONLY);
      read (io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

void SxRadBasis::read (const SxBinIO &io)
{
   SxArray<SxVector<double> > radFunc_;
   SxArray<double> logDr_;

   try {
      //get dimensions
      int nSpecies = io.getDimension ("nSpecies");
      radFunc_.resize(nSpecies);
      logDr_.resize(nSpecies);
      for(int iSpecies = 0; iSpecies < nSpecies; iSpecies++)   {
         SxString dimRadName = "dimRad-" + SxString(iSpecies);
         int dimRad = io.getDimension (dimRadName);
         radFunc_(iSpecies).resize(dimRad);
         SxVector<double> &Vec = radFunc_(iSpecies);
         SxString radialName = "radFunc-" + SxString(iSpecies); 
         io.read (radialName,&Vec,dimRad);
         SxString logDrName = "logDr-"+SxString(iSpecies);
         double value;
         io.read(logDrName, &value);
         logDr_(iSpecies) = value;
      }

   } catch (SxException e)   {
      e.print ();
      SX_EXIT;
   }


   set(radFunc_, logDr_);
}

SxRadBasis SxRadBasis::readMesh(const SxString &file)
{
   SxRadBasis result;
   try  {
      SxBinIO io (file, SxBinIO::BINARY_READ_ONLY);
      result = readMesh (io);
      io.close ();
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   return result;
}

SxRadBasis SxRadBasis::readMesh (const SxBinIO &io)
{
   SxArray<SxVector<double> > radFunc_;
   SxArray<double> logDr_;

   try {
      //get dimensions
      int nSpecies = io.getDimension ("nSpecies");
      radFunc_.resize(nSpecies);
      logDr_.resize(nSpecies);
      for(int iSpecies = 0; iSpecies < nSpecies; iSpecies++)   {
         SxString dimRadName = "dimRad-" + SxString(iSpecies);
         int dimRad = io.getDimension (dimRadName);
         radFunc_(iSpecies).resize(dimRad);
         SxVector<double> &Vec = radFunc_(iSpecies);
         SxString radialName = "radFunc-" + SxString(iSpecies); 
         io.read (radialName,&Vec,dimRad);
         SxString logDrName = "logDr-"+SxString(iSpecies);
         double value;
         io.read(logDrName, &value);
         logDr_(iSpecies) = value;
      }

   } catch (SxException e)   {
      e.print ();
      SX_EXIT;
   }


   SxRadBasis result(radFunc_, logDr_);
   return result;
}

int SxRadBasis::addMesh (const SxVecRef<double> &radFuncIn,
                               double             logDrIn)
{
   int is = int(radFunc.getSize ());
   int nSpecies = is + 1;
   radFunc.resize (nSpecies, true);
   logDr.resize (nSpecies, true);

   radFunc(is) = radFuncIn;
   logDr(is) = logDrIn;

   cashedRl.resize (0);
   gBases.resize (0);

   return is;
}

double SxRadBasis::getRMax () const
{
   SxVector<double> rMax (getNSpecies ());

   for (int is = 0; is < getNSpecies (); is++) rMax(is) = getRMax(is);

   return rMax.maxval ();
}

SxVector<SxRadBasis::CoeffType>
SxRadBasis::changeRadBasis (const SxRadBasis *basisPtr,
                            const SxVecRef<CoeffType> &vec) const
{
   SX_CHECK(vec.getBasisPtr() == this);
   
   const SxRadBasis &basis = *basisPtr;
   
   int is = vec.auxData.is;
   int l = vec.auxData.l;

   // Check for identical basis
   if (basis.radFunc(is).getSize() == radFunc(is).getSize())
      if ((radFunc(is) - basis.radFunc(is)).norm() < 1e-10) return (1.0 * vec);

   SX_CHECK(is < basis.radFunc.getSize(), is, basis.radFunc.getSize());
   // Fourierinterpolation
   // result = (newRad|radG|vec);
   /*
   SxRadialBasis radGBasis(0.0, 30.0, 3000, false, SxRadialBasis::Linear);
   SxVector<CoeffType> result = *basisPtr | ( radGBasis | vec);
   */
   // SplineInterpolation
   ssize_t dim = radFunc(is).getSize ();
   ssize_t dimNew = basis.radFunc(is).getSize();
   SxVector<CoeffType> x,y;
   SxCubicSpline spline;
   // Check for boundary
   if ((basis.radFunc(is)(dimNew-1) - radFunc(is)(dim-1)) > 1e-4)  {
      x = radFunc(is);
      y = vec;

      double xLast = radFunc(is)(dim-1);
      double yLast = vec(dim-1);
      double xPreLast = radFunc(is)(dim-2);
      double yPreLast = vec(dim-2);
      double slope = (yLast - yPreLast)/(xLast - xPreLast);

      double xLastNew = basis.radFunc(is)(dimNew-1);
      double slopeRight = 0.0;

     if (fabs(yLast) < 1e-6 || (slope * yLast) > 0)  { 
        // force to zero
        x.resize(dim+1,true);
        y.resize(dim+1,true);
        if (xLastNew - xLast > 5.0)
           x(dim) = xLastNew;
        else
           x(dim) = xLastNew + 5.0;
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
        spline 
           = SxCubicSpline (
         x,y, 
         SxCubicSpline::NaturalHermite,
         slopeRight);
     } else  {
        spline 
           = SxCubicSpline (
         x,y, 
         SxCubicSpline::Hermite,
         0.0,
         slopeRight);
     }
   } else  {
      x = radFunc(is);
      y = vec;
      if (l < 2) {
         spline = SxCubicSpline(
               x,y, 
               SxCubicSpline::NaturalHermite);
      } else  {
         spline  = SxCubicSpline (
               x,y, 
               SxCubicSpline::Hermite);
      }
   }
   
   SxVector<CoeffType> result = spline.getY(basis.radFunc(is));

   cout << SX_SEPARATOR;
   cout << "WARNING: Change of RADBASIS by interpolation!" << endl;
   SxString file = "ORIG-" + SxString(is) + SxString(l) + ".dat";
   cout << "Check " << file << endl;
   SxTextIO (file).writeXYPlot(radFunc(is), vec);

   file = "INTERPOL-" + SxString(is) + SxString(l) + ".dat";
   cout << "Check " << file << endl;
   SxTextIO(file).writeXYPlot(basis.radFunc(is), result);
   cout << SX_SEPARATOR;
   
   result.setBasis(basisPtr);
   result.auxData.is = vec.auxData.is;
   result.auxData.ia = vec.auxData.ia;
   result.auxData.n  = vec.auxData.n;
   result.auxData.l  = vec.auxData.l;
   result.auxData.m  = vec.auxData.m;

   return result;
}


SxVector<SxComplex16>
SxRadBasis::toPWBasis (const SxGBasis *basis,
                       const SxVecRef<CoeffType> &vec) const
{
   SX_CLOCK (Timer::RadTotal);
   SX_CHECK (basis->structPtr);
   SX_CHECK (basis->structPtr->cell.volume > 0.);
   double rNorm = FOUR_PI / sqrt (basis->structPtr->cell.volume);
   SxVector<SxComplex16> res;

   int is = vec.auxData.is;  // TODO: ugly
   int ia = vec.auxData.ia;  // TODO: ugly
   int l  = vec.auxData.l;   // TODO: ugly
   int m  = vec.auxData.m;   // TODO: ugly

   bool useCache =  cashedVec.getSize () == vec.getSize () 
                 && cashedVec.auxData.l  != l
                 && cashedVec.auxData.is != is;
   if (useCache)  {
      // --- check cached vector elements
      for (int i = 0; i < vec.getSize (); ++i)  {
         if (vec(i) != cashedVec(i))  {
            useCache = false;
            break;
         }
      }
   }

   if (!useCache) {
      // --- clear cache
      for (int jk = 0; jk < cashedRl.getSize (); ++jk)
         cashedRl(jk) = SxVector<double> ();
      cashedVec = vec; // copy
      cashedVec.auxData.l = char(l);
      cashedVec.auxData.is = is;
   }

   int ik = int(gBases.findPos (basis));
   if (ik == -1)  {
      // new G basis
      gBases.append (basis);
      int nk = int(gBases.getSize ());
      ik = nk - 1;
      cashedRl.resize (nk, true);

      UPDATE_MEMORY (cashedRl);
      UPDATE_MEMORY (cashedVec);
   }

   if (cashedRl(ik).getSize () > 0) {
      // use cashed result
      // cout << "Use cash: ik =" << ik << "; l=" << l << "; m=" << m << endl;
   } else {
      // cout << "Set cash: ik =" << ik << "; l=" << l << "; m=" << m << endl;
      cashedRl(ik) = toPWBasis (radFunc(is), vec, *basis, l, logDr(is));
   }

   // multiply with Ylm
   res = cashedRl(ik) * basis->getYlm (l,m);
   // normalization factors
   res *= SxYlm::getYlmNormFactor(l,m) * rNorm;

   // --- do we calculate <G|r><r|Psi> or <G|r><r|T*Psi>?
   if (ia >= 0)  {
      SX_CLOCK (Timer::Phase);
      res *= basis->getPhaseFactors(is,ia);
   }

   res.setBasis (basis);

   res.auxData.is = is; // TODO: ugly
   res.auxData.ia = ia; // TODO: ugly
   res.auxData.l  = char(l);  // TODO: ugly
   res.auxData.m  = char(m);  // TODO: ugly

   return res;

}

SxVector<double> SxRadBasis::toRadialBasis (const SxRadialBasis *radBasisPtr,
                                          const SxVecRef<CoeffType> &vec) const
{
   SX_CHECK (radBasisPtr);
   if (radBasisPtr->realSpace)  {
      SX_CLOCK(Timer::rad2radR);

      // Get aux data
      int is = vec.auxData.is;
      int ia = vec.auxData.ia;
      int n  = vec.auxData.n;
      int l  = vec.auxData.l;
      int m  = vec.auxData.m;

      const SxVector<double> &radRFunc = radBasisPtr->getRadFunc();

      ssize_t dim = radFunc(is).getSize ();
      ssize_t newDim = radRFunc.getSize ();
      SxVector<CoeffType> x = radFunc(is),
                          y = vec;
      SxCubicSpline spline;

      // Check for boundary
      if ((radRFunc(newDim-1) - radFunc(is)(dim-1)) > 1e-4)  {
         double xLast = radFunc(is)(dim-1);
         double yLast = vec(dim-1);
         double xPreLast = radFunc(is)(dim-2);
         double yPreLast = vec(dim-2);
         double slope = (yLast - yPreLast)/(xLast - xPreLast);

         double xLastNew = radRFunc(newDim-1);
         double slopeRight = 0.0;

        if (fabs(yLast) < 1e-6 || (slope * yLast) > 0)  {
           // force to zero
           x.resize(dim+1,true);
           y.resize(dim+1,true);
           if (xLastNew - xLast > 5.0)
              x(dim) = xLastNew;
           else
              x(dim) = xLastNew + 5.0;
           y(dim) = 0.0;
        } else {
           // exponential decay
           double alpha = log(yLast/yPreLast) / (xLast - xPreLast);
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
            spline = SxCubicSpline (x, y, SxCubicSpline::NaturalHermite);
         else
            spline = SxCubicSpline (x, y, SxCubicSpline::Hermite);
      }

      TPsi result = spline.getY(radRFunc);

      // set aux data
      result.auxData.is = is;
      result.auxData.ia = ia;
      result.auxData.n  = char(n);
      result.auxData.l  = char(l);
      result.auxData.m  = char(m);
      result.setBasis (radBasisPtr);

      return result;
   }
   // else rad -> G
   SX_CLOCK (Timer::rad2radG);

   // Get aux data
   int is = vec.auxData.is;
   int ia = vec.auxData.ia;
   int n  = vec.auxData.n;
   int l  = vec.auxData.l;
   int m  = vec.auxData.m;

   const SxVector<double> &radGFunc = radBasisPtr->getRadFunc();
   int ng = (int)radGFunc.getSize ();
   SxVector<double> result(ng);
   SxVector<CoeffType> vec_rad3 = vec * radFunc(is).cub();

   for (int ig = 0; ig < ng; ig++)  {
      SxVector<double> jl = SxRadialBasis::jsb (l, radGFunc(ig) * radFunc(is));
      result(ig) = (vec_rad3 * jl).integrate(logDr(is));
   }
   result *= sqrt(2./PI);

   // set aux data
   result.auxData.is = is;
   result.auxData.ia = ia;
   result.auxData.n  = char(n);
   result.auxData.l  = char(l);
   result.auxData.m  = char(m);
   result.setBasis (radBasisPtr);

   return result;
}

inline double weightSimpson (int i, int n)
{
   SX_CHECK (i >= 0 && i < n, i, n);
   if (i == 0) return 1. / 3.;
   int d = n - i;
   if (d > 4) return ((i & 1) ? 4. : 2. ) / 3.;
   if (n & 1)  {
      if (d == 1) return 1./3.;
      if (d & 1) return 2./3.;
      return 4. / 3.;
   } // else
   if (d == 1) return 3./8.;
   if (d == 4) return 17./24.;
   return 9. / 8.;
}


//------------------------------------------------------------------------------
// This function projects the radial vector onto plane waves.
//------------------------------------------------------------------------------
SxVector<double> SxRadBasis::toPWBasis (const SxVecRef<double> &rad,
                                        const SxVecRef<double> &psi,
                                        const SxGBasis &pwBasis,
                                        int l, Real8 logDrIn) const
{
   SX_CHECK (rad.getSize() == psi.getSize() && rad.getSize() > 0,
             rad.getSize(), psi.getSize());
   SX_CLOCK (Timer::JsbInt);

   ssize_t ig, ng = pwBasis.ng;
   int nr = (int)psi.getSize ();
   SxVector<double> psi_rad3, jl, g(pwBasis);
   Real8                        gLast = 0., gAbsLast;

   // ---     radFunc * radX
   psi_rad3 = psi     * rad.cub();

   gAbsLast = -1.;
#ifdef USE_OPENMP
#pragma omp parallel
#pragma omp for firstprivate(gLast,gAbsLast)
#endif
   for (ig = 0; ig < ng; ig++)  {
      Real8 gAbs  = sqrt (pwBasis.g2(ig));
      if ( ig == 0 || fabs(gAbsLast - gAbs) > 1e-10 )  {
         gLast = 0.;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:gLast)
#endif
         for (int ir = 0; ir < nr; ++ir)
            gLast += weightSimpson (ir, nr) * psi_rad3(ir)
                     * SxYlm::jsb(l, gAbs * rad(ir));
         gLast *= logDrIn;
         g(ig) = gLast;

         gAbsLast = gAbs;
      }  else
         g(ig) = gLast;  // integrals for equal |G|s are equal!
   }

   g.setBasis (&pwBasis);
   return g;
}

Real8 SxRadBasis::tr (const SxVecRef<double> &vec) const
{
   int is = vec.auxData.is;
   SX_CHECK (is >= 0 && is < radFunc.getSize (), is, radFunc.getSize ());
   // TODO: check if this is time critical
   // explicit for loop might be faster
   return (vec * radFunc(is).cub ()).integrate (logDr(is));
}



//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
// return LegendreP [l, x]; x = cos(theta)
Real8 SxRadBasis::pl (int l, Real8 x)
{
   SX_CHECK (l >= SxQuantumNumbers::s && l <= SxQuantumNumbers::f, l);

   switch ( l )  {
      case SxQuantumNumbers::s :
         return 1.;
         break;
      case SxQuantumNumbers::p :
         return x;
         break;
      case SxQuantumNumbers::d :
         return 0.5 * ( 3.*x*x - 1. );
         break;
      case SxQuantumNumbers::f :
         return 0.5 * ( 5.*x*x*x - 3.*x );
         break;
      default:
         SX_EXIT;
   }
   return 0.;
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Real8 SxRadBasis::cosTheta (int ig, int jg, const SxGBasis &g)
{
   Real8 giTgj = sqrt ( g.g2(ig) * g.g2(jg)  );
   if (giTgj <  1.e-7) return 1.;
   else return ( ( g.getG(ig) ^ g.getG(jg) ) / giTgj ) ;
}


// see numrec eq. 4.1.14
Real8 SxRadBasis::integrate (const SxVecRef<double> &vec) const
{
   SX_EXIT;
   SX_CHECK (vec.getSize() >= 10, vec.getSize());

   double res = 0.;

   const double C0=-31./48.;
   const double C1=+11./48.;
   const double C2= -5./48.;
   const double C3= +1./48.;


   double r1  = radFunc(0)(0);
   double dex = logDr(0);

   const double xexp = exp(dex);
   int nr = (int)vec.getSize();


   // --- integrate from 0 to first grid point
   double r2=r1*xexp;
   double r3=r2*xexp;
   double s21=(vec(2)-vec(1))/(r2-r1);
   double s31=(vec(3)-vec(1))/(r3-r1);
   double g1=(s21-s31)/(r2-r3);
   double b=s21-g1*(r1+r2);
   double a=vec(1)-b*r1-g1*r1*r1;
   res = a*r1 + .5*b*r1*r1 + 1./3.*g1*r1*r1*r1;


   // --- map onto linear grid
   SxVector<double> f (nr);
   double ri=r1/xexp;
   int ir;
   for (ir=0; ir < nr; ir++)  {
      ri=ri*xexp;
      f(ir)=vec(ir)*dex*ri;
   }

   // --- summation
   res=res+C0*f(1)+C1*f(2)+C2*f(3)+C3*f(4);
   for (ir=0; ir < nr; ir++)  {
      res=res+f(ir);
   }
   if (nr <= 4)
      res=res+C0*f(nr-1)+C1*f(nr-2)+C2*f(nr-3)+C3*f(nr-4);
   else
      res=res-0.5*f(nr-1);

   return res;


}

SxVector<double> SxRadBasis::realYlm (int lmax, 
                                        const SxGBasis &G,
                                        const Coord &dG)
{
   SX_CHECK (lmax >= 0, lmax);
   int nl = lmax + 1;
   int nLm = nl*nl;
   int ng = G.ng;
   SX_CHECK (ng >= 0, ng);
   SxVector<double> Ylm(nLm);

   SxVector<double> result (ng, nLm);

   // set up Ylm for all G vectors
   for (int ig = 0; ig < ng; ig++)  {
      Coord Gk = G.getG (ig) + dG;
      if (Gk.normSqr () < 1e-12)  {
         result (ig, 0) = 1.;
         for (int lm = 1; lm < nLm; lm++)
            result(ig,lm) = 0.;
      } else {
         SxYlm::getYlmArray(lmax, Gk(0), Gk(1), Gk(2), &Ylm);
         for (int lm = 0; lm < nLm; lm++)
            result(ig,lm) = Ylm(lm);
      }
   }
   result.setBasis (&G);

   return result;
}


void SxRadBasis::registerMemoryObservers ()
{
   TRACK_MEMORY (radFunc);
   TRACK_MEMORY (cashedRl);
   TRACK_MEMORY (cashedVec);
   TRACK_MEMORY (cashYlm);
   TRACK_MEMORY (YlmNormFactors);
}

// --- Stand alone interpolation function
#include <SxNaturalCubicSpline.h>
SxVector<double> interpolateRad(const SxVecRef<double> &psi,
                                  double r0, double logDr,
                                  const SxVecRef<double> &newRad)
{
   SX_CHECK (logDr > 0., logDr);
   SX_CHECK (r0 > 0., r0);
   // natural cubic spline interpolation
   SxNaturalCubicSpline cubSpline (psi, true);
   // 1st derivative at r0
   double dPsi = (psi(1) - psi(0)) / ((exp(logDr) - 1.) * r0);

   int n = (int)newRad.getSize ();
   SxVector<double> res(n);
   for (int i = 0; i < newRad.getSize (); ++i)  {
      double r = newRad(i);
      if (r <= r0)  {
         // linear extrapolation below r0
         res(i) = dPsi * (r - r0) + psi(0);
      } else {
         double x = log(r/r0) / logDr;
         SX_CHECK ((x - (double)psi.getSize ()) < 1e-10, x,psi.getSize ());
         res(i) = cubSpline.getVal (x);
      }
   }
   SX_VALIDATE_VECTOR(res);
   res.auxData = psi.auxData;
   return res;
}

