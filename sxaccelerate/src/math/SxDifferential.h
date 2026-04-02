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

#ifndef _SX_DIFFERENTIAL_H_
#define _SX_DIFFERENTIAL_H_
#include <SxMath.h>
#include <SxVector.h>

/** \brief Numerical differential operator d/d n


    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_MATH SxDifferential
{
   public:
      /// Prefactors for (2N+1)-point formulas
      SxArray<SxVector<double> > preFac;

      /// Constructor
      SxDifferential ()
      {
         // empty
      }

      /// Constructor
      SxDifferential (int n);

      /// Setup routine
      void setup (int n);

      /// Apply at specific point
      double apply (const SxVector<double>& f, ssize_t i) const;

      /// Apply to function
      SxVector<double> apply (const SxVector<double> &f) const;

      /** For xc potentials, we need to integrate the variational
          derivative of the numerical gradient, i.e.,
          \f[
          \frac{\delta E[\dots, \rho'}{\delta \rho(r)}
          = \dots + \int dr' 
          \frac{\delta E[\dots, \rho'}{\delta \rho'(r')}
          \frac{\delta \rho'(r')}{\delta \rho(r)}
          \f]
          This function adds the contribution of one specific point.
          @param fPtr      the resulting potential
          @param i         point r'
          @param variation variational derivative of functional wrt gradient
                           at r'
        */
      void addVariation (double variation, SxVector<double>& f, ssize_t i, ssize_t col = 0) const;
};

inline double SxDifferential::apply (const SxVector<double>& f, ssize_t i) const
{
   SX_CHECK (i >=0 && i < f.getSize (), i, f.getSize ());
   ssize_t n = (preFac.getSize () - 1) >> 1;
   ssize_t o;
   if (i < n)                      // high-side shifted formula
      o = i;
   else if (i + n >= f.getSize ()) // low-side shifted formula
      o = i + preFac.getSize () - f.getSize ();
   else                            // symmetric formula 
      o = n;
   SX_CHECK (o >= 0 && o < preFac.getSize (), o, preFac.getSize ());

   double res = 0.;
   const SxVector<double> &c = preFac(o);
   //const double *fio = &f(i-o);
   for (int j = 0; j < c.getSize (); ++j)
      res += f(i+j-o) * c(j);
      //res += fi[j] * c(j);

   return res;
}

inline SxVector<double> SxDifferential::apply (const SxVector<double> &f) const
{
   SxVector<double> res(f.getSize ());
   for (ssize_t j = 0; j < f.getSize (); ++j)
      res(j) = apply (f, j);
   return res;
}

inline void SxDifferential::addVariation (double variation, SxVector<double> & f, ssize_t i, ssize_t col)
const
{
   ssize_t N = f.getNRows (); // size of f's columns
   SX_CHECK (i >=0 && i + col * N < f.getSize (), i, f.getSize ());
   ssize_t n = (preFac.getSize () - 1) >> 1;
   ssize_t o;
   if (i < n)           // high-side shifted formula
      o = i;
   else if (i + n >= N) // low-side shifted formula
      o = i + preFac.getSize () - N;
   else                 // symmetric formula 
      o = n;
   SX_CHECK (o >= 0 && o < preFac.getSize (), o, preFac.getSize ());

   // Go to right location
   i += col * N - o;

   const SxVector<double> &c = preFac(o);
   for (int j = 0; j < c.getSize (); ++j)
      f(i+j) += c(j) * variation;
}

#endif /* _SX_DIFFERENTIAL_H_ */
