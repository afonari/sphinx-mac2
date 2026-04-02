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

#include <SxError.h>
#include <SxNaturalCubicSpline.h>

//------------------------------------------------------------------------------
// SxNaturalCubicSpline-class
//------------------------------------------------------------------------------
// Constructors
// -- Constructor
SxNaturalCubicSpline::SxNaturalCubicSpline ()
{
   // empty
}

SxNaturalCubicSpline::SxNaturalCubicSpline (const SxVector<double> &x, const SxVector<double> &y)
   : xVals (x)
{
   compute (x,y);
}

// -- Destructor
SxNaturalCubicSpline::~SxNaturalCubicSpline ()
{
   // empty
}

SxNaturalCubicSpline::SxNaturalCubicSpline (const SxVector<double> &y,
                                            bool estimateD2)
{
   compute(y, estimateD2);
}

void SxNaturalCubicSpline::compute (const SxVector<double> &y, bool estimateD2)
{
   int n = (int)y.getSize () - 1;
   SX_CHECK (n > 2, n);

   // prepare internal state
   xVals.resize (0);
   polyCoeff.resize (2);

   // set up tridiagonal matrix for coefficients
   SxVector<double> dl(n),d(n+1),du(n);
   d(0) = 1.;
   du(0) = 0.;
   // h_i := x_{i+1} - x_i
   for (int i = 1; i < n; ++i)   {
      dl(i-1) = 1.; /* h_{i-1} */
      d (i)   = 4.; /* 2(h_{i-1} + h_i) */
      du(i)   = 1.; /* h_i */
   }
   dl(n-1) = 0.;
   d(n) = 1.;

   // --- set up right-hand side
   SxVector<double> rhs(n+1);
   if (estimateD2)
      rhs(0) = 2. * y(0) - 5. * y(1) + 4. * y(2) - y(3);
      //rhs(0) = 3. * y(0) - 9. * y(1) + 10. * y(2) - 5. * y(3) + y(4);
   else
      rhs(0) = 0.;
   for (int i = 1; i < n; ++i)  {
      // rhs(i) = 6 * ((y_{i+1}-y_i)/h_i - (y_i - y_{i-1})/h_{i-1})
      rhs(i) = 6. * (y(i+1) + y(i-1) - 2. * y(i));
   }
   if (estimateD2)
      rhs(n) = 2. * y(n) - 5. * y(n-1) + 4. * y(n-2) - y(n-3);
   else
      rhs(n) = 0.;
   // solve tdm * x = rhs
   polyCoeff(1) = symmetricTridiagonalGauss (du, d, dl, rhs);
   polyCoeff(1) /= 6.; // prefactor

   // --- copy original y-data
   SxVector<double> &a0 = polyCoeff(0);
   a0.resize (y.getSize ());
   for (int i = 0; i < a0.getSize (); ++i)
      a0(i) = y(i);

}

SxVector<double> SxNaturalCubicSpline::getVals (const SxVecRef<double> &x) const
{
   SxVector<double> res;
   ssize_t n = x.getSize ();
   res.resize (n);
   SxVecConstIt<double,Compact> xIt = x.begin ();
   SxVecIt<double,Compact> resIt = res.begin ();
   if (xVals.getSize () == 0)  {
      for ( ; resIt != res.end (); ++xIt, ++resIt)
         *resIt = getValY (*xIt);
   } else {
      for ( ; resIt != res.end (); ++xIt, ++resIt)
         *resIt = getValXY (*xIt);
   }
   return res;
}


void SxNaturalCubicSpline::print (const SxVector<double> &diag1, const SxVector<double> &diag2, const SxVector<double> &diag3, const SxVector<double> &b)
{
  ssize_t n = diag2.getSize ();
  
   // TESTOUTPUT
  for (int i = 0; i < n; ++i) {
     for (int iCol = 0; iCol < i - 1; ++iCol)  {
        cout << setw (9) << 0.;
     }
     // lower
     if (i-1 >= 0)  {
        cout << setw (9)<< diag3(i-1);
     }
     // mid
     cout << setw (9) << diag2(i);
     // upper
     if (i < n - 1) {
        cout << setw (9) << diag1(i);
     }
     for (int iCol = i + 2; iCol < n; ++iCol)  {
        cout << setw (9) << 0.;
     }
     cout << setw (9) << b(i) << endl;		
   }
}

//extern "C" {
//#include <f2c.h>
//#include <clapack.h>
//}
SxVector<double> SxNaturalCubicSpline::symmetricTridiagonalGauss (
      const SxVector<double> & upper,
      const SxVector<double> &mid,
      const SxVector<double> &lower,
      const SxVector<double> &rhs)
{
   /* --- LAPACK version
   integer N = mid.getSize (), nrhs=1, info = 0;
   //dgtsv_(&N, &nrhs, lower.elements, mid.elements, upper.elements, 
   //       rhs.elements, &N, &info);

   char fact='N', trans='N';
   SxVector<double> dlf(N-1), df(N), duf(N-1), du2(N-2), res(N), work(3*N);
   integer *ipiv = new integer[N],
           *iwork = new integer[N];
   double rcond, ferr, berr;
   dgtsvx_ (&fact, &trans, &N, &nrhs,
            lower.elements, mid.elements, upper.elements,
            dlf.elements, df.elements, duf.elements, du2.elements, ipiv,
            rhs.elements, &N, res.elements, &N,
            &rcond, &ferr, &berr, work.elements, iwork, &info);
   delete iwork;
   delete ipiv; 

   if (info != 0) {
      cout << "dgtsv failed: info=" << info << endl;
      SX_EXIT;
   }
   return res;
   // return rhs;
   */
   SxVector<double> diag1 (upper), diag2 (mid), diag3 (lower), b (rhs);
   ssize_t n = diag2.getSize ();
  
  // Overview:
  // ---------  
  // diag2(0) diag1(0) 0        0
  // diag3(0) diag2(1) diag1(0) 0
  // 0        ...
  // ...
  // 0        0        0        diag3(n-2) diag2(n-1) diag1(n-1)
  // 0        0        0        0          diag3(n-1) diag2 (n)
  
//  cout << "Before Gaussian Elimination: " << endl; 
//  print (diag1, diag2, diag3, b);
  	
  // Gaussian Elimination
  for (ssize_t i = 0; i < n - 1; ++i)  {
     diag2(i + 1) -= diag1 (i) * diag3 (i) / diag2(i);
     b(i + 1) -= b(i) * diag3(i) / diag2(i);
     diag3(i) = 0;
  }

//  cout << "Gaussian Elimination: " << endl;
//  print (diag1, diag2, diag3, b);
  

  // Backward substitution
  b(n-1) /= diag2(n-1);
  diag2(n-1) = 1;
  for (ssize_t i = n - 2; i >= 0; --i)  {
     b(i) -= diag1(i) * b(i + 1);
     diag1(i) = 0;
     b(i) /= diag2(i);
     diag2(i) = 1;
  } 
//  cout << "Backward substitution: " << endl;
//  print (diag1, diag2, diag3, b);

  return b;
}

double SxNaturalCubicSpline::getValXY (double x) const
{
   // Introducing epsilon which is supposed to guarantee that even if
   // floating point errors occured one would receive the corresponding
   // y-value for the last x-value
   SX_CHECK ((xVals(0) <= x && xVals(xVals.getSize () - 1) + 1e-8>= x), xVals(0), x, xVals(xVals.getSize () - 1));
   
   for (int i = 0; i< xVals.getSize (); ++i){
      if (i == xVals.getSize () - 1 || (xVals(i)<= x && x < xVals(i+1)))  {
         SxVector<double> power (4);
         double dx = x - xVals(i);
         power(0) = 1;
         power(1) = dx;
         power(2) = dx * dx;
         power(3) = dx * dx * dx;
         //cout << "Hello!!!" << endl;         
         return dot (polyCoeff(i), power);
      }     
   }

   SX_EXIT;
   
}


// Methods
void SxNaturalCubicSpline::compute (const SxVector<double> &x,const SxVector<double> &y)
{
   // Testing if the sizes are the same
   SX_CHECK (x.getSize () == y.getSize ());
   
   ssize_t n (x.getSize () - 2);
  
   SxVector<double> h (n + 1), gamma (n);
   for (int i = 0; i < n + 1; ++i)  {
      h(i) = x(i + 1) - x(i);
      if (i > 0 && i < n + 1) {
         gamma(i-1) = 6. * ((y(i + 1) - y(i)) / h(i) - (y(i) - y(i - 1)) / h(i-1));
      }
   }
   
   SxVector<double> diag1 (n-1), diag2 (n), diag3 (n-1), b (n);
   
   // Init vectors
  for (int i = 0; i < n - 1; ++i)  {
     diag1(i) = h(i + 1);
     diag2(i) = 2. * (h(i) + h (i + 1));
     diag3(i) = h(i + 1);
     b(i) = gamma(i);
  }
  b(n - 1) = gamma(n-1);
  diag2(n - 1) = 2. * (h(n - 1) + h(n));

  SxVector<double> beta (n+2);
  beta(0) = 0.;
  beta(n + 1) = 0.;
  b = symmetricTridiagonalGauss (diag1, diag2, diag3, b);
  for (int i = 0; i < b.getSize (); ++i)  {
     beta(i+1) = b(i);
  }

/*  cout << "Beta: " << endl;
  for (int i = 0; i < beta.getSize (); ++i){
     cout << beta(i) << " | ";
  }
  cout << endl;
*/
  SxVector<double> alpha (beta.getSize () - 1);

//  cout << "Alpha:" << endl;
  for (int i = 0; i < alpha.getSize (); ++i)  {
     alpha (i) = (y(i + 1) - y(i)) / h(i) - 1. / 3. * beta(i) * h(i) - 1. / 6. * beta(i+1) * h(i);
//     cout << alpha(i) << " | ";
  }
//  cout << endl;

  polyCoeff.resize (y.getSize ()); 
  //polyCoeff.resize (alpha.getSize ());
  //cout << "polyCoeff"<< polyCoeff.getSize () << endl;
  //cout << "y" << y.getSize () << endl;
  for (int i = 0; i < polyCoeff.getSize (); ++i)  {
     if (i < alpha.getSize ())  {
        polyCoeff(i).resize (4);
        polyCoeff(i)(0) = y(i);
        polyCoeff(i)(1) = alpha(i);
        polyCoeff(i)(2) = beta(i)/2.;
        polyCoeff(i)(3) = (beta(i + 1) - beta(i)) / (6. * h(i));   
     } else  {
        // Necessary to maintain all y-values
        polyCoeff(i).resize (1);
        polyCoeff(i)(0) = y(i);
     }  
  }
}
