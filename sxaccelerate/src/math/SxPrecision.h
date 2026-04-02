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

#ifndef _SX_PRECISION_H_
#define _SX_PRECISION_H_

#include <SxComplex.h>

// --- precision of low level types
typedef int                Int4;
typedef long int           Int8;
typedef float              Real4;
typedef double             Real8;
typedef SxComplex8         Cmplx8;
typedef SxComplex16        Cmplx16; 

// --- precision of several physical variables
//     beware of mixing different precisions. mixing may cause
//     slow down the code due to many type cast commands!
typedef Real8              PrecG;
typedef Int4               PrecFFTIdx;
                        
//#define WAVES_SINGLE_PREC  1
#ifdef  WAVES_SINGLE_PREC
   typedef Cmplx8          PrecCoeffG;
   typedef Cmplx8          PrecPotNl;
   typedef Cmplx16         PrecCoeffR;
#else  /* WAVES_SINGLE_PREC */
   typedef Cmplx16         PrecCoeffG;
   typedef Cmplx16         PrecPotNl;
   typedef Cmplx16         PrecCoeffR;
#endif /* WAVES_SINGLE_PREC */
   
typedef Cmplx16         PrecPhase;
typedef Real8           PrecPhi;
 
typedef Cmplx16            PrecPhaseG;
typedef Real8              PrecPhiG;
typedef Cmplx16            PrecRhoG;
typedef Real8              PrecRhoR;
typedef Cmplx16            PrecEffPotG;   
typedef Real8              PrecEffPotR;
                        
typedef Real8              PrecEps;
typedef Real8              PrecWeights;
typedef Real8              PrecFocc;
typedef Real8              PrecEnergy;    

typedef Real8              PrecTauR;
                        
typedef Real8              PrecForcesR;   
typedef Cmplx16            PrecForcesG;   

#ifdef FFT_SINGLE_PREC
   typedef float         SX_FFT_REAL;
   typedef SxComplex8    SX_FFT_COMPLEX;
#else
   typedef double        SX_FFT_REAL;
   typedef SxComplex16   SX_FFT_COMPLEX;
#endif /* FFT_SINGLE_PREC */


#endif  // _SX_PRECISION_H_

