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

#ifndef _SX_PDOS_H_
#define _SX_PDOS_H_

#include <SxExt.h>
#include <SxPartialWaveBasis.h>
#include <SxPAWSet.h>
#include <SxFermi.h>
class SxHamSolver;

/** \brief Compute partial DOS for PAW

    \author Siyuan Zhang, C. Freysoldt */
class SX_EXPORT_EXT SxPDOS
{
   public:
      SxConstPtr<SxPartialWaveBasis> pBasisPtr;
      SxConstPtr<SxPAWPot> pawPotPtr;

      double broadening;
      double eMin, eMax;
      SxString outFile;
      SxList<int> atoms;
      int nPerPeak;
      bool verbose;

      SxPDOS (const SxConstPtr<SxPartialWaveBasis> &pBasisIn,
              const SxConstPtr<SxPAWPot> &pawPotIn);

      void compute (const SxPAWSet &waves, const Eps &eps);
      void getRange (const Eps &eps, double shift);

};

/// Interface for input file
SX_EXPORT_EXT
void computePDOS (SxHamSolver *pot, const SxSymbolTable *table, bool calc);

#endif /* _SX_PDOS_H_ */
