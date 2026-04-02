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

#include <SxCliStandard.h>
#include <SxParser.h>
#include <SxStickyFilter.h>
#include <SxTextIO.h>
#include <SxPseudoPot.h>
#include <SxPAWPot.h>
#include <SxAtomicOrbitals.h>
#include <SxPartialWaveBasis.h>
#include <SxPAWOverlap.h>

SxAtomicStructure structureFromInput (const SxString &inFile, bool sxbFile,
                                      bool withMovable)
{
   SxAtomicStructure structure;
   if (sxbFile)  {
      try {
         SxBinIO io (inFile, SxBinIO::BINARY_READ_ONLY);
         structure.read (io);
         io.close ();
      } catch (SxException e)  {
         e.print ();
         SX_EXIT;
      }
   } else  {
      SxParser parser;
      SxConstPtr<SxSymbolTable> tree;
      tree = parser.read (inFile, "std/structure.std");
      structure = SxAtomicStructure (&*tree);
      structure.readElements (&*tree);
      if (withMovable)  {
         SxPtr<SxStickyFilter> stickyFilter;
         try {
            stickyFilter
               = SxPtr<SxStickyFilter>::create(tree->getGroup ("structure"));
         } catch (const SxException &e)  {
            e.print ();
            SX_EXIT;
         }
         // --- validate the sticky filter
         SX_CHECK (structure.cell.symGroupPtr);
         const SxSymGroup &S = *structure.cell.symGroupPtr;
         SxVector<int> equivalentIdx (structure.getNAtoms());
         equivalentIdx.set (0);
         int           nSym = S.getSize ();
         cout << "| Validating sticky filter ....\n\n";
         for (int iSym = 0; iSym < nSym; iSym++) {
            structure.isEqual (S(iSym) ^ structure, &equivalentIdx);
            stickyFilter->validate (S(iSym).rot, equivalentIdx);
         }
         structure.atomInfo->meta.attach (SxAtomicStructure::StickyFilter,
                                          stickyFilter->getStickyArray ());
      }
   }
   return structure;
}

void writeStructure (const SxAtomicStructure &structure,
                     const SxString &outFile,
                     bool outSxb)
{
   if (outSxb)  {
      structure.write (outFile);
   } else if (outFile.getSize () > 0) {
      structure.fprint (SxTextIO (outFile).getFp ());
   } else {
      structure.fprint (stdout);
   }
}

SxStructOut::SxStructOut (SxCLI &cli, bool sxbFile)
{
   outFile = cli.option ("-o","filename",
                         "output file name (screen otherwise)")
             .toString ("");

   outSxb = cli.option ("--outsxb", "output file is binary "
                        "rather than S/PHI/nX input file").toBool ();
   cli.last ().defaultValue = "default: same format as input file";
   if (sxbFile) outSxb = true;
   if (outSxb && outFile.getSize () == 0)
      outFile = "structure.sxb";
}

SxPtr<SxAOBasis> aoBasisFromInput (const SxParser::Table &table,
                                   const SxString &basisFile,
                                   const SxGkBasis &gkBasis)
{
   const SxAtomicStructure &structure = gkBasis.getTau ();
   SxArray<SxVector<double> > pseudoRad;
   SxArray<SxArray<SxVector<double> > > pseudoPsi;

   SxConstPtr<SxRadBasis> radBasisPtr;
   SxPtr<SxOverlapBase> SPtr;
   SxPtr<SxAOBasis> aoBasisPtr;

   // Normconserving or PAW Potential ?
   if (table->containsGroup("pseudoPot"))   {
      SPtr = SxPtr<SxPWOverlap>::create ();

      if (basisFile.getSize () == 0)  {
         SxPtr<SxPseudoPot> potPtr = SxPtr<SxPseudoPot>::create (&*table);
         radBasisPtr = radBasisPtr.create (potPtr->rad, potPtr->logDr);
         pseudoRad = potPtr->rad;
         pseudoPsi = potPtr->getPseudoPsi ();
         aoBasisPtr = aoBasisPtr.create (gkBasis, *radBasisPtr, pseudoPsi,
                                         SPtr);
      } else {
         SxParser aoParser;
         SxConstPtr<SxSymbolTable> aoTable
            = aoParser.read(basisFile, "std/aoBasis.std");
         SxSymbolTable *aoGroup = aoTable->getGroup("AOBasis");
         SxAtomicOrbitals orbs;
         orbs.setup(aoGroup);
         radBasisPtr = orbs.getRadBasisPtr ();
         aoBasisPtr = aoBasisPtr.create (gkBasis, *radBasisPtr,
                                         orbs.getMuSet ());
      }
   } else if (table->containsGroup("pawPot"))   {
      SxPtr<SxPAWPot> potPtr =SxPtr<SxPAWPot>::create (&*table);
      SxPtr<SxPartialWaveBasis> pBasis
         = SxPtr<SxPartialWaveBasis>::create (potPtr, structure);
      pBasis->createProjBasis (gkBasis);
      SPtr = SxPtr<SxPAWOverlap>::create (pBasis);
      if (basisFile.getSize () == 0)  {
         radBasisPtr = radBasisPtr.create (potPtr->rad, potPtr->logDr);
         pseudoRad = potPtr->rad;
         pseudoPsi = potPtr->getPhiPS ();
         aoBasisPtr = aoBasisPtr.create (gkBasis, *radBasisPtr, pseudoPsi,
                                         SPtr);
      } else {
         SxParser aoParser;
         SxConstPtr<SxSymbolTable> aoTable = aoParser.read(basisFile);
         SxSymbolTable *aoGroup = aoTable->getGroup("AOBasis");
         SxAtomicOrbitals orbs;
         orbs.setup(aoGroup);
         radBasisPtr = orbs.getRadBasisPtr ();
         aoBasisPtr = aoBasisPtr.create (gkBasis, *radBasisPtr,
                                         orbs.getMuSet (), SPtr);
      }
   } else   {
      cout << "No known Potential Group found !" << endl;
      SX_QUIT;
   }
   return aoBasisPtr;
}

