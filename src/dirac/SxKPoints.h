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

#ifndef _SX_K_POINTS_H_
#define _SX_K_POINTS_H_

#include <SxVector3.h>
#include <SxVector.h>
#include <SxArray.h>
#include <SxList.h>
#include <SxDiracLib.h>
#include <SxAtomicStructure.h>
#include <SxSymbolTable.h>

/** \brief k-point container class

    \b SxClass = SPHInX k-points

    \todo Documentation, notably of Monkhorst-Pack

    \author Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_DIRAC SxKPoints
{
   public:

      /** List of k vectors in absolute coordinates */
      SxArray<SxVector3<PrecG> >            kVec;
      /** List of weights of all the k points */
      SxVector<PrecWeights>                 weights;
      /** Unsymmetrized Monkhorst-Pack k-points */
      SxAtomicStructure                      kVecMP;
      /** Weights of unsymmetrized Monkhorst-Pack k-points */
      SxVector<PrecWeights>                 weightsMP;
      /** k-vectors that enter Monkhorst-Pack folding */
      SxArray<SxVector3<PrecG> >            kVecPrimitive;
      /** Weights of primitve k-points */
      SxVector<PrecWeights>                 weightsPrimitive;
      /** Labels of special k-points (may have 0 size if no labels present) */
      SxList<SxString>                       kLabels;
      /** \brief Number of k vectors */
      int                                    nk;
      /** \brief Return number of k vectors */
      inline int getNk () const { return nk; }
      /** Monkhorst-Pack folding (along reciprocal basis vectors) */
      SxVector3<int>                         foldingMP;
      /** uses the additional time-reversal symmetry (default value:
          true */
      bool                                   useInvSymmetry;

      /** Empty constructor */
      SxKPoints () : nk(-1), useInvSymmetry(true) { /* empty */}
      /** Copy constructor */
      SxKPoints (const SxKPoints &in);
      /** Constructor */
      SxKPoints (const SxCell &cell,
                 const SxList<SxVector3<PrecG> > &kp_,
                 const SxList<double> &weights_,
                 const SxVector3<int> &folding_,
                 bool useInvSymmetry_=true);
      /** Constructor */
      SxKPoints (const SxCell &cell,
                 const SxSymbolTable *,
                 bool useInvSymmetry_=true);

      /// Assignment operator
      SxKPoints& operator= (const SxKPoints &) = default;
    
      /** \brief Read from symbol table
          @param recCell reciprocal Cell
          @param table the symbol table
          @param pointSpec "k" or "q"
          @param kPointId  - name of the MonkhorstPack-type group
          @param kPointsId - name of the bandstructure-type group
        */
      void read (const SxCell &recCell, const SxSymbolTable *table,
                 const SxString &pointSpec, const SxString &kPointId, 
                 const SxString &kPointsId);

      /** Destructor */
      ~SxKPoints ();

      /** Returns a copy of the specified k point. */
      SxVector3<PrecG> getK (int ik) const;

      /**
        \brief Initialize from binary file.
        \param io  netcdf binary file, in mode SxBinIO::BINARY_READ_ONLY

        */
      void read (const SxBinIO &io);

   protected:

      /** \brief Internal initialization routine
          @param recCell reciprocal cell
          @param kp_      initial k-point before folding
          @param weights_ initial weights before folding
          @param folding_ Monkhorst-Pack folding parameters
        */
      void init (const SxCell &recCell,
                 const SxList<SxVector3<PrecG> > &kp_,
                 const SxList<double> &weights_,
                 const SxVector3<int> &folding_);

      /** \brief Generate Monkhorst-Pack mesh and symmetrize it.
          \param recCell Reciprocal cell with symmetries (i.e. from structure)
        */
      void generateMonkPack (const SxCell &recCell);

      /** \brief Symmetry-reduce the Monkhorst-Pack points
          \param recCell Reciprocal cell with symmetries (i.e. from structure)
          
          \note Whether time reversal (k/-k) symmetry is used depends on 
                #useInvSymmetry.
          \todo The search routine scales quadratically with the number
                of k-points and becomes slow for large numbers. Use
                a better algorithm (e.g. subcell partitioning) for the
                kS <-> kVecMP mapping.
        
        */
      void symmetryReduce (const SxCell &recCell);

      /** \brief Check if k-point is high-symmetry point
          \note
          In practice, ensure that this is Gamma (0,0,0) or lies
          on Wigner-Seitz boundary within 1e-4.
          Check that there is at least one symmetry that leaves k-point
          unchanged.
        */
      static bool isHighSymPoint (const SxCell &recCell, const Coord &k);

};

#endif /* _SX_K_POINTS_H_ */
