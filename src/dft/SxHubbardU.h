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

#ifndef _SX_HUBBARD_U_H_
#define _SX_HUBBARD_U_H_

#include <SxDFT.h>
#include <SxVector.h>
#include <SxBlockDensityMatrix.h>
//#include <SxHubbardMO.h>
class SxPAWPot;
class SxAtomicStructure;
class SxSymbolTable;
class SxHubbardMO;


/** \brief Hubbard U extension in the "rotationally invariant" formulation

 * Reference 1: SPHInX Hubbard U implementation notes
 * Reference 2: M. Cococcioni, S. de Gironcoli, Phys. Rev. B 71, 035105 (2005)

    \author C. Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxHubbardU
{
   public:
      /// Constructor
      SxHubbardU ();

      /// The effective U at each site
      SxArray<double> Ueff;

      /// energy
      double energy;

      /// double counting correction energy
      double eDoubleCounting;

      /// Control debug output
      bool verbose;

      /// Set energy etc. to zero (before computeIncr)
      void setZero ();
      /** \brief Incremental compute
          @param iSite  which site
          @pram  nij    occupation matrix

          @return \f$\partial E/\partial n_ij
        */
      SxVector<SxComplex16> computeIncr(int iSite,
                                      const SxVecRef<SxComplex16> &nij);

      /** \brief Site bias for fitting/constraints
        */
      SxArray<double> alphaBias;

      /// The occupation matrix at each site
      SxBlockDensityMatrix siteOccupation;

      /// Get total site occupation for one site
      double getSiteTotal (ssize_t iSite) const;

      /// Get number of sites
      inline ssize_t getSize ()  const
      {
         return Ueff.getSize ();
      }
      /// Resize (keeping current values)
      void resize (ssize_t nSite);

      /** Abstract class for conversions between PAW projectors
          and Hubbard orbitals */
      class HubbardPAWMapper {
         public:
            /// l quantum number for the orbital type
            int l;

            /// Map the PAW density matrix Dij to the occupancy matrix
            virtual SxVector<SxComplex16>
            mapPAWtoHubbard (const SxVecRef<double> &Dij) = 0;

            /** \brief Add the Hubbard contribution to PAW projector Hamiltonian
               @param Hij the Hubbard Hamiltonian in orbital space
               @param Vij the PAW onsite Hamiltonian
               */
            virtual void addHubbardToPAW (const SxVecRef<SxComplex16> &Hij,
                                          SxVecRef<double> &Vij) = 0;

            /// Destructor
            virtual ~HubbardPAWMapper () {}
      };

      /** \brief atomic site map

          This container class maps atoms to Hubbard U sites.
          This map must be set up from the outside.

        */
      class AtomicSite {
         public:
            /// Site ID in SxHubbardU
            int iSite;
            /// Mapping from PAW projectors to Hubbard orbitals
            SxPtr<HubbardPAWMapper> map;
      };

      /// Atomic site map, listing for each atom the Hubbard U sites (:iTlAtom,iSite (at atom)
      SxArray<SxArray<AtomicSite> > atomicSiteMap;

      /** \brief List of MO sites (iTypeMO)

          This list contains the types of MO sites. Each type can consist
          of several sites
        */
      SxArray<SxPtr<SxHubbardMO> > moSite;

      /** \brief Compute contribution of atom iTl
          @param iTl    which atom
          @param Dij    PAW density matrix
          @param fFull  value of full occupation: 1 for spin-polarized
                        calculations, 2 otherwise
          @return contribution to Vij
        */
      SxVector<double> computeAtom(int iTl, const SxVecRef<double> &Dij,
                                     double fFull);

      /// Compute MO sites
      void computeMO (const SxArray<SxPtr<SxBlockDensityMatrix> > &Pij,
                      const SxAtomicStructure &structure);

      /// Read atomic sites from symbol table
      void read (const SxSymbolTable *table,
                 const SxPtr<SxPAWPot> &pawPotPtr,
                 const SxAtomicStructure &structure);

      /// Synchronize occupations
      void syncMPI ();

      /// Write occupations
      void writeOccupations(const SxString &fileName) const;
};

#endif /* _SX_HUBBARD_U_H_ */
