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
#ifndef _SX_TB_HAM_SOLVER_H_
#define _SX_TB_HAM_SOLVER_H_

#include <SxPotential.h>
#include <SxAtomicStructure.h>
#include <SxTBSpecies.h>
#include <SxTBAtomAtom.h>
#include <SxSymMatrix.h>
#include <SxEigensystem.h>
#include <SxRhoMixer.h>
#include <SxGrid.h>

namespace Timer {
   enum tbHamSolverTimer {
      initialization,
      tbElMinim,
      diagonalization,
      Forces,
      mullikenRho,
      gammaMatrix,
      ewald,
      hamOvl,
      hamMod
   };
}

SX_REGISTER_TIMERS (Timer::tbHamSolverTimer)
{
   using namespace Timer;
   regTimer (initialization, "Initialization");
   regTimer (tbElMinim,      "Electronic Loop");
   regTimer (diagonalization,"Diagonalization");
   regTimer (Forces,         "Forces");
   regTimer (mullikenRho,    "Mulliken Charges");
   regTimer (gammaMatrix,    "Gamma Matrix");
   regTimer (ewald,          "Ewald Summation");
   regTimer (hamOvl,         "H & S Mat. Calc.");
   regTimer (hamMod,         "H Mat. SCF-Corr.");
}

/** \brief Tight-binding Hamiltonian solver

    \b SxTBHamSolver = S/PHI/nX Tight-Binding Hamiltonian Solver

    \author Hazem Abu-Farsakh, h.farsakh@mpie.de */

class SX_EXPORT_DFT SxTBHamSolver : public SxPotential
{
   public:

      SxTBHamSolver (const SxAtomicStructure &str,
                     const SxSymbolTable *table);
      virtual ~SxTBHamSolver ();

      virtual bool isRegistered (const SxSymbolTable *) const;
      virtual void execute (const SxSymbolTable *, bool calc=true);
      virtual PrecEnergy getEnergy () const;
      virtual SxAtomicStructure getForces (const SxAtomicStructure &tau,
                                           const SxSymbolTable * =NULL);
      virtual SxSpeciesData getSpeciesData () const;
   
   protected:

//----------------------------------------------------------------------------
//    1) functions
//----------------------------------------------------------------------------

      void tbInitialGuess (const SxSymbolTable *, bool calc=true);
      void diagonalize    (const SxSymbolTable *, bool calc=true);
      
//----------------------------------------------------------------------------
      /**@name Gamma k-point calculations */
//----------------------------------------------------------------------------
      //@{ 
      void diagonalizeGamma ();
     
      void calcForcesGamma ();

      /** \brief calculates the non-SCF Hamiltonian (HG0) and overlap (SG)
                 matrices  */
      void calcHG0AndSG ();
   
      /** \brief calculates the correction to the non-SCF Hamiltonian matrix 
          \param rho      vector of mulliken charges on atoms.:(iAtom) */
      void calcHG1 (SxArray<SxVector<double> >  &rhoIn);

      /** \brief  diagonalization process (generalized). 
          \return the eigensystem
          \param  HIn Hamiltonian matrix 
          \param  SIn Overlap matrix */
      SxSymEigensystem<double> diag (const SxVector<double> &HIn,
                                     const SxVector<double> &SIn) const;

      /** \brief set the fermi function */ 
      void setFermi (const SxVector<double> &eigenvalues);

      /** \brief calculates and returns the density matrix */
      SxVector<double> getPDensMat () const;

      /** \brief calculates and returns the energy density matrix */
      SxVector<double> getPErgDensMat () const;

      /** \brief  Mulliken atomic charges. : (iAtom)(iOrb)
          \return Mulliken atomic charges. : (iAtom)(iOrb) */
      SxArray<SxVector<double> > getMullikenRhoG (); //:iAtom, iOrb
      //@}

//----------------------------------------------------------------------------

      /** calculates and returnes the additional contribution 
        to the electronic energy in the SC calculations (eCorrectoin) */
      double getECorrection () const;
      
      /** calculates and returns additional contribution to forces 
          due to self-consistent scheme, for non-periodic systems.:(iAtom)*/
      SxArray<SxVector3<double> > getDRhoForcesNonPeriodic ();

      /** calculates and returns additional contribution to forces 
          due to self-consistent scheme, for periodic systems.:(iAtom)*/
      SxArray<SxVector3<double> > getDRhoForcesPeriodic ();

//----------------------------------------------------------------------------
      /**@name Many k-points calculations */
//----------------------------------------------------------------------------
      //@{ 
      void diagonalizeCmplx ();

      void calcForcesCmplx ();

      /** \brief calculates the non-SCF Hamiltonian (H0) and overlap (S)
                 matrices */
      void calcH0AndS ();

      /** \brief calculates the correction to the non-SCF Hamiltonian matrix 
          \param rho vector of mulliken charges on atoms.:(iAtom) */
      void calcH1 (SxArray<SxVector<double> > &rhoIn);

      /** \brief  diagonalization process (generalized). 
          \return the eigensystem
          \param  HIn Hamiltonian matrix 
          \param  SIn Overlap matrix */
      SxSymEigensystem<SxComplex16> diag (const SxVector<SxComplex16> &HIn,
                                          const SxVector<SxComplex16> &SIn) const;
      /** \brief set the fermi function */ 
      void setFermi ();
      
      /** calculates and returns the density matrix */
      SxVector<SxComplex16> getPDensMat (int ik) const;

      /** calculates and returns the energy density matrix */
      SxVector<SxComplex16> getPErgDensMat (int ik) const;

      /** \brief  Mulliken atomic charges. : (iAtom)(iOrb)
          \return Mulliken atomic charges. : (iAtom)(iOrb) */
      SxArray<SxVector<double> > getMullikenRho (); //:iAtom, iOrb
      //@}

//----------------------------------------------------------------------------
      /**@name General */
//----------------------------------------------------------------------------
      //@{ 
      /** \brief initialize an indexing system for the Hamiltonian and
                 the overlap, such that each atom has an index 
          \sa idx  */ 
      void initIdx (const SxAtomicStructure &str);
      
      /** \brief  calculates the repulsive energy part 
          \return the repulsive energy */
      double getRepulsiveEnergy ();

      /** \brief  isolated atoms charges. :(iAtom)(iOrb) 
          \return isolated atoms charges */
      SxArray<SxVector<double> > getAtomicRho ();
      //@}
      
//----------------------------------------------------------------------------
      /**@name calculating gamma values for non-periodic structures */
//----------------------------------------------------------------------------
      //@{ 
      /** \brief  calculates gamma values (as a matrix) for non-periodic systems.
          \return gamma values, used for calculating the modification
                  to the non-SCF Hamiltonian in the SCF calculations 
          \sa     gammaSCF */
      SxVector<double> getGammaMatNonPeriodic () const;

      /** \return a value called gamma, used for calculating the modification
           to the non-SCF Hamiltonian in the SCF scheme*/
      double gammaSCF (int iA, int jA, int iOrb, int jOrb, 
                       double distance) const;

      /** a function (which is short range w.r.t. distance between atoms)
        used in the calculation of Gamma value */
      double shortRange (double tA, double tB, double distance) const; 

      /** a function used in the calculation of the short-range function, 
          \sa shortRange */
      double K (double tA, double tB, double distance) const; 
      
      /** \brief  calculates the derivative of gamma values w.r.t. distance,
                  as a matrix, and for non-periodic systems. 
          \return derivative of gamma values, used in calculating the 
                  SCF forces. */
      SxVector<double> getDGammaMatNonPeriodic () const;

      /** \return derivative of gammaSCF */
      double dGammaSCF (int iA, int jA, int iOrb, int jOrb, 
                        double distance) const; 

      /** calculates and returnes derivative of shortRange with respect 
         to distance */
      double dShortRange (double tA, double tB, double distance) const;

      /** derivative of K w.r.t. distance */
      double dK (double tA, double tB, double distance) const;
      //@}

//----------------------------------------------------------------------------
      /**@name calculating gamma values for periodic systems */
//----------------------------------------------------------------------------
      //@{     
      /** initialization of a G-space sphere up to certain cutoff to produce 
         a list of reciprocal space translation vectors */
      void initGSphere();

      /** initialization of an R-space sphere up to certain cutoff to produce 
         a list of real space translation vectors */
      void initRSphere();

      /** \brief calculates R and G vectors that are enough for convergence 
                 of ewald sum parts. \sa longRangeSum */
      void initRadii();
      
      /** \brief calculates the best eta value for a better convergence 
                 in both real and reciprocal space summations of Ewald sum */
      double getEta ();

      /** \brief a function used by getEta ()
          \sa    getEta () */
      double diffReciprocalReal (double sR, double sG, double etaVal);

      /** \brief  calculates gamma values for periodic systems, as a matrix. 
          \return gamma values, used for calculating the modification
                  to the non-SCF Hamiltonian in the SCF calculations 
          \sa     gammaSCF */
      SxVector<double> getGammaMatPeriodic () const;

      /** long range sum ( SUM_R {1/|r - R|} ) for periodic systems, 
          using Ewald method */
      double longRangeSum (int iA, int jA, 
                           SxVector3<double> distVec) const;

      /** sum of short range funtion( SUM_R {short range function} ), 
          for periodic systems */
      double shortRangeSum (int iA, int jA, int iOrb, int jOrb, 
                            SxVector3<double> distVec) const;

      /** \brief  calculates the derivative of gamma values w.r.t. distance,
                  as a matrix, for periodic systems. 
          \return derivative of gamma values, used in calculating the 
                  SCF forces. */
      SxArray<SxVector<double> > getDGammaMatPeriodic () const;

      /**  sum of short range funtion derivative, for periodic systems */
      SxVector3<double> dShortRangeSum (int iA, int jA, int iOrb, int jOrb,
                                        SxVector3<double> distVec) const;

      /** long range sum ( SUM_R {1/|r - R|} ) derivative (for periodic 
          systems), using Ewald method */
      SxVector3<double> dLongRangeSum (int iA, int jA,
                                       SxVector3<double> distVec) const;
      //@}
   
//----------------------------------------------------------------------------
      /**@name Service routines */
//----------------------------------------------------------------------------
      //@{ 
      /** \brief  extracts a sub-matrix from a matrix.
          \return a sub-matrix (from the big one). 
          \param  pMat   the complete "big" matrix
          \param  index1 row index of the first entry in the sub-matrix
          \param  index2 column index of the first entry in the sub-matrix
          \param  nOrb1  number of rows-1    in the sub-matrix
          \param  nOrb2  number of columns-1 in the sub-matrix.
          i.e. row index of the last entry in the sub-matrix = index1 + nOrb1 
          & column index of the last entry in the sub-matrix = index2 + nOrb2 */
      SxVector<double> extractMat (const SxVector<double> &pMat,
                                   int index1, int index2,
                                   int nOrb1,  int nOrb2) const
      {
         return pMat.blockRef (index1, index2, nOrb1, nOrb2);
      }

     /** extracts a sub-matrix from a complex matrix \sa extractMat */
      SxVector<SxComplex16> extractMat (const SxVector<SxComplex16> &pMat,
                                      int index1, int index2, 
                                      int nOrb1,  int nOrb2) const
      {
         return pMat.blockRef (index1, index2, nOrb1, nOrb2);
      }

      
      /** utility to put rho in a vector form */
      SxVector<double> putRhoInVec
                       (SxArray<SxVector<double> > &rhoArray) const;

      /** utility to put rho from vector form to array of vectors,
          :(iAtoms) */
      SxArray<SxVector<double> > putRhoInArrayOfVec
                                (SxVector<double> &rhoVec) const;
      
      /** check if the matrix is Hermitian */
      bool isHermitian (const SxVector<SxComplex16> &matIn) const;

      /** \brief return the minimum of three numbers */
      inline double minimum (double n1, double n2, double n3)
      { double n4 = (n1 < n2) ? n1 : n2; return (n4 < n3) ? n4 : n3; }
      
      /** utility to print mulliken charges into a file */
      void printMullikenRho() const;
      //@}

//----------------------------------------------------------------------------
//    2) variables
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
      /**@name General */
//----------------------------------------------------------------------------
      //@{ 
      /** Atom-atom Tight binding objects. :(iSpecies)(jSpecies) */
      SxArray<SxArray<SxTBAtomAtom> >  tbAtomAtom;
      
      SxAtomicStructure structure;
      SxTBSpecies       data; 
      SxKPoints         kPoints;
      SxFermi           fermi;
      double            ekt;

      /** \brief  gird for calculating neighbors */
      SxGrid            grid;

      /** \brief number of species */
      int nSpecies; 

      /** \brief total number of atoms */
      int nTlAtoms; 

      /** \brief total number of electrons */
      int nTlElect;

      /** \brief number of atoms of each species. :(iSpecies) */
      SxVector<int> nAtoms;

      /** \brief number of electrons of each species. :(iSpecies)*/
      SxVector<int> nElect;
      
      /** \brief true for Gamma-point calculations */
      bool isGammaPoint; 

      /** \brief true for periodic systems */
      bool periodic; 

      /** \brief true will claculate total energy using zero repulsive energy

        This is used if the repulsive energy coefficients are not available
        in the slater-koster files.  */
      bool withoutERepulsive;

      /** \brief the sum of the number of orbitals of all atoms */
      int lOrbMax; 
      
      /** \brief index for each atom, used to put the atom - atom matrix 
                 (sub matrix) into the complete matrix. :(iAtom)  */
      SxVector<int> idx;

      /** \brief band structure energy */
      double eBand;
      
      /** \brief repulsive energy */
      double eRepulsive;
      
      /** \brief total energy */
      double  eTotal;

      /** an additional contribution to the electronic energy 
          (in case of SCF calculations only) */
      double eCorrection;

      /** \brief cutoff radius (for neighbors only) */
      double cutoff;

      /** \brief true for self-consistent calculations*/ 
      bool   isSCF;

      double rhoMixing;
      double dCharge;
      int    maxSteps;
      int    printSteps;
      int    nMixingSteps;

      /** \brief charges mixer for SCF loop */ 
      SxRhoMixer::MixerType mixerType;

      /** \brief gamma values stored in a matrix. 
          \sa    getGammaMatNonPeriodic, getGammaMatPeriodic, gammaSCF */
      SxVector<double> gammaMatrix;

      /** \brief true if R & G spaces spheres are initialized */ 
      bool   spheresInit;

      /** \brief cutoff radius for real space part of ewald sum */
      double rEwald;

      /** \brief cutoff radius for reciprocal space part of ewald sum */
      double gEwald;

      /** \brief parameter in Ewald sum. 
        
          eta value affects only the convergence in real and 
          reciprocal spaces: larg eta -> slow conv. in G-space, 
                            small eta -> slow conv. in R-space */
      double eta;

      /** \brief List of real space translation vectors */
      SxList<SxVector3<Real8> > rLat;

      /** \brief List of reciprocal space translation vectors */
      SxList<SxVector3<Real8> > gLat;
      //@}

//----------------------------------------------------------------------------
      /**@name Gamma k-point calculations */
//----------------------------------------------------------------------------
      //@{ 
      /** \brief the non-SCF Hamiltonian matrix */
      SxVector<double> HG0;

      /** \brief the overlap matrix */
      SxVector<double> SG;

      /** \brief eigen system  */
      SxSymEigensystem<double> eigG;
      
      /** \brief correction to the non-SCF Hamiltonian matrix 
        (in SCF calculations) */
      SxVector<double> HG1;
      
      /** \brief part of HG1
          \sa    HG1  */
      SxVector<double> HG1Part1 ;
      //@}

//----------------------------------------------------------------------------
      /**@name Many k-points calculations */
//----------------------------------------------------------------------------
      //@{ 
      /** \brief the Hamiltonian matrices. :(ik) */
      SxArray<SxVector<SxComplex16> > H0;
      
      /** \brief the Overlap matrices. :(ik) */
      SxArray<SxVector<SxComplex16> > S;

      /** \brief eigen systems */
      SxArray<SxSymEigensystem<SxComplex16> > eig;

      /** \brief correction to the non-SCF Hamiltonian matrix 
        (in SCF calculations). :(ik) */
      SxArray<SxVector<SxComplex16> > H1;
      
      /** \brief part of H1
          \sa    H1  */
      SxVector<double> H1Part1;
      //@}

//----------------------------------------------------------------------------
      /**@name Forces */
//----------------------------------------------------------------------------
      //@{ 
      /** \brief Hamiltonian contribution to the forces. :(iAtom)*/
      SxArray<SxVector3<double> > hamForces;
      
      /** \brief repulsive forces. :(iAtom) */
      SxArray<SxVector3<double> > repulsiveForces;

      /** \brief total forces. :(iAtom) */
      SxAtomicStructure fTotal;  
      //@}

//----------------------------------------------------------------------------
      /**@name Charges  */
//----------------------------------------------------------------------------
      //@{ 
      /** \brief mulliken atomic charges. :(iAtom)(iOrb) */
      SxArray<SxVector<double> > mullikenRho;

      /** \brief isolated atoms charges. :(iAtom)(iOrb) */
      SxArray<SxVector<double> > atomicRho;

      /** memory for charges */
      SxArray<SxVector<double> > rhoMemory;
      /** false causes to get charges from previous memory */
      bool firstTime;
      //@}

};

#endif /* _SX_TB_HAM_SOLVER_H_ */

