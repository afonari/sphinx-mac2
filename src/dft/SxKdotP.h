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
#ifndef _SX_KDOTP_H_
#define _SX_KDOTP_H_

#include <SxRBasis.h>
#include <SxBinIO.h>
#include <SxVector3.h>
#include <SxRho.h>
#include <SxHamiltonian.h>
#include <SxFermi.h>
#include <SxTypes.h>

/** \brief SxKdotP n band Hamiltonian

    \b SxKdotP = SFHIngX implementation of k dot p theory
    for general, user-defined n-band Hamiltonians

    \author Oliver Marquardt, marquardt@mpie.de
  */
class SX_EXPORT_DFT SxKdotP : public SxHamiltonian,
                   public SxRho
{
   public:
   
   SxKdotP ();
   /// Constructor
   SxKdotP (const SxPWSet &,
               const RhoR &);
   
   enum Preconditioner { Payne, Arias };

   /// Destructor
   virtual ~SxKdotP () {/*empty*/}

   static const SxSymbolTable * getHamiltonianGroup (const SxSymbolTable *);

   const SxPWSet *wavesPtr;
   
   virtual void read (const SxSymbolTable *); // --- read input parameter
   virtual void writeUpdatedParameters (); // --- write updated parameter files after fit
   virtual void update (const SxFermi &) { /* empty */ }
   virtual void compute (const SxFermi &, bool, bool) { /* empty */ }
   virtual SxRho &getRho ();
   virtual void normalizeRho ();

   virtual void computeRho (const SxPsiSet &, const SxFermi &); 
   virtual PrecEnergy getEnergy (const SxPsiSet &, const SxFermi &);
   virtual void writeRho (const SxString &) const;
   virtual void readRho (const SxBinIO &);
   virtual PrecEnergy getETrial ();
   virtual void printEnergies () const;
   virtual double getEnergy () const { return 0.; }

 
   PsiG operator* (const PsiRef &);
   /// \brief Apply NxN Hamiltonian
   // hook in for SxHamiltonian
   virtual PsiG apply (const PsiRef &psi) const
   {
      // note: can operator* be turned into const?
      return const_cast<SxKdotP*>(this)->SxKdotP::operator* (psi);
   }
 
   PsiR zero, one;
   /// contains j-th psi component while evaluating tree(i,j) in R-space
   PsiRef psiR;
   /// contains j-th psi component while evaluating tree(i,j) in G-space
   PsiRef psiGcomp;
   /// contains 1st derivative of j-th psi component while evaluating tree(i,j)
   SxArray<PsiR>  fdCache;
   /// contains 2nd derivative of j-th psi component while evaluating tree(i,j)
   SxArray<PsiRef> sdCache;
   /// returns 1st derivative of j-th psi component while evaluating tree(i,j)
   const PsiR& fdC(int i);
   /// returns 2nd derivative of j-th psi component while evaluating tree(i,j)
   const PsiRef& sdC(int i, int j);

   SxVector<PrecCoeffG> gZero, gOne;
   const SxGBasis *gBasisPtr;
   SxArray<SxMeshR> materials;
   SxArray<SxMeshR> eIJ; // strain fields
   SxMeshR vP, vExt, outMesh, chargePotential;


   // --- for bandstructure calculations
   SxVector3<double> bsStart, bsEnd, wgt;
   SxVector3<int> plLayer;
   SxVector<double> bandPriority, kPriority, kPriorities;
   SxList<double> kPriorityWeights, kPriorityFWHMs;
   SxArray<SxVector3<double> > bsPts; // list of all points of bandstructure
   bool moreIO, accurateInterfaces;

   //void formDerivatives ();
   bool derivatives, kDerivFirstRound; // belong to formDerivatives()
   /** \brief
      used to calculate <psiI|H|psiJ> for optical momentum matrix
      */    
//   double psiIHPsiJ (const PsiRef &, const PsiRef &);
   void writeMeshASCII (const SxString &name, const PsiR &data) const;

   void showBandstructure(SxString, SxVector3<int>);
   void potentialLandscape(SxString, SxVector3<double>);
   /** \brief Evaluate parse tree for preconditioner.
     @param col which columnn within Hamiltonian matrix
     @param row which row within Hamiltonian matrix
     @param pos position within sub-expression list
     */
   PsiRef evaluateTreePrecond(int col, int row, int pos);

   /// Evaluate partial tree that yields a constant value
   SxComplex16 evaluateTreeNumber (int col, int row, int pos) const;

   /// Evaluate tree to generate bytecode
   void evaluateTreeByteCode(int col, int row, int pos,
                             SxStack<char> &code) const;
   /// Generate bytecode for Hamiltonian matrix
   SxArray2<SxArray<char> > generateByteCode () const;

   /// Hamiltonian in bytecode
   SxArray2<SxArray<char> > byteCodeHam;

   /// Evaluate byte code
   PsiRef evaluateByteCode (const SxArray<char> &code);

   /** \brief Evaluate byte code for band structure
       @param code the bytecode
       @param kVec k-vector
       @param strain vector of 6 strain components
       @param vPol value of polarization potential
       @param vExtern value of external potential
       @param vCharge value of charge potential
       @param paramVals vector of parameter values (nParam)
       @return Hamiltonian matrix element according to bytecode
     */
   SxComplex16 evaluateByteCode (const SxArray<char> &code,
                                 const SxVector3<double> &kVec,
                                 const SxVector<double> &strain,
                                 double vPol, double vExtern, double vCharge,
                                 const SxVecRef<double> &paramVals,
                                 SxStack<SxComplex16> &stack) const;
   /// Print what byte code does
   void printByteCode (const SxArray<char> &code) const;

   /** \brief
     Transform Hamiltonian element to parse tree
     */
   void buildTree(int, int, SxString, SxString);
   bool containsOperator(int, int, int) const;
   bool isRealNumber(int, int, int) const;
   void printElement(int, int, int) const;
   void countOperatorMultiplications(SxString);
   void determineDerivatives(int);

   // resolve expression and return position of right branch in tree list
   int resolveExpr(SxString, int, int);
   static int resolveOpKey(const SxString &);
   /** \brief
     A specified parameter can be separately extracted in an .sxb file,
     e.g. for viewing using the Phinax viewer.
     */
   int nOpMult, iMult;
   SxString outPar;
   SxArray<int> opKey;
   SxArray<PsiR> kIPar, kJPar, par, kKPar;

   /** \brief
     Routines that mark elements as operator or potential.
     To make sure that no element in the parse tree is
     counted twice (i.e. double multiplication of the
     element with the wave function, unnecessary marks
     are removed after parsing the tree.
     */
   void correctOperators(int, int);
   void removeOpMarks(int, int, int);
   /** \brief
     Determine wether element is an operator or a potential term
     */
   int whatIsElement(int, int, int);
   int firstOutsideBracket(SxString, char);
   /** \brief
     The parse tree: expression, left branch, right branch
     */
   SxArray<SxArray<SxList<SxString> > >expression;
   SxArray<SxArray<SxList<int> > > leftPtr, rightPtr;
   /** \brief
     routines search for unknown expressions in the hamiltonian
     and replace these if given previously in the input file
     */
   SxString containsUnknown(SxList<SxString>);
   SxString replaceUnknown(SxString, SxString);
   /** \brief create Hamiltonian matrix for band structure plot
       @note This uses the fast bytecode algorithm
     */
   SxVector<SxComplex16> hMatrixBS(const SxVector3<double> &kVec,
                                   const SxVector<double> &strain,
                                   double vPol, double vExtern, double vCharge,
                                   const SxVecRef<double> &params,
                                   SxStack<SxComplex16> &stack) const;

   SxArray<SxMeshR> pMem; // for speed-optimized calculation: parameters in memory
   bool speedOpt;
   SxVecRef<double> parameters (int iParam);
   SxVector<double> matParam; // all parameters of all materials
   SxVector<double> matParamOriginal; // all parameters of all materials
   SxVector<double> fitParam; // all fit ranges of all materials
   SxVector<double> bowParam; // all bowings of involved materials
   SxVector<double> nonzeros;
   int nMat, nParam, nComp, precMaterial;
   SxList<SxString> matNames; // names of involved materials
   SxList<SxString> paramNames; // names of involved parameters
   SxList<SxString> chargeFiles; // names of external charge files
   SxList<SxVector<double>> kPrioritiesList, bandPrioritiesList;
   SxList<int> fitParams, nFitSteps, nFitPoints, BZSamplingList; // parameters that request fit, nIterations for fit, Points in Sobol sequence
   SxList<SxString> fitFileNames; // file names of input band structures for fit

   SxList<double> eVBList, eCBList, gapPriceList;

   SxVector<PrecCoeffG> D;
   
   PrecEnergy eTrial;
   inline const SxPWSet &getWavesRef () const {
      SX_CHECK (wavesPtr);
      return *wavesPtr;
   }
 
   virtual void set (const SxPsiSet &, const SxFermi &);
   /** \brief The general n-band k.p-Hamiltonian
    */
   PsiR hamNxN (const PsiRef &);
   PsiR returnAccurate(int);
   
   void fitBandstructure (SxVector<double>, int, int);
   
   virtual SxVector<double> preconditioner (const PsiRef &psi) const
   {
      return preconditioner (psi, Payne);
   }

   SxVector<double> preconditioner (const PsiRef &,
                                    Preconditioner type) const;
};

#endif /* _SX_KDOTP_H_ */
