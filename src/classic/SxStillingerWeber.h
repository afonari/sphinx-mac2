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

#ifndef _SX_STILLINGER_WEBER_H_
#define _SX_STILLINGER_WEBER_H_

#include <SxClassic.h>
#include <SxPotential.h>

/** \brief Stillinger-Weber empirical potential

    \b SxStillingerWeber = S/PHI/nX Stillinger-Weber Empirical Potential

    ...

    \author Lars Ismer, ismer@fhi-berlin.mpg.de */
class SX_EXPORT_CLASSIC SxStillingerWeber : public SxPotential
{
   public:

      SxStillingerWeber ();
      SxStillingerWeber (const SxAtomicStructure &, const  SxSymbolTable *);
      virtual ~SxStillingerWeber ();

      virtual bool isRegistered (const SxSymbolTable *) const;
      virtual void execute (const SxSymbolTable *, bool calc=true);
      
      virtual SxAtomicStructure getForces (const SxAtomicStructure &,
                                           const SxSymbolTable * =NULL);

      virtual SxSpeciesData getSpeciesData () const;
      virtual PrecEnergy getEnergy () const;
	
      SxArray<SxVector3<double> > coord;

	  SxMatrix3<double> superCell;
	  SxArray<double>  directEnergyContrib;
	  SxArray<double>  indirectEnergyContrib;

	  SxArray<SxList<int> > neighbours;
	   SxArray<SxList<int> > sNeighbours;

	  SxArray<SxList<int> >  borderType;
	  SxArray<SxList<int> > supercellIds;
	  SxArray<SxList<int> > sSupercellIds;


	  SxArray<SxList<int> > borderCrossers;
	  SxArray<SxVector3<double> > directForces;
	  SxArray<SxVector3<double> > indirectForces;
	  SxArray<SxVector3<double> > Forces;

	  SxArray<SxList<double> > dist;
	  SxArray<SxList<double> > expTerm;
	  
	SxSpeciesData speciesData;


	  bool output;
	  int nAtoms;
	  double scalingFactor;
	  double totalEnergy;

	  
	  double getDist(int, int, int);
	  
	  void resize(int, double, SxMatrix3<double>);
	  void update(SxArray<SxVector3<double> >);
	  void updateDistances ();
	  void updateNeighbours ();
	  void updateSecondNeighbours ();

	  void updateBondType ();
	  void updateCoord (SxArray<SxVector3<double> >);
	  void updateBorderCrossers ();
	  
	  void updateExponentialTerms (); 
	  
	  SxVector3<double>  getLatticeVec(int);
	  SxVector3<double>  getIndirectCenterForceOnAtom (int);
	  SxVector3<double>  getDirectForceOnAtom (int);
	  void updateDirectForces (); 
	  void updateIndirectForces (); 
	  SxArray<SxVector3<double> > getForces ();
	  SxArray<SxVector3<double> > getNumericalForces (double);

	  double getParam (SxString, int, int); 

	  void printParameters ();
	  void printNeighbours ();
	  
	  void printSecondNeighbours ();
	  void printCoord ();


};

#endif /* _SX_STILLINGER_WEBER_H_ */
