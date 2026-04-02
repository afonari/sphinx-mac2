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
// Authors: Lars Ismer (ismer@fhi-berlin.mpg.de)
//          Michael Ashton (ashton@mpie.de)

#include <SxVDW.h>
#include <SxElemDB.h>
#include <d3Data.h>
#include <SxNeighbors.h>
#include <SxTimer.h>

SxVDW::SxVDW () {}

SxVDW::SxVDW (const SxAtomicStructure &t, const SxSymbolTable *table)
{
   const SxElemDB &elemDB = SxElemDB::getElemDB ();

   int nAtoms = t.nTlAtoms;

   update (t);
   speciesNum     .resize (nAtoms);
   polarizability .resize (nAtoms);
   C6D2           .resize (nAtoms);
   vdwRadiusD2    .resize (nAtoms);
   C6TS           .resize (nAtoms);
   vdwRadiusTS    .resize (nAtoms);
   covalentRadius .resize (nAtoms);
   effectiveVolume.resize (nAtoms);
   d3Coordination .resize (nAtoms);
   dist           .resize (nAtoms);
   neighbors      .resize (nAtoms);

   SxSpeciesData speciesData = SxSpeciesData(table->topLevel());

   for (int i = 0; i < nAtoms; i++) {
      const SxString &species = speciesData.chemName(t.getISpecies (i));
      // atomic number in the periodic table
      int Z = elemDB.getAtomicNumber(species);
      speciesNum(i) = Z;
      polarizability(i) = elemDB.getPolarizability(Z);

      // It's not in their manual, but VASP
      // uses different C6 values for D2 and TS.
      // So I guess we will too.
      C6D2(i) = elemDB.getC6D2(Z);
      vdwRadiusD2(i) = elemDB.getVdwRadiusD2(Z);
      C6TS(i) = elemDB.getC6TS(Z);
      vdwRadiusTS(i) = elemDB.getVdwRadiusTS(Z);

      effectiveVolume(i) = 1.;

      // only relevant for D3 method
      d3Coordination(i) = 0.;
      covalentRadius(i) = elemDB.getCovalentRadius(Z);
      // from Grimme's D3 paper, metallic CN's get overestimated
      // so he reduces their covalent radii by 10%.
      SxArray<SxString> nonMetals = {
            "H", "He", "B", "C", "N", "O", "F",
            "Ne", "Si", "P", "S", "Cl", "Ar", "Ge",
            "As", "Se", "Br", "Kr", "Sb", "Te",
            "I", "Xe", "Po", "At", "Ra"
      };
      if (! nonMetals.contains(species)) {
         covalentRadius(i) *= 0.9;
      }
   }

   // Read in vdW group from input file
   SxSymbolTable *vdwGroup = table -> getGroup ("vdwCorrection");
   SxString methodString = vdwGroup -> get("method") -> toString ();
   if (methodString == "D2") {
      method = Vdw_D2;
   } else if (methodString == "D3") {
      method = Vdw_D3;
   } else if (methodString == "TS") {
      method = Vdw_TS;
   } else if (methodString == "D3+TS") {
      method = Vdw_D3_TS;
   } else if (methodString == "TSI") {
      method = Vdw_TSI;
   } else if (methodString == "D3+TSI") {
      method = Vdw_D3_TSI;
   } else {
      cout << "Unknown van-der-Waals method '" << methodString << "'." << endl;
      SX_QUIT;
   }

   if (vdwGroup -> contains("combinationRule")) {
      SxString combiRule = vdwGroup -> get("combinationRule") -> toString ();
      if (combiRule == "GB")  {
         combinationRule = Vdw_GouldBucko;
      } else if (combiRule == "Tang")  {
         combinationRule = Vdw_Tang;
      } else {
         cout << "Unknown combination rule '" << combiRule << "'." << endl;
         SX_QUIT;
      }
   } else {
      combinationRule = Vdw_Default;
   }

   if (vdwGroup -> contains("damping")) {
      SxString dampingString = vdwGroup -> get("damping") -> toString ();
      if (dampingString == "BJ")  {
         damping = Vdw_BJ;
      } else if (dampingString == "Fermi") {
         damping = Vdw_Fermi;
      } else {
         cout << "Unknown damping rule '" << dampingString << "'." << endl;
         SX_QUIT;
      }
   } else {
      damping = Vdw_Fermi;
   }
}

SxVDW::~SxVDW ()
{
   // empty
}

void SxVDW::update (const SxAtomicStructure &newcoord)
{
   coord = newcoord;
   updateNeighbors ();
   if (method == Vdw_D3) {
      updateD3Coordination ();
      cout << "d3Coordination: " << d3Coordination << endl;
   }
}

void SxVDW::update (const SxArray<double> &effectiveV)
{
   SX_CHECK (effectiveV.getSize () > 0, effectiveV.getSize ());
   effectiveVolume = effectiveV;
}

enum VdwTimers { UpdateNeighbors, VdwCompute, SearchNeighbors };

SX_REGISTER_TIMERS(VdwTimers)
{
   regTimer (UpdateNeighbors, "VDW neighbors");
   regTimer (SearchNeighbors, "VDW n search");
   regTimer (VdwCompute, "VDW compute");
}

void SxVDW::updateNeighbors ()
{
   SX_CLOCK (UpdateNeighbors);
   /*
      Create an array for each atom of the indices
      of all its neighbors within the cutoff
    */

   int nAtoms = coord.getNAtoms ();
   unitVectors.resize (nAtoms);
   dist.resize (nAtoms);
   neighbors.resize (nAtoms);
   nPairs = 0;

   double cutoff = getVdwCutoff ();

   SxGrid grid (coord, 10);
   SxNeighbors nn;
   int mode = SxNeighbors::StoreRel | SxNeighbors::StoreDistSqr;
   for (int i = 0; i < nAtoms; i++) {
      SX_START_TIMER(SearchNeighbors);
      nn.compute (grid, coord, coord(i), cutoff, mode);
      SX_STOP_TIMER(SearchNeighbors);
      int nN = nn.getSize ();
      neighbors(i).resize (nN);
      dist(i).resize (nN);
      unitVectors(i).resize (nN);
      SX_LOOP(in)  {
         neighbors(i)(in) = nn.idx(in);
         double d = sqrt(nn.distSqr(in));
         dist(i)(in) = d;
         unitVectors(i)(in) = nn.relPositions.ref(in) / d;
      }
      nPairs += nN;
   }

}

void SxVDW::updateD3Coordination () {
   /*
      vdW-D3 updates C6 and C8 values based on the local
      coordination number (CN) of each atom. These CN
      values are not necessarily integers.

      NOTE: Constants k1 and k2 are optimized for PBE
    */

   double R, k1, k2, covRad, denominator;
   for (int i = 0; i < neighbors.getSize (); i++) {
      d3Coordination(i) = 0.;
      for (int j = 0; j < neighbors(i).getSize (); j++) {
         R = dist(i)(j);
         int neighj = neighbors(i)(j);
         k1 = 16.;
         k2 = 4./3.;
         covRad = covalentRadius(i) + covalentRadius(neighj);
         denominator = 1 + exp(-k1 * ( k2 * covRad / R - 1));
         d3Coordination(i) += 1. / denominator;
      }
   }
}

double SxVDW::getDampingFunction (double R, double Rij, double *deriv)
{
   /*
    A damping function is required to smoothly turn vdW
    "off" for very close distances dominated by covalent
    interactions. Most vdW methods just use Fermi for
    this.
    */
   SX_CHECK (deriv);

   // --- Fermi damping with cutoff (set by default to~ 100 Bohr)
   double vdwCutoff = getVdwCutoff ();
   if (R >= vdwCutoff) {
      *deriv = 0.;
      return 0.;
   }

   double s6 = (method == Vdw_D2) ? 0.75 : 1.0;
   double d = 20.0;
   double sR = (method == Vdw_D2) ? 1.0 : 0.94;
   double scale = 1.0 / (sR * Rij);
   double chi = exp(-d * (R * scale - 1.0));
   double chi1 = 1./(1. + chi);

   double fd = s6 * chi1;

   //*deriv = s6 * d * scale * chi * chi1 * chi1;
   *deriv = fd * d * scale * chi * chi1;

   return fd;
}

void SxVDW::compute () {
   /*
    This is the main function in SxVDW.
    Here we update the only externally-used
    attributes:
    `totalEnergy` (double)
    `forces` (SxAtomicStructure)
    */
   SX_CLOCK(VdwCompute);
   double R, Rij, C6ij, fd, fdPrime;
   forces = coord.getNewStr ();

   // Reset vdW energy to 0
   double e6 = 0.;
   double e8 = 0.;
   totalEnergy = 0.;
   bool useD3 = (method == Vdw_D3 || method == Vdw_D3_TS);
   for (int i = 0; i < coord.getNAtoms (); i++) {
      // Reset force array for atom i to 0

      forces.ref(i).set (0.);
      for (int j = 0; j < neighbors(i).getSize (); j++) {

         int neighj = neighbors(i)(j);
         R = dist(i)(j);
         Rij = getRij (i, neighj);
         C6ij = getC6ij (i, neighj);

         fd = getDampingFunction (R, Rij, &fdPrime);

         double R6 = R*R*R; R6 = R6 * R6;
         if (useD3) {
            double C8ij = 3.*C6ij*r2r4(speciesNum(i))*r2r4(speciesNum(neighj));
            double R0ij = sqrt(C8ij/C6ij);
            double a1 = 0.4289;
            double a2 = 4.4407;
            double s8 = 0.7875;

             // Default damping in the original D3 paper:
             //fd = 1./(1.+6.*::pow(R/(1.217*::pow(C8ij/C6ij, 0.5)), -14.));

             // Becke-Johnson style damping
             double Q = a1 * R0ij + a2;
             double Q2 = Q*Q, Q6 = Q2*Q2*Q2;
             fd = R6 / (R6 + Q6);
             double R8 = R6*R*R;
             double fd8 = (s8*R8) / (R8 + Q6*Q*Q);

             e8 += -0.5*fd8*C8ij/R8;
             totalEnergy += -0.5*fd8*C8ij/R8;
             //fd8prime = ; // not sure if r8 terms are really necessary for forces
             cout << "Missing D3 forces" << endl;
             // CF 2019-12-19: need consistent implementation of fd and fdPrime !
             // could be put into getDampingFunction (make "if (D3)" there...)
             SX_EXIT;
         }
         double C6R6 = C6ij/R6;
         double dE = -0.5*fd*C6R6;
         e6 += dE;
         totalEnergy += dE;

         // Update forces
         //double derivative = fdPrime * C6ij/::pow(R, 6.) -
         //                    6. * fd * C6ij/::pow(R, 7.);
         //double derivative = 6. * fd * C6ij/::pow(R, 7.) -
         //   fdPrime * C6ij/::pow(R, 6.);
         double derivative = 6. * fd * C6R6/R - fdPrime * C6R6;
         forces.ref(i) += derivative * unitVectors(i)(j);
      }
   }
}

double SxVDW::getRij (int atom1, int atom2)
{
   /*
    Return the combined vdw radius of atom1 and atom2.

    If the TS correction is chosen, this value is modified
    according to each atom's Hirshfeld effective volume.
    */

   double Ri, Rj;
   if (method == Vdw_TS) {
      //Ri = vdwRadiusTS(atom1) * ::pow(effectiveVolume(atom1), 1./3.);
      //Rj = vdwRadiusTS(atom2) * ::pow(effectiveVolume(atom2), 1./3.);
      Ri = vdwRadiusTS(atom1) * cbrt(effectiveVolume(atom1));
      Rj = vdwRadiusTS(atom2) * cbrt(effectiveVolume(atom2));
   } else {
      //Ri = vdwRadiusD2(atom1) * ::pow(effectiveVolume(atom1), 1./3.);
      //Rj = vdwRadiusD2(atom2) * ::pow(effectiveVolume(atom2), 1./3.);
      Ri = vdwRadiusD2(atom1) * cbrt(effectiveVolume(atom1));
      Rj = vdwRadiusD2(atom2) * cbrt(effectiveVolume(atom2));
   }
   double Rij = Ri + Rj;
   return Rij;
}

double SxVDW::getC6ij (int atom1, int atom2)
{
    /*
    Return the combined C6 dispersion coefficient between atom1 and atom2.

    If the TS correction is chosen, this value is modified
    according to each atom's Hirshfeld effective volume.

    If the D3 correction is chosen, this value is based on
    the local coordination of each atom, combined with a
    database of lookup C6 values (see d3Data.h).

    */

   double C6i, C6j, C6ij, alphai, alphaj;

   if (method == Vdw_TS) {
      C6i = C6TS(atom1) * sqr(effectiveVolume(atom1));
      C6j = C6TS(atom2) * sqr(effectiveVolume(atom2));
      alphai = polarizability(atom1) * effectiveVolume(atom1);
      alphaj = polarizability(atom2) * effectiveVolume(atom2);
      // from Tang,K. T. Phys. Rev. 1969, 177, 108
      C6ij = (2 * C6i*C6j /
            ( C6j * alphai / alphaj + C6i * alphaj / alphai ) );
   } else if (method == Vdw_D2) {
      C6i = C6D2(atom1);
      C6j = C6D2(atom2);
      alphai = polarizability(atom1);
      alphaj = polarizability(atom2);
      if (combinationRule == Vdw_Tang) {
         C6ij = (2 * C6i*C6j /
                ( C6j * alphai / alphaj + C6i * alphaj / alphai ) );
      } else if (combinationRule == Vdw_GouldBucko) {
         // from Gould & Bucko (https://arxiv.org/pdf/1604.02751.pdf) page 18
         double alphaij = sqrt( alphai * alphaj);
         C6ij = 1.43 * ::pow (alphaij, 1.45);
      }
      else  C6ij = sqrt(C6i * C6j); // Default
   } else if (method == Vdw_D3 || method == Vdw_D3_TS) {
       ssize_t nIMols = refMolecules(speciesNum(atom1)-1).getSize ();
       ssize_t nJMols = refMolecules(speciesNum(atom2)-1).getSize ();
       double numerator = 0;
       double denominator = 0;
       for (ssize_t im=0; im < nIMols; im++) {
          for (ssize_t jm=0; jm < nJMols; jm++) {
             int iMol = refMolecules(speciesNum(atom1)-1)(im);
             double iCN = refCoordination(iMol);
             int jMol = refMolecules(speciesNum(atom2)-1)(jm);
             double jCN = refCoordination(jMol);
             double L = exp ( -4 * (sqr(( iCN - d3Coordination(atom1) ))
                      + sqr(( jCN - d3Coordination(atom2) )) ) );
             numerator += refC6(iMol)(jMol) * L;
             denominator += L;
          }
       }
       C6ij = numerator / denominator;
       if (method == Vdw_D3_TS) {
          C6ij *= sqr(effectiveVolume(atom1) * effectiveVolume(atom2));
       }
    } else {
       SX_EXIT;
    }

    return C6ij;
}

