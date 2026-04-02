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

#include <SxCLI.h>
#include <SxCell.h>
#include <SxFFT1d.h>
#include <SxGBasis.h>
#include <SxPreconditioner.h>
#include <SxRBasis.h>
#include <SxNeighbors.h>
#include <SxTimer.h>
#include <SxDefectAlignUtil.h>
#include <SxTextIO.h>
#include <SxIonicCoreFit.h>

using namespace SxDefectAlignUtil;

double rhoRec (double g2, double beta2, double gamma2, double x)
{
   return         x / sqr(1. + gamma2 * g2)
          + (1. - x) * exp( - 0.25 * beta2 * g2);
}

double rhoRecLimit0 (double beta2, double gamma2, double x)
{
   // rhoRec(G->0) -> 1 + rhoRecLimit0 * G^2
   return       x  * (-2. * gamma2)
          + (1.-x) * (-0.25 * beta2);
}


/** \brief Description of a single model charge
  */
class SxCharge {
   public:
      /// Charge
      double charge;
      /// Position
      Coord pos;
      /// Beta of Gaussian
      double beta2;
      /// Relative weight of exponential tail
      double expNorm;
      /// Decay of exponential
      double gamma2;

      /** \brief Constructor
          \param Q     charge
          \param where position
          \param betaIn broadening for Gaussian
          \param gammaIn decay of exponential tail (if any)
          \param expNorm relative weight of exponential tail
        */
      SxCharge (double Q = 0., const Coord &where = Coord(0.,0.,0.),
                double beta = 1.,
                double gamma = -1., double expNormIn = 0.)
         : charge(Q), pos(where), beta2(beta * beta),
           expNorm(expNormIn), gamma2(gamma * gamma)
      {
         // empty
      }

      inline double rhoRec (double g2) const
      {
         return charge * ::rhoRec(g2, beta2, gamma2, expNorm);
      }

      inline SxComplex16 rhoRec (double g2, const Coord &g) const
      {
         return rhoRec(g2) * SxComplex16::phase (-(g ^ pos));
      }

      inline double rhoRecLimit0 () const
      {
         return charge * ::rhoRecLimit0(beta2, gamma2, expNorm);
      }

};


class SxAtomAverages {
   public:
      /// G-basis for defect and reference
      SxGBasis defG, refG;
      /// potential in G space: defect, reference, long-range model
      SxVector<PrecCoeffG> vDefG, vRefG, vLrDef;

      /// Constructor
      SxAtomAverages(const SxAtomicStructure & defStr,
                     const SxAtomicStructure & refStr,
                     const SxMeshR &vDef, const SxMeshR &vRef,
                     const SxMesh3D &defMesh, const SxMesh3D &refMesh,
                     double eCut)
      : defG(defMesh, defStr, eCut), refG(refMesh, refStr, eCut),
        vDefG ( defG | vDef ), vRefG ( refG | vRef ), vLrDef (defG)
      {
         vLrDef.set (0.);
      }

      /// Set average of long-range model
      void setV0 (double V0)
      {
         vLrDef(0) = V0;
      }

   protected:
      /// Calculate potential for one G-vector
      static PrecCoeffG getPot (double g2, const Coord &g,
                                const SxArray<SxCharge> &chargeList,
                                const SxMatrix3<double> &epsTensor)
      {
         SX_CHECK (g2 > 0., g2);
         PrecCoeffG rho = 0.;
         for (auto &it : chargeList)
            rho += it.rhoRec (g2) * SxComplex16::phase (-(g ^ it.pos));
         return FOUR_PI / (g ^ epsTensor ^ g) * rho;
      }

   public:
      /// Add potential from several charges
      void addPot (const SxArray<SxCharge> &charges,
                   const SxMatrix3<double> eps)
      {
         for (int ig = 1; ig < defG.ng; ++ig)
            vLrDef(ig) += getPot (defG.g2(ig), defG.getG (ig), charges, eps);
      }

      /// Compute atomic averages from potential and print to file
      void evaluate (const Coord &center, double betaAvg, bool field);
};

void SxAtomAverages::evaluate (const Coord &center, double betaAvg, bool field)
{
   const SxAtomicStructure &defStr = defG.getTau ();
   const SxAtomicStructure &refStr = refG.getTau ();
   double betaAvg2 = betaAvg * betaAvg;

   // evaluate potential at atom positions (with Gaussian) for defect cell
   double normD = 1./sqrt(defStr.cell.volume);
   SxVector<PrecCoeffG> gaussDef = exp((-0.5*betaAvg2) * defG.g2);
   gaussDef *= normD;
   SxVector<PrecCoeffG> gaussRef = exp((-0.5*betaAvg2) * refG.g2);
   gaussRef /= sqrt(refStr.cell.volume);

   ofstream outFile("vAtoms.dat");
   SxGrid grid (refStr, 10);
   for (int is = 0; is < refStr.getNSpecies (); ++is)  {
      for (int defIa = 0; defIa < defStr.getNAtoms(is); ++defIa)  {
         // map to defect cell
         Coord atomPos = defStr(is,defIa);
         int refIdx = refStr.find(atomPos, grid);
         int refIa = -1;
         if (refIdx >= 0 && refStr.getISpecies (refIdx, &refIa) == is)  {
            SxVector<PrecCoeffG> defPhase = gaussDef * defG.getPhaseFactors (is, defIa);
            SxVector<PrecCoeffG> refPhase = gaussRef * refG.getPhaseFactors (is, refIa);
            double vDefAtom = dot(defPhase, vDefG).re;
            double vLrAtom = dot(defPhase, vLrDef).re * normD ;
            double vRefAtom = dot(refPhase, vRefG).re;
               outFile << defStr.cell.getMapped (defStr(is,defIa) - center,
                     SxCell::WignerSeitz).norm ()
               << '\t' << vLrAtom * HA2EV
               << '\t' << (vDefAtom - vRefAtom) * HA2EV
               << '\t' << (vDefAtom - vRefAtom - vLrAtom) * HA2EV;
            if (field)  {
               for (int idir = 0; idir < 3; ++idir)
                  outFile << '\t' <<
                     (-dot(defG.getPhaseFactors (is, defIa),
                           vLrDef * defG.gVec.colRef(idir)).im) * normD;
            }
            SX_LOOP(iDir) outFile << '\t' << atomPos(iDir);
            outFile << endl;
         }
      }
      outFile << endl;
   }
}

class SxIsolatedEnergy {
   public:
      // cutoff energy
      double eCut;
      // integration step
      double gStep;

      /// Constructor
      SxIsolatedEnergy (double cutoff, double step)
         : eCut(cutoff), gStep (step)
      { /* empty */ }

      /** \brief Compute unscreened isotropic interaction for isolated case
          \param chargeList list of charges
          \param dr         distances
      */
      double compute (const SxArray<SxCharge> &chargeList,
                      const SxArray<double> &dr) const;
      /** \brief Compute unscreened isotropic interaction for isolated case
                 for a pair of charge distributions
          \param chargeList1 list of charges
          \param chargeList2 list of charges
      */
      double compute (const SxArray<SxCharge> &charges1,
                      const SxArray<SxCharge> &charges2) const;

      void setupAniso (const SxMatrix3<double> &epsTensor,
                       double dielecConstant);
      /// Compute anisotropy of charge-charge interactions
      double computeAniso (const SxArray<SxCharge> &chargeList,
                           const SxArray<double> &dr) const;
      /// Compute anisotropy of charge-charge interactions
      double computeAniso (const SxArray<SxCharge> &charges1,
                           const SxArray<SxCharge> &charges2) const;
   protected:
      SxVector<double> angInt;
};

double SxIsolatedEnergy::compute (const SxArray<SxCharge> &chargeList,
                                  const SxArray<double> &dr) const
{
   int nCharges = (int)chargeList.getSize ();
   SxArray<double> rho(nCharges);
   double eIso = 0., dE = 0.;
   double wSimpson = 4.;
   // G=0
   double qTot = 0.;
   for (int ic = 0; ic < nCharges; ++ic)
      qTot += chargeList(ic).charge;
   eIso += qTot * qTot;
   // G > 0
   for (double g = gStep; g*g < eCut || wSimpson > 3.; g += gStep)  {
      double g2 = g * g;
      dE = 0.;
      // single charges
      for (int ic = 0; ic < nCharges; ++ic)  {
         rho(ic) = chargeList(ic).rhoRec (g2);
         dE += sqr(rho(ic));
      }
      // charge-charge interactions
      for (int ic = 0, ij = 0; ic < nCharges; ++ic)  {
         for (int jc = ic + 1; jc < nCharges; ++jc, ++ij)  {
            dE += 2. * rho(ic) * rho(jc) * SxYlm::jsb (0, g * dr(ij));
         }
      }

      eIso += wSimpson * dE;
      wSimpson = (wSimpson > 3.) ? 2. : 4.;
   }
   eIso -= dE; // final Simpson weight is 1 rather than 2.
   return eIso * gStep / (3. * PI) ;
}

double SxIsolatedEnergy::compute (const SxArray<SxCharge> &charges1,
                                  const SxArray<SxCharge> &charges2) const
{
   SxArray2<double> dr((int)charges1.getSize (), (int)charges2.getSize ());
   SX_LOOP2(jc,ic)
      dr(ic, jc) = (charges1(ic).pos - charges2(jc).pos).norm ();
   double eIso = 0., dE = 0.;
   // G=0
   {
      double qTot1 = 0., qTot2 = 0.;
      for (auto it : charges1) qTot1 += it.charge;
      for (auto it : charges2) qTot2 += it.charge;
      eIso += qTot1 * qTot2; // Simpson weight = 1.
   }
   // G > 0
   double wSimpson = 4.;
   SxArray<double> rho1(charges1.getSize ());
   for (double g = gStep; g*g < eCut || wSimpson > 3.; g += gStep)  {
      double g2 = g * g;
      dE = 0.;
      // precompute rho1
      SX_LOOP(ic) rho1(ic) = charges1(ic).rhoRec (g2);

      SX_LOOP2(jc, ic)
         dE += rho1(ic) * charges2(jc).rhoRec (g2)
             * SxYlm::jsb (0, g * dr(ic, jc));

      eIso += wSimpson * dE;
      wSimpson = (wSimpson > 3.) ? 2. : 4.;
   }
   eIso -= dE; // final Simpson weight is 1 rather than 2.
   return 2. * eIso * gStep / (3. * PI) ;
}

void SxIsolatedEnergy::setupAniso (const SxMatrix3<double> &epsTensor,
                                   double dielecConstant)
{
   const int lMax = 12;
   SxSphereGrid grid(SxSphereGrid::Grid_230);
   grid.setupYlm (lMax);

   angInt.resize (sqr(lMax+1));
   angInt.set (0.);
   // angular integral in reciprocal space
   // int dOmega Ylm(g) / (g ^ eps ^ g)
   for (int i = 0; i < grid.getSize (); ++i)  {
      Coord xyz = grid.getXyz((int)i);
      double epsInv = 1. / (xyz ^ epsTensor ^ xyz);
      for (int l = 2; l <= lMax; l += 2)  {
         for (int m = -l; m <= l; ++m)  {
            int lm = SxYlm::combineLm(l,m);
            angInt(lm) += grid.weights(i) * epsInv * grid.Ylm(i,lm);
         }
      }
   }
   // apply normalization factor for YlmReal (see below) here
   for (int l = 2; l <= lMax; l += 2)  {
      double sign = (l % 4) == 0 ? 1. : -1.; // i^l
      for (int m = -l; m <= l; ++m)  {
         angInt(SxYlm::combineLm(l,m)) *= SxYlm::getYlmNormFactor (l,m)
                                        * sign;
      }
   }
   angInt *= FOUR_PI * dielecConstant;
   // cout << angInt << endl;
}

/// Compute anisotropy of charge-charge interactions
double SxIsolatedEnergy::computeAniso (const SxArray<SxCharge> &chargeList,
                                       const SxArray<double> &dr) const
{
   // G > 0
   const int lMax = 12;
   int nCharges = (int)chargeList.getSize ();
   SxVector<double> YlmReal(sqr(lMax+1));
   SxVector<double> integral(lMax);
   double eIso = 0.;
   // sum over charge pairs
   for (int ic = 0, ij = 0; ic < nCharges; ++ic)  {
      for (int jc = ic + 1; jc < nCharges; ++jc, ++ij)  {
         integral.set (0.);
         double wSimpson = 4.;
         for (double g = gStep; g*g < eCut || wSimpson > 3.; g += gStep)  {
            double g2 = g * g;
            double rho2 = chargeList(ic).rhoRec(g2)
                        * chargeList(jc).rhoRec(g2);
            for (int l = 0; l < lMax; l+=2)  {
               integral(l) = SxYlm::jsb(l + 2, g * dr(ij)) * rho2;
               integral(l+1) += wSimpson * integral(l);
            }
            wSimpson = (wSimpson > 3.) ? 2. : 4.;
         }
         // get Ylm in real space
         SxYlm::getYlmArray (lMax,
                             chargeList(ic).pos - chargeList(jc).pos,
                             &YlmReal);
         for (int l = 2; l <= lMax; l += 2)  {
            integral(l-1) -= integral(l-2); // final wSimpson = 1, not 2
            for (int m = -l; m <= l; ++m)  {
               int lm = SxYlm::combineLm(l,m);
               eIso += angInt(lm) * YlmReal(lm) * integral(l-1);
            }
         }
      }
   }
   return eIso * 2. * gStep / (3. * PI) ;
}

/// Compute anisotropy of charge-charge interactions
double SxIsolatedEnergy::computeAniso (const SxArray<SxCharge> &charges1,
                                       const SxArray<SxCharge> &charges2) const
{
   const int lMax = 12;
   SxVector<double> YlmReal(sqr(lMax+1));
   SxVector<double> integral(lMax);
   ssize_t ng = 0;
   for (double g = gStep; g*g < eCut; g += gStep) ng++;
   if (ng % 2 == 1) ng--;
   // --- precompute rho1 density
   SxVector<double> rho1 (ng, charges1.getSize ());
   SX_LOOP (ig) {
      double g2 = sqr(double(ig+1) * gStep);
      SX_LOOP(ic) rho1(ig, ic) = charges1(ic).rhoRec (g2);
   }
   double eIso = 0.;
   // sum over charge pairs
   SX_LOOP2(jc,ic)  {
      Coord dr = charges1(ic).pos - charges2(jc).pos;
      double rij = dr.norm ();
      if (rij < 1e-10) continue;
      integral.set (0.);
      SX_LOOP (ig) {
         double wSimpson = (ig.i & 1) ? 2. : 4.;
         double g = double(ig+1) * gStep;
         double g2 = g * g;
         double rho2 = rho1(ig, ic) * charges2(jc).rhoRec(g2);
         for (int l = 0; l < lMax; l+=2)  {
            integral(l) = SxYlm::jsb(l + 2, g * rij) * rho2;
            integral(l+1) += wSimpson * integral(l);
         }
      }
      // get Ylm in real space
      SxYlm::getYlmArray (lMax, dr, &YlmReal);
      for (int l = 2; l <= lMax; l += 2)  {
         integral(l-1) -= integral(l-2); // final wSimpson = 1, not 2
         for (int m = -l; m <= l; ++m)  {
            int lm = SxYlm::combineLm(l,m);
            eIso += angInt(lm) * YlmReal(lm) * integral(l-1);
         }
      }
   }
   return eIso * 2. * gStep / (3. * PI) ;
}

SxMatrix3<double> toTensor (SxCLI::CliArg &option)
{
   SxList<double> epsVals = option.toDoubleList ();
   SxMatrix3<double> epsTensor (0.);
   switch (epsVals.getSize ())  {
      case 1:
         epsTensor(0,0) = epsTensor(1,1) = epsTensor(2,2) = epsVals(0);
         break;
      case 3:
         for (int i = 0; i < 3; ++i) epsTensor(i,i) = epsVals(i);
         break;
      case 6:
         for (int i = 0; i < 3; ++i) epsTensor(i,i) = epsVals(i);
         epsTensor(1,2) = epsTensor(2,1) = epsVals(3);
         epsTensor(0,2) = epsTensor(2,0) = epsVals(4);
         epsTensor(1,0) = epsTensor(0,1) = epsVals(5);
         break;
      case 9:
         epsTensor = SxMatrix3<double>(epsVals);
         break;
      default:
         if (!option.getParent ().error)  {
            cout << "Failed to interpret ";
            option.printMark ();
            cout << " with " << epsVals.getSize ()
                 << " values - must be 1, 3, 6, or 9" << endl;
            option.getParent ().setError ();
         }
   }
   return epsTensor;
}

void printOut (const SxArray<SxCharge> &chargeList)
{
   for (auto it : chargeList)  {
      cout << "Excess Electrons = " << it.charge
           << " located at " << it.pos;
      if (fabs(it.beta2 - 1.) > 1e-10)
         cout << " (beta=" << sqrt(it.beta2) << ")";
      cout << endl;
      if (fabs(it.expNorm) > 1e-10)
         cout << "Exponential tail: norm=" << it.expNorm
              << " gamma=" << sqrt(it.gamma2) << endl;
   }
}

int main (int argc, char** argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   SxCLI cli(argc, argv);
   cli.authors = "C. Freysoldt";
   cli.version ("3.0");

   double eCut, vAlign, dielecConstant, avgWidth;
   SxMatrix3<double> epsTensor(1,0,0,0,1,0,0,0,1);
   SxString refPotFile, defPotFile;

   int sxInput = cli.newGroup ("sxdefectalign input file");
   SxString inFile = cli.option ("-i|--input", "file",
                                 "sx-style sxdefectalign input file")
                     .toString ();
   
   int clInput = cli.newGroup ("cell from command line");
   cli.excludeGroup(sxInput);
   Coord col1 (cli.option ("--axis1","vector", "a1 (units: bohr)")
               .toList3 ()),
         col2 (cli.option ("--axis2","vector", "a2 (units: bohr)")
               .toList3 ()),
         col3 (cli.option ("--axis3","vector", "a3 (units: bohr)")
               .toList3 ());
   
   cli.newGroup ("parameters from command line");
   cli.excludeGroup(sxInput);
   if (!cli.groupAvailable (sxInput))
      cli.currentGroup->status = SxCLI::Required;
   eCut = cli.option ("--ecut", "energy (Ry)", "cutoff energy in Ry")
                 .toDouble ();

   double step = cli.option ("--gstep", "step size", "|G| step size")
                 .toDouble (1e-4, 1e-10);

   bool relative = cli.option ("--relative", "pos is in relative coordinates")
                   .toBool ();
   
   bool printRho
      = cli.option ("--printRho", "print model rho to rhoModel.dat")
        .toBool ();
   
   cli.newGroup ("CSRB screening");
   cli.excludeGroup(sxInput);
   double q2tf = cli.option ("--csrb", "(screening vector)^2", "square of "
         "Thomas-Fermi screening vector for CSRB screening").toDouble (0.); 
   
   cli.newGroup("dielectric specification");
   cli.excludeGroup(sxInput);
   bool useTensor = cli.option ("--tensor", "tensor",
   "dielectric tensor. Specify as\n"
   "eps_xx,eps_yy,eps_zz (3 values) for a diagonal tensor\n"
   "eps_xx,eps_yy,eps_zz,eps_yz,eps_xz,eps_xy (6 values) or ...\n"
   "xx,xy,xz,xy,yy,yz,xz,yz,zz                (9 values) for a full tensor")
                 .exists ();
   cli.last ().optional = true;
   if (useTensor && !cli.error)
      epsTensor = toTensor (cli.last ());

   cli.option ("-e|--eps", "constant", "dielectric constant");
   if (useTensor && cli.last ().exists ())  {
      cout << "Both dielectric constant and tensor are given!" << endl;
      cout << "Dielectric constant will be ignored" << endl;
   }
   dielecConstant = cli.last ().toDouble (1.);


   int vGroup = cli.newGroup ("read potentials");
   cli.excludeGroup(sxInput);
   refPotFile = cli.option ("--vref", "potential file",
                                 "reference potential")
                         .toString ();
   defPotFile = cli.option ("--vdef", "potential file",
                                 "defect potential")
                         .toString ();
   bool defIsHartree = cli.option ("--defInHartree",
                               "defect potential file in Hartree").toBool ();
   bool refIsHartree = cli.option ("--refInHartree",
                               "reference pot file in Hartree").toBool ();
   FileType fileType = getFileType (cli);
   bool printOrigPot
      = cli.option ("--printPot", "print defect and reference potential")
        .toBool ();
   SxString structDefFile
      = cli.option ("--structDef", "SPHInX structure file", 
                    "evaluate potentials at atomic coordinates given in the SPHInX file").toString ("");
   SxString structRefFile
      = cli.option ("--structRef", "SPHInX structure file", 
                    "evaluate potentials at atomic coordinates given in the SPHInX file").toString ("");
   double betaAvg = cli.option ("--atomAverage","length",
                                "Gaussian broadening for atomic-sphere averages").toDouble (1.5, 0.);
   bool atomField
      = cli.option ("--field", "compute long-range fields at atoms")
        .toBool ();

   double fitCoreCut = 3.75;
   int nCorePts = 15, coreOrder = 3;
   bool antiCore = cli.option ("--anticore", "rcut[,nPoints[,order]]",
                               "remove ionic cores from potential").toBool ();
   cli.last ().defaultValue = "Defaults: rcut=3.75 nPoints=15 order=3";
   if (cli.last ().hasValue ())  {
      SxList<SxString> vals = cli.last ().getValue ().tokenize (',');
      try {
         fitCoreCut = vals(0).toDouble ();
         if (vals.getSize () > 1)
            nCorePts = vals(1).toInt ();
         if (vals.getSize () > 2)
            coreOrder = vals(2).toInt ();
      } catch (SxException e)  {
         cout << "Failed to interpret --anticore parameters rcut[,nPoints[,order]] = " << vals(0);
         for (int i = 1; i < vals.getSize (); ++i) cout << ',' << vals(i);
         cli.setError ();
      }
   }

   enum { OldStyle, AllColumn, Separate } outStyle = OldStyle;
   cli.option ("--format", "matrix|separate",
         "output format for vline-eV-a{0,1,2}.dat.\n"
         "'matrix'   - <x> <Vlr> <def-ref> <Vsr>\n"
         "'separate' - <x> <Vlr>\n\n"
         "             <x> <def-ref>\n\n"
         "             <x> <Vsr>\n"
         "'xmgrace'  - <x> <Vlr>\n"
         "             &\n"
         "             <x> <def-ref> <Vsr>\n");
   cli.last ().defaultValue = "default: xmgrace";
   if (cli.last ().toBool () && cli.last ().hasValue ())  {
      SxString style = cli.last ().getValue ();
      if (style == "matrix")
         outStyle = AllColumn;
      else if (style == "separate")
         outStyle = Separate;
      else if (style != "xmgrace")  {
         cout << "Unknown format style '" << style << "'" << endl;
         cli.setError ();
      }
   }

   cli.setGroup (cli.generalGroup);
   avgWidth = cli.option ("--average", "length", "local average (bohr)")
              .toDouble (0., 0.);
   vAlign = cli.option ("-C|--align", "align",
                        "potential alignment constant (eV)")
            .toDouble (0.) / HA2EV;


   cli.setLoopSeparator ("--pos|--center|--vertical-eps", true);
   int chargeGroup = cli.newGroup("model charge specifications");
   cli.excludeGroup(sxInput);

   SxMatrix3<double> eps2Tensor(1,0,0,0,1,0,0,0,1);
   double eps2Inv = 0.;
   bool mixedScreening = false;
   // --- collect charge positions
   SxArray<SxCharge> chargeList, elChargeList;
   while (cli.looping ())  {
      // first capture "sticky" options charge ...
      double Q = cli.option ("-q|--charge", "excess electrons",
            "number of additional excess electrons of the defect: a defect "
            "charge of -1 corresponds to 1 excess electrons, while a defect "
            "charge of +1 corresponds to -1 excess electrons")
                       .toDouble (1.);
      // ... and beta
      double beta = cli.option ("--beta", "length",
                                "Gaussian decay exp(-r^2/beta^2)")
                      .toDouble (1.);
      // --- switch to electronic screening
      if (cli.option ("--vertical-eps", "value(s)",
                      "mixed screening: electron-only dielectric tensor.\n"
                      "1 value  = scalar.\n"
                      "3 values = xx yy zz\n"
                      "6 values = xx yy zz yz xz xy\n"
                      "9 values = xx xy xz yx yy yz zx zy zz\n")
          .sticky (false).exists () && !cli.error)
      {
         if (mixedScreening)  {
            cli.setError ();
            cout << "ERROR: --vertical-eps can be used only once!" << endl;
         }
         if (chargeList.getSize () == 0)  {
            cli.setError ();
            cout << "ERROR: no (ionic) charges given before --vertical-eps"
                 << endl;
         }
         mixedScreening = true;
         eps2Tensor = toTensor (cli.last ());
         continue;
      }
      //cout << cli.arguments(chargeList.getSize ()) << endl;
      SxArray<SxCharge> &chList = mixedScreening ? elChargeList : chargeList;
      chList.resize (chList.getSize () + 1, true);
      SxCharge &current = chList(chList.getSize () - 1);

      cli.setGroup (chargeGroup);
      cli.option ("--pos|--center", "vector", "defect center position");
      cli.last ().sticky (false);
      cli.last ().optional = true;
      if (cli.last ().exists ())  {
         current.pos = Coord(cli.last ().toList3 ());
      }
      current.charge = Q;
      current.beta2 = beta * beta;

      cli.newGroup ("exponential model");
      cli.excludeGroup(sxInput);
      current.gamma2 = cli.option ("--gamma", "length",
                                   "exponential decay exp(-r/gamma)")
                       .sticky (false)
                       .toDouble ();
      current.gamma2 *= current.gamma2;
      current.expNorm = cli.option ("--expnorm", "<0..1>",
                                "relative norm of exponential part")
                       .sticky (false)
                       .toDouble ();
   }

   cli.finalize ();

   ssize_t nCharges = chargeList.getSize ();

   if (useTensor || mixedScreening)  {
      if ((epsTensor - epsTensor.transpose ()).absSqr ().sum () > 1e-12)  {
         cout << "Dielectric tensor must be symmetric" << endl;
         cout << "Tensor is " << epsTensor << endl;
         SX_QUIT;
      }
      if (q2tf > 1e-12)  {
         cout << "CSRB screening cannot be used with dielectric tensors"
              << endl;
         SX_QUIT;
      }
      SxSphereGrid grid(SxSphereGrid::Grid_110);
      double epsInv = 0.;
      SX_LOOP(i) {
         Coord xyz = grid.getXyz((int)i);
         epsInv  += grid.weights(i) / (xyz ^ epsTensor  ^ xyz);
         eps2Inv += grid.weights(i) / (xyz ^ eps2Tensor ^ xyz);
      }
      if (useTensor) dielecConstant = 1./epsInv;
   }

   // --- read cell
   bool printFullVLine = false;
   SxCell defCell, refCell;
   SxMesh3D defMesh, refMesh;
   // vDef and vRef in Hartree
   SxMeshR defPot, refPot;
   SxAtomicStructure structRef, structDef;

   if (cli.groupAvailable (sxInput)) {
      SxParser parser;
      SxConstPtr<SxSymbolTable> table = parser.read (inFile, "std/defectalign.std");
      SxSymbolTable *mainGroup = table->getGroup("DefectAlign");
      
      // --- cell
      SxSymbolTable *subGroup = mainGroup->getGroup("cellGroup");
      if (subGroup->contains("fromVElStat")) 
         cout << "Cells are taken from electrostatic potentials." << endl;
      else if (subGroup->contains("file"))  {
         SxString cellFile = subGroup->get("file")->toString();
         SxParser cellParser;
         SxParser::Table cellTable = cellParser.read (cellFile,"std/structure.std");
         try {
            defCell = SxCell(&*cellTable);
            cout << "cell defect = " << defCell << endl;
         } catch (SxException e)  {
            e.print ();
            SX_EXIT;
         }
      } else if (subGroup->contains("cell"))  {
         defCell = CellMat(subGroup->get("cell")->toList ()).transpose ();
         cout << "cell defect = " << defCell << endl;
      } else {
         cout << "No valid option in cell group!" << endl;
         SX_QUIT;
      }
      
      // --- electrostatic potentials
      SxSymbolTable *vTable = mainGroup->getGroup("vElStat");
      if (vTable->contains("fileType"))  {
         SxString fileMode = vTable->get("fileType")->toString();
         if (fileMode == "VASP") fileType = VASP_LOCPOT;
         else if (fileMode == "SPHINX") fileType = sxb;
         else {
            cout << "Unknown filetype for electrostatic potential." << endl;
            SX_QUIT;
         }
      } else fileType = sxb;
      subGroup = vTable->getGroup("vDefect");
      defPotFile = subGroup->get("file")->toString();
      defPot = getPot (defCell, defMesh, defPotFile, fileType, defIsHartree, &structDef);
      cout << "cell defect = " << defCell << endl;
      subGroup = vTable->getGroup("vReference");
      refPotFile = subGroup->get("file")->toString();
      refPot = getPot (refCell, refMesh, refPotFile, fileType, refIsHartree, &structRef);
      cout << "cell bulk = " << refCell << endl;
      printFullVLine = true;

      // --- parameters
      subGroup = mainGroup->getGroup("parameters");
      eCut = subGroup->get("eCut")->toReal ();

      // --- model charge
      subGroup = mainGroup->getGroup("modelCharge");
      if (subGroup->contains("dielecConstant")) 
         dielecConstant = subGroup->get("dielecConstant")->toReal ();
      else 
         dielecConstant = 1.0;

      chargeList.resize (subGroup->getNItems ("gaussian"));
      int iCharge = 0;
      for (SxSymbolTable *gaussian = subGroup->getGroup ("gaussian");
           gaussian != NULL;
           gaussian = gaussian->nextSibling ("gaussian"))
      {
         chargeList(iCharge).charge = gaussian->get("electrons")->toReal ();
         chargeList(iCharge).pos = Coord(gaussian->get("pos")->toList ());
      }
   } else  {
      if (cli.groupAvailable (vGroup))  {
         printFullVLine = true;
         defPot = getPot (defCell, defMesh, defPotFile, fileType, defIsHartree, &structDef);
         cout << "cell defect = " << defCell << endl;
         refPot = getPot (refCell, refMesh, refPotFile, fileType, refIsHartree, &structRef);
         cout << "cell bulk = " << refCell << endl;
#ifdef SXDA_STRUCTURE_PRINTSX
         // --- print structure in sx format (for debugging)
         structDef.fprint (SxTextIO ("defectStruct.sx").getFp ());
         structRef.fprint (SxTextIO ("refStruct.sx").getFp ());
#endif
      } else if (inFile.getSize () > 0)  {
         SxParser parser;
         SxParser::Table table = parser.read (inFile,"std/structure.std");
         try {
            defCell = SxCell(&*table);
         } catch (SxException e)  {
            e.print ();
            SX_EXIT;
         }
      } else if (cli.groupAvailable (clInput))  {
         defCell = SxCell (col1, col2, col3);
      } else  {
         cout << "No cell provided, specify electrostatic potentials or SPHinX structure/inputfile, or"
            " use command line option!" << endl;
         SX_QUIT;
      }
   }

   SxFFT::quickFFTPlanner ();

   // --- handle relative coordinates (after obtaining defCell)
   if (relative)  {
      for (auto &it : chargeList)
         it.pos = defCell.relToCar (it.pos);
      if (mixedScreening)
         for (auto &it : elChargeList)
            it.pos = defCell.relToCar (it.pos);
   }

   SxPtr<SxRBasis> RDef, RRef;
   if (defMesh.product () != 0)  {
      RDef = RDef.create (defMesh, defCell);
      defPot.setBasis(&*RDef);
   }
   if (refMesh.product () != 0)  {
      RRef = RRef.create (refMesh, refCell);
      refPot.setBasis(&*RRef);
   }

   if (mixedScreening) cout << "Ionic charge:" << endl;
   printOut (chargeList);
   if (mixedScreening) {
      cout << "-------------" << endl << "Electronic charge:" << endl;
      printOut (elChargeList);
   }

   // --- setup atomic structures
   cout << "Atomic structure ";
   if (structDefFile.getSize () > 0 && structRefFile.getSize () > 0)  {
      SxParser parser;
      structRef = SxAtomicStructure (&*parser.read(structRefFile, "std/structure.std"));
      structDef = SxAtomicStructure (&*parser.read(structDefFile, "std/structure.std"));
   } else if (structDef.getNAtoms () == 0 || structRef.getNAtoms () == 0) {
      if (antiCore)  {
         cout << "is missing, but needed for ionic core removal!" << endl;
         SX_QUIT;
      }
      cout << "not ";
   }
   cout << "specified." << endl;

   if (printFullVLine && antiCore)  {
      cout << "Fitting core potentials with rcut=" << fitCoreCut
           << " nPoints=" << nCorePts << " spline order=" << coreOrder
           << "..." << endl;
      SxIonicCoreFit ionFit;
      ionFit.set (nCorePts, fitCoreCut, coreOrder);
      ionFit.addFitData (refPot, refMesh, structRef);
      cout << "elements=" << ionFit.chemNames << endl;
      ionFit.addFitData (defPot, defMesh, structDef);
      cout << "elements=" << ionFit.chemNames << endl;
      ionFit.computeFit (true);
      cout << "Removing cores..." << endl;
      refPot -= ionFit.getFunc (structRef, refMesh);
      defPot -= ionFit.getFunc (structDef, defMesh);
   }

   // --- compute energies
   double eIso = 0.;

   // --- isolated

   if (q2tf < 1e-10)  {
      // --- precompute intercharge distances
      SxArray<double> dr(nCharges * (nCharges - 1) / 2);
      for (int ic = 0, ij = 0; ic < nCharges; ++ic)
         for (int jc = ic + 1; jc < nCharges; ++jc, ++ij)
            dr(ij) = (chargeList(ic).pos - chargeList(jc).pos).norm ();

      // --- unscreened charge self-energy
      SxIsolatedEnergy isoEn (eCut, step);
      eIso = isoEn.compute (chargeList, dr);

      // --- anisotropic corrections to charge-charge terms
      if (useTensor) isoEn.setupAniso (epsTensor, dielecConstant);
      if (useTensor && nCharges > 1)
         eIso += isoEn.computeAniso (chargeList, dr);
      if (useTensor || mixedScreening) eIso /= dielecConstant;

      if (mixedScreening)  {
         // ionic-electronic interaction
         eIso += isoEn.compute (chargeList, elChargeList) / dielecConstant;
         if (useTensor)
            eIso += isoEn.computeAniso (chargeList, elChargeList)
                  / dielecConstant;

         // --- get dr for electronic charges
         ssize_t nCharges2 = elChargeList.getSize ();
         dr.resize (nCharges2 * (nCharges2 - 1) / 2);
         for (int ic = 0, ij = 0; ic < nCharges2; ++ic)
            for (int jc = ic + 1; jc < nCharges2; ++jc, ++ij)
               dr(ij) = (elChargeList(ic).pos - elChargeList(jc).pos).norm ();

         // electronic charge: isotropic self-energy
         eIso += isoEn.compute (elChargeList, dr) * eps2Inv;

         if (nCharges2 > 1)  {
            // anisotropic electronic charge self-energy
            isoEn.setupAniso (eps2Tensor, 1./eps2Inv);
            eIso += isoEn.computeAniso (elChargeList, dr) * eps2Inv;
         }
      }
   } else {
      if (nCharges > 1)  {
         cout << "CSRB for multiple charges not implemented" << endl;
         SX_QUIT;
      }
      // G = 0;
      eIso = sqr(chargeList(0).charge)
           * SxPreconditioner::invCSRB(0.,q2tf,dielecConstant);
      double lastVal = 0.;
      for (double g = step; g*g < eCut; g += 2. * step)  {
         eIso += 4. * sqr(chargeList(0).rhoRec(sqr(g)))
              * SxPreconditioner::invCSRB(sqr(g),q2tf,dielecConstant);
         eIso += 2. * (lastVal = sqr(chargeList(0).rhoRec(sqr(g+step)))
              * SxPreconditioner::invCSRB(sqr(g+step),q2tf,dielecConstant));
      }
      eIso -= lastVal;
      eIso *= step / (3. * PI) ;
   }

   // --- periodic
   if (mixedScreening && !useTensor)  {
      // mixed screening works only for tensorial code
      useTensor = true;
      SX_LOOP(i) epsTensor(i,i) = dielecConstant;
   }
   SxCell recCell = defCell.getReciprocalCell ();

   SxMesh3D mesh = SxGBasis::getMeshSize (eCut, defCell);
   if (defMesh.product () == 0) defMesh = mesh;
   int meshSize = mesh.product (), ng = 0;

   double ePeriodic = 0.;
   for (int i = 1; i < meshSize; ++i)  {
      Coord g = recCell.relToCar(mesh.getMeshVec(i, SxMesh3D::Origin));
      double g2 = g.normSqr ();
      if (g2 < eCut)  {
         SxComplex16 rhoPeriodic = 0.;
         for (auto it : chargeList)  {
            rhoPeriodic += it.rhoRec(g2, g);
         }
         if (q2tf < 1e-10)
            ePeriodic += rhoPeriodic.absSqr () / (g ^ epsTensor ^ g);
         else
            ePeriodic += rhoPeriodic.absSqr () / g2
                       * SxPreconditioner::invCSRB(sqr(g2),q2tf,dielecConstant);
         ng++;
      }
   }
   if (mixedScreening)  {
      for (int i = 1; i < meshSize; ++i)  {
         Coord g = recCell.relToCar(mesh.getMeshVec(i, SxMesh3D::Origin));
         double g2 = g.normSqr ();
         if (g2 < eCut)  {
            SxComplex16 rhoIon = 0., rhoElectron = 0.;
            for (auto it : chargeList)   rhoIon += it.rhoRec(g2, g);
            for (auto it : elChargeList) rhoElectron += it.rhoRec(g2, g);

            // ion-electron interaction
            ePeriodic += 2. * (rhoIon.conj () * rhoElectron).re / (g ^ epsTensor ^ g);
            // electronic charge self-energy
            ePeriodic += rhoElectron.absSqr () / (g ^ eps2Tensor ^ g);
         }
      }
   }
   cout << "ng=" << ng << endl;
   ePeriodic *= TWO_PI / defCell.volume;
   double V0 = 0., totQ = 0.;
   for (auto it : chargeList)  {
      V0 += it.rhoRecLimit0 ();
      totQ += it.charge;
   }
   V0 *= FOUR_PI;
   if (useTensor)
      V0 /= dielecConstant;
   if (mixedScreening)  {
      double V0el = 0., totEl = 0.;
      for (auto it : elChargeList)  {
         V0el += it.rhoRecLimit0 ();
         totEl += it.charge;
      }
      // correction for width of electronic charge (screened with dielecConstant)
      // this accounts for the interaction of the ionic screening response due to
      // the presence of the background for the ionic charge with the electronic
      // finite charge distribution (rather than a point charge)
      ePeriodic += FOUR_PI * V0el * totQ * (1./ dielecConstant - eps2Inv)
                 / defCell.volume;
      V0el *= FOUR_PI * eps2Inv;
      V0 += V0el;
      totQ += totEl;
   }

   if (q2tf >= 1e-10)  {
      // switch to "unscreened"
      eIso *= dielecConstant;
      ePeriodic *= dielecConstant;
   }
   ePeriodic += V0 * totQ / defCell.volume;

   if (structDef.getNAtoms () > 0 && structRef.getNAtoms () > 0)  {
      structRef.epsEqual = 0.5;
      SxAtomAverages avg(structDef, structRef, defPot, refPot, defMesh,
                         refMesh, eCut);
      avg.addPot (chargeList, epsTensor);
      if (mixedScreening)
         avg.addPot (elChargeList, eps2Tensor);
      avg.setV0 (V0);
      if (epsTensor == SxMatrix3<double> (1.))
         avg.vLrDef /= dielecConstant;
      avg.evaluate (chargeList(0).pos, betaAvg, atomField);
   }

   // --- compute potential
   SxTextIO rhoFile;
   if (printRho) rhoFile.open ("rhoModel.dat");
   SxArray<Coord> posRel(nCharges);
   SX_LOOP(iCharge)
      posRel(iCharge) = defCell.carToRel (chargeList(iCharge).pos);
   for (int idir = 0; idir < 3; idir++)  {
      int nz = defMesh(idir);
      double dz = defCell.basis(idir).norm () / nz;
      SxFFT1d fft1d (SxFFT::Forward, nz);
      double dg = recCell.basis(idir).norm ();
      double epsEff;
      if (useTensor)
        epsEff = recCell.basis(idir) ^ epsTensor ^ recCell.basis(idir)
               / (dg * dg);
      else
        epsEff = dielecConstant;
      SxVector<SxComplex16> vInG(nz), vInR(nz);
      vInG(0) = V0;
      if (idir == 0)  {
         cout << "V average: " <<  (vInG(0).re/defCell.volume * HA2EV) << " eV"
              << endl;
      }
      for (int z = 1; z < nz; ++z)  {
         double g = ((2 * z < nz) ?  z : (z - nz)) * dg;
         double g2 = g * g;
         vInG(z) = 0.;
         for (int iCharge = 0; iCharge < nCharges; ++iCharge)  {
            const SxCharge &it = chargeList(iCharge);
            vInG(z) += it.rhoRec(g2)
                       * exp ( - TWO_PI * I * g/dg * posRel(iCharge)(idir));
         }
         vInG(z) *= FOUR_PI / ( epsEff * g2 );
         if (mixedScreening)  {
            double epsEffEl = recCell.basis(idir) ^ eps2Tensor ^ recCell.basis(idir)
               / (dg * dg);
            for (const SxCharge &it : elChargeList)  {
               Coord rel = defCell.carToRel (it.pos);
               vInG(z) += it.rhoRec(g2)
                          * SxComplex16::phase ( -TWO_PI * g/dg * rel(idir))
                          * FOUR_PI / (epsEffEl * g2);
            }
         }
         if (q2tf > 1e-10)
            vInG(z) *= dielecConstant
               * SxPreconditioner::invCSRB(g2,q2tf,dielecConstant);
      }
      if (nz % 2 == 0) vInG(nz/2) = 0.;
      fft1d.fftForward (nz, vInG.elements, vInR.elements);
      vInR /= defCell.volume;
      vInR += vAlign;

      if (printRho)  {
         SxVector<SxComplex16> rhoInG(nz), rhoInR(nz);
         for (int z = 0; z < nz; ++z)  {
            double g = ((2 * z < nz) ?  z : (z - nz)) * dg;
            rhoInG(z) = 0.;
            for (int iCharge = 0; iCharge < nCharges; ++iCharge)  {
               const SxCharge &it = chargeList(iCharge);
               rhoInG(z) += it.rhoRec(sqr(g))
                          * exp (-TWO_PI * I * g/dg * posRel(iCharge)(idir));
            }
         }
         if (mixedScreening)  {
            cout << "WARNING: rhoModel.dat contains ionic charge only" << endl;
         }
         if (nz % 2 == 0) rhoInG(nz/2) = 0.;
         fft1d.fftForward (nz, rhoInG.elements, rhoInR.elements);
         rhoInR /= defCell.volume;
         rhoFile.writeXYPlot (dz, rhoInR.real ());
         rhoFile.print ("\n");
      }

      if (avgWidth > 1e-16)  {
         cout << "Averaging (" << (avgWidth / dz) << " points)" << endl;
         vInR = average (vInR, avgWidth / dz);
      }

      SxTextIO out("vline-eV-a"+SxString(idir)+".dat");
      SxVector<double> vLR = vInR.real () * HA2EV;
      // --- read potential files
      // at this point only vRef has to be read in
      if (printFullVLine)  {
         SxVector<double> vRef, vDef;
         vRef = readLine (refCell, refMesh, refPot, idir, nz, defCell,refPotFile);
         vDef = readLine (defCell, defMesh, defPot, idir, nz, defCell,defPotFile);
         if (avgWidth > 1e-16)  {
            vRef = average (vRef, avgWidth / dz);
            vDef = average (vDef, avgWidth / dz);
         }
         vDef *= HA2EV;
         vRef *= HA2EV;
         SxVector<double> deltaV = vDef - vRef;
         SxVector<double> vSR = deltaV - vLR;

         switch (outStyle)  {
            case Separate:
               out.writeXYPlot (dz, vLR, "");
               out.writeXYPlot (dz, deltaV, "");
               out.writeXYPlot (dz, vSR);
               if (printOrigPot)  {
                  out.print ("\n");
                  out.writeXYPlot (dz, vDef, "");
                  out.writeXYPlot (dz, vRef);
               }
               break;
            case OldStyle:
               out.writeXYPlot (dz, vLR, "&");
               out.writeMultiCol (6) << dz << deltaV << vSR;
               if (printOrigPot)  {
                  out.print("&\n");
                  out.writeMultiCol (6) << dz << vDef << vRef;
               }
               break;
            case AllColumn:
               if (printOrigPot)
                  out.writeMultiCol (6) << dz << vLR << deltaV << vSR
                                        << vDef << vRef;
               else
                  out.writeMultiCol (6) << dz << vLR << deltaV << vSR;
         }
      } else {
         out.writeXYPlot (dz, vLR);
      }

   }
   if (printRho) rhoFile.close ();

   cout << "vAlign=" << (vAlign *HA2EV) << " eV" << endl;

   // --- report energies
   cout << SX_SEPARATOR;
   cout << "=== Intermediate results (" << (useTensor ? "" : "un")
        << "screened) ===" << endl;
   cout << "Isolated energy       : " << eIso << endl;
   cout << "Periodic energy       : " << ePeriodic << endl;
   cout << "Difference (Hartree)  : " << ePeriodic - eIso << endl;
   cout << "Difference (eV)       : " << (ePeriodic - eIso) * HA2EV << endl;

   cout << SX_SEPARATOR;
   if (useTensor)  {
      cout << "Calculation performed with epsilon = " << epsTensor << endl
           << "Spherical harmonic average         = ";

   } else
      cout << "Calculation performed with epsilon = ";
   cout << dielecConstant << endl;
   if (mixedScreening)  {
      cout << "Electronic epsilon                 = " << eps2Tensor << endl
           << "Electronic isotropic average       = " << (1./eps2Inv) << endl;
   }

   cout << SX_SEPARATOR;
   cout << "Defect correction (eV): ";
   if (useTensor)
      cout << ((eIso - ePeriodic)  - totQ * vAlign) * HA2EV;
   else
      cout << ((eIso - ePeriodic) / dielecConstant - totQ * vAlign) * HA2EV;
   cout << " (incl. screening & alignment)" << endl;
}
