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

#include <SxIonicCoreFit.h>
#include <SxNeighbors.h>
#include <SxEESGProj.h>

void SxIonicCoreFit::set (int nPts, double rMaxIn, int order)
{
   rMax = rMaxIn;
   splineOrder = order;
   nPoints = nPts + splineOrder;
   dr = rMax / double(nPts);
   fitA.resize (0);
   fitb.resize (0);
   param.resize (0);
}

SxVector<double>
SxIonicCoreFit::getParamDeriv ()
{
   ssize_t nr = mesh.getSize ();
   SxVector<double> res(nPoints * structure.getNSpecies (), nr);
   int mode = SxNeighbors::StoreIdx | SxNeighbors::StoreDistSqr
            | SxNeighbors::IncludeZeroDistance;
   double invDr = 1. / dr;
#ifdef USE_OPENMP
#  pragma omp parallel
#endif
   {
      SxNeighbors nn;
      SxVector<double> bspline(splineOrder + 1);
#ifdef USE_OPENMP
#     pragma omp for
#endif
      for (ssize_t ir = 0; ir < nr; ++ir)  {
         Coord r = meshCell.relToCar (mesh.getMeshVec (ir, SxMesh3D::Positive));
         nn.compute (grid, structure, r, rMax, mode);
         res.colRef (ir).set (0.);
         for (ssize_t in = 0; in < nn.idx.getSize (); ++in) {
            int is = chemId(structure.getISpecies (nn.idx(in)));
            double rad = sqrt(nn.distSqr(in));
            double xn = rad * invDr;
            int jN = int(xn);
            double x = xn - double(jN);

            SxEESGProj::cardBSplineAll (x, splineOrder + 1, bspline);
            for (int j = max(0, jN - nPoints + splineOrder + 1); j <= splineOrder; ++j)
               res(is * nPoints + jN + splineOrder - j, ir) += bspline(j);
         }
      }
   }
   return res.transpose ();
}

void SxIonicCoreFit::set (const SxAtomicStructure &structureIn,
                          const SxMesh3D &meshIn)
{
   structure = structureIn;
   grid = SxGrid (structure, 3);
   mesh = meshIn;
   meshCell = SxCell(structure.cell.basis(0) / mesh (0),
                     structure.cell.basis(1) / mesh (1),
                     structure.cell.basis(2) / mesh (2));
   const SxArray<SxString> &elemNew = structureIn.getElements ();
   chemId.resize (structure.getNSpecies ());
   SX_LOOP(is)  {
      int js = (int)chemNames.findPos (elemNew(is));
      if (js < 0)  {
         js = (int)chemNames.getSize ();
         chemNames.append (elemNew(is));
      }
      chemId(is) = js;
   }
   ssize_t nParam = chemNames.getSize () * nPoints;
   if (fitA.getSize () > 0)  {
      SxVector<double> A(nParam, nParam);
      A.set (0.);
      A.getRef<SubMatrix> (0, fitA.getNRows (), fitA.getNCols ()) = fitA;
      fitA = std::move (A);
      fitb.resize (nParam, true);
   } else {
      fitA.reformat (nParam, nParam);
      fitA.set (0.);
      fitb.resize (nParam);
      fitb.set (0.);
   }
}

void SxIonicCoreFit::addFitData (const SxVector<double> &data,
                                 const SxMesh3D &meshIn,
                                 const SxAtomicStructure &structureIn,
                                 double weight)
{
   SX_CHECK (meshIn.getSize () == data.getSize (),
             meshIn.getSize (), data.getSize ());
   set (structureIn, meshIn);
   SxVector<double> fitP = getParamDeriv ();

   fitb.plus_assign_ax (weight, fitP.overlap (data));
   fitA.plus_assign_ax (weight, fitP.overlap (fitP));
}

void SxIonicCoreFit::computeFit (bool zeroAvg)
{
   for (int is = 0; is < chemNames.getSize (); ++is)  {
      ssize_t off = is * nPoints;
      for (int i = 0, j = splineOrder - 1; i < j; ++i, --j)  {
         // symmetrize: C(-i) = C(i-p-1) => Cn(i) = Cn(p-1-i)
         for (ssize_t k = 0; k < fitA.getNRows (); ++k)  {
            fitA(k, off+i) += fitA(k, off+j);
            fitA(k, off+j) = 0.;
         }
         for (ssize_t k = 0; k < fitA.getNCols (); ++k)  {
            fitA(off+i, k) += fitA(off+j, k);
            fitA(off+j, k) = 0.;
         }
         fitb(off+i) += fitb(off+j);
         fitb(off+j) = 0.;
         fitA(off+j, off+i) = -1.;
         fitA(off+j, off+j) = 1.;
      }
   }
   param = fitA.solve (fitb);
   param.reshape (nPoints, chemNames.getSize ());
   // align fit shape to zero at rMax
   if (zeroAvg)  {
      SxVector<double> C1(nPoints);
      C1.set (1.);
      for (int is = 0; is < chemNames.getSize (); ++is)  {
         double avg = 0., one = 0.;
         const SxVecRef<double> Ci = param.getRef (nPoints, nPoints * is);
         for (double r = 0.; r < dr * nPoints; r += 0.01 * dr)
         {
            double r2 = r * r;
            avg += r2 * getVal (r, Ci);
            one += r2 * getVal (r, C1); // = 1 for r <= rMax
         }
         param.getRef (nPoints, nPoints * is)-= avg/one;

      }
   }
}

SxVector<double>
SxIonicCoreFit::getFunc () const
{
   SX_CHECK (param.getSize () == nPoints * chemNames.getSize (),
             param.getSize (), nPoints, chemNames.getSize ());
   ssize_t nr = mesh.getSize ();
   SxVector<double> res(nr);
   int mode = SxNeighbors::StoreIdx | SxNeighbors::StoreDistSqr
            | SxNeighbors::IncludeZeroDistance;
#ifdef USE_OPENMP
#  pragma omp parallel
#endif
   {
      SxVector<double> bspline(splineOrder + 1);
      SxNeighbors nn;
#ifdef USE_OPENMP
#     pragma omp for
#endif
      for (ssize_t ir = 0; ir < nr; ++ir)  {
         Coord r = meshCell.relToCar (mesh.getMeshVec (ir, SxMesh3D::Positive));
         nn.compute (grid, structure, r, rMax + splineOrder*dr, mode);
         res(ir) = 0.;
         for (ssize_t in = 0; in < nn.idx.getSize (); ++in) {
            int is = chemId(structure.getISpecies (nn.idx(in)));
            double rad = sqrt(nn.distSqr(in));
            res(ir) += getVal (rad, param.getRef (nPoints, nPoints * is), bspline);
         }
      }
   }
   return res;
}

double SxIonicCoreFit::getVal (double x, const SxVecRef<double> &Cn,
                               SxVector<double> &work) const
{
   double sigmoid = 1.;
   if (x > 0.8 * rMax)  {
      // --- smooth truncation from 80%-100% of rMax
      //     via Cardinal B spline integral
      if (x < rMax) {
         double xx = (x - 0.8 * rMax) / (0.2 * rMax);
         xx *= splineOrder + 1;
         int ix = int(xx);
         xx -= ix;
         SxEESGProj::cardBSplineAll (xx, splineOrder + 1, work);
         for ( ; ix >= 0; --ix)
            sigmoid -= work(ix);
      } else {
         return 0.;
      }
   }
   x /= dr;
   int jN = int(x);
   if (jN >= nPoints) return 0.;
   x -= jN;
   double res = 0.;
   SxEESGProj::cardBSplineAll (x, splineOrder + 1, work);
   for (int j = max(0, jN - nPoints + splineOrder + 1); j <= splineOrder; ++j)
      res += Cn(jN + splineOrder - j) * work(j);
   return res * sigmoid;
}

