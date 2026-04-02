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

#include <SxHessian.h>
#include <SxNeighbors.h>

SxHessian::SxHessian (const SxVecRef<double> &full)
{
   int nAtoms = (int)full.getNCols () / 3;
   SX_CHECK (full.getNCols () == 3 * nAtoms, full.getNCols (), nAtoms);
   SX_CHECK (full.getNRows () == 3 * nAtoms, full.getNRows (), nAtoms);

   hessian.resize (nAtoms);
   neighbors.resize (nAtoms);
   for (int ia = 0; ia < nAtoms; ++ia)  {
      hessian(ia).resize (nAtoms);
      neighbors.resize (nAtoms);
      for (int ja = 0; ja < nAtoms; ++ja)  {
         neighbors(ia)(ja) = ja;
         SxMatrix3<double> &D = hessian(ia)(ja);
         for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
               D(i,j) = full(3 * ia + i, 3*ja + j);
      }
   }
}

SxVector<double> SxHessian::getFull () const
{
   int nAtoms = (int)hessian.getSize ();
   SxVector<double> res(3 * nAtoms, 3*nAtoms);
   res.set (0.);

   for (int ia = 0; ia < nAtoms; ++ia)  {
      for (int in = 0; in < nAtoms; ++in)  {
         int ja = neighbors(ia)(in);
         const SxMatrix3<double> &D = hessian(ia)(ja);
         for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
               res(3 * ia + i, 3*ja + j) = D(i,j);
      }
   }
   return res;
}

void SxHessian::writeFull (SxBinIO &io,
                           const SxVecRef<double> &full)
{
   SX_CHECK (full.getNCols () == full.getNRows (), full.getNCols (), full.getNRows ());
   int nDof = (int)full.getNCols ();
   SX_CHECK (nDof > 0, nDof);
   try  {
      io.addDimension ("nDoF", nDof);
      io.write ("hessian", full, "nDoF", "nDoF");
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

SxVector<double> SxHessian::readFull (const SxBinIO &io)
{
   int nDoF = io.getDimension ("nDoF");
   SxVector<double> res (nDoF, nDoF);
   try {
      io.read ("hessian", &res, nDoF, nDoF);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   return res;
}

void SxHessian::write (const SxString &fileName,
                       const SxVecRef<double> &hessianFull,
                       const SxAtomicStructure &str)
{

   SX_CHECK (hessianFull.getNCols () == hessianFull.getNRows (),
             hessianFull.getNCols (), hessianFull.getNRows ());
   SX_CHECK (3 * str.getNAtoms () == hessianFull.getNCols (),
             str.getNAtoms (), hessianFull.getNCols ());
   try {
      SxBinIO io(fileName, SxBinIO::BINARY_WRITE_ONLY);
      writeFull (io, hessianFull);
      str.write (io);
      io.setMode (SxBinIO::WRITE_DATA);
      writeFull (io, hessianFull);
      str.write (io);
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

