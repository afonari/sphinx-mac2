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

#include <SxUtil.h>
#include <SxCLI.h>
#include <SxCliStandard.h>
#include <SxAtomicStructure.h>
#include <SxElemDB.h>
#include <SxSpeciesData.h>
#include <SxRho.h>
#define a0 0.5291772083

int main (int argc, char **argv)
{
   // --- parse command line options
   SxCLI cli (argc, argv);
   
   // --- input and output files
   // standard sxb option
   bool sxbFile = sxbOption (cli);
   // standard input file option
   SxString inFile = inputOption (cli, sxbFile);

   SxString outFile = cli.option("-o","file","XSF output file to be written")
                      .toString("input.xsf");
   SxString rhoFile = cli.option ("-r|--rho","file","rho file")
                        .toString("");

   cli.finalize ();
   
   // --- read from input file
   initSPHInXMath ();

   SxAtomicStructure structure = structureFromInput (inFile, sxbFile);
   const SxArray<SxString> &chemName = structure.getElements ();
   SxElemDB elemDB;

   cout << inFile << " -> " << outFile << endl; 

   // --- write down structure
   ofstream file;
   file.open(outFile.ascii ());
   file << "CRYSTAL" << endl;
   file << "PRIMVEC" << endl;
   for(int i = 0; i < 3; i++)   {
      file << fixed << structure.cell(0,i) * a0 
           << " "   << structure.cell(1,i) * a0
           << " "   << structure.cell(2,i) * a0
           << endl;
   }
   
   file << "CONVVEC" << endl;
   for(int i = 0; i < 3; i++)   {
      file << fixed << structure.cell(0,i) * a0
           << " "   << structure.cell(1,i) * a0
           << " "   << structure.cell(2,i) * a0
           << endl;
   }

   file << "PRIMCOORD" << endl;
   int nTlAtoms = structure.nTlAtoms;
   file << nTlAtoms << " " << 1 << endl;
   for (int i = 0; i < nTlAtoms; i++)   {
      int iSpecies = structure.getISpecies(i);
      int idx = elemDB.getIdx(chemName(iSpecies));
      file << fixed << idx 
           << " " << structure(i)(0) * a0 
           << " " << structure(i)(1) * a0
           << " " << structure(i)(2) * a0
           << endl;
   }
   file << endl;
   
   // add rho
   if (rhoFile.getSize () > 0)   {
      SxMatrix3<double> cell;
      SxVector3<int>    dim;
      SxBinIO io;
      io.open (rhoFile, SxBinIO::BINARY_READ_ONLY);
      RhoR rhoIn = io.readMesh (&cell, &dim);
      file << "BEGIN_BLOCK_DATAGRID_3D" << endl;
      file << rhoFile.ascii () << endl;
      file << "BEGIN_DATAGRID_3D_" << rhoFile.ascii() << endl;
      file << dim(0) << " " << dim(1) << " " << dim(2) << endl;
      file << fixed << 0. << " " << 0. << " " << 0. << endl;
      for(int i = 0; i < 3; i++)   {
         file << fixed << structure.cell(0,i) * a0
            << " "   << structure.cell(1,i) * a0
            << " "   << structure.cell(2,i) * a0
            << endl;
      }
      file << endl;
      int meshSize = dim(0)*dim(1)*dim(2);
      for (int i = 0; i < meshSize; i++)   {
         file << scientific << rhoIn(0)(i) << " ";
         if ((i+1) % 8 == 0) file << endl;
      }
      file << endl;
      file << "END_DATAGRID_3D" << endl;
      file << "END_BLOCK_DATAGRID_3D" << endl;
   }

   file.close();

   return 0;
}

