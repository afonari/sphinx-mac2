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

#include <SxLocpot.h>
#include <SxConstants.h>
#include <SxString.h>
#include <SxSort.h>
#include <SxFileParser.h>

SxLocpot::SxLocpot ()
{
   // empty
}

SxLocpot::SxLocpot (const SxString &file_)
{
   setFilename (file_);
}

SxLocpot::~SxLocpot ()
{
   // empty
}

// Specifies the Locpot filename.
void SxLocpot::setFilename (const SxString &file_)
{
   filename = file_;
}

void SxLocpot::read ()
{

   SxArray<SxString> chemNames;

   SxFileParser fp (filename);
   fp.nextLine ();
   fp.topic ("scale");
   double scale = fp.getDouble ();
   fp.topic ("cell");
   for (int iBasis = 0; iBasis < 3; ++iBasis)
      for (int xyz = 0; xyz < 3; ++xyz)
         fp >> cell(xyz,iBasis);
   cell*= scale * A2B;
   cell.setup ();
   fp.nextLine ();
   fp.topic ("chemical species info");
   // read first line
   SxString line = fp.getLine ();
   // --- optional vasp comment line for species information chem names
   SxList<SxString> nAt = line.left('\n').simplifyWhiteSpace ().tokenize(' ');
   if (nAt.getSize () == 0)  {
      cout << "Problem while reading VASP file " << filename << endl;
      cout << "Syntax error: cannot read structure information" << endl;
      cout << "Please check file format or report this error." << endl;
      SX_QUIT;
   }
   bool names = false;
   nAt(0).toInt (&names); // try to convert; otherwise set names=true
   // --- if chem names read next line, these are the atom numbers
   if (names) {
      SxList<SxString>::ConstIterator it;
      for (it = nAt.begin(); it != nAt.end(); ++it)
         chemNames.append(*it);
      line = fp.getLine ();
      nAt = line.left('\n').simplifyWhiteSpace ().tokenize(' ');
   } else  {
      int number = 0;
      SxList<SxString>::ConstIterator it;
      for (it = nAt.begin(); it != nAt.end(); ++it, number++)
         chemNames.append(SxString("POTCAR-ELEMENT-") + number);
   }
   if (nAt.getSize () == 0)  {
      cout << "Problem while reading VASP file " << filename << endl;
      cout << "Syntax error: number of atoms line is empty" << endl;
      cout << "Please check file format or report this error." << endl;
      SX_QUIT;
   }

   // --- number of atoms
   SxList<int> nAtomList;
   SxVector<int> nAtoms(chemNames.getSize ());
   for (SxList<SxString>::Iterator it = nAt.begin ();
         it != nAt.end (); ++it)
   {
      try {
         if ((*it).getSize () > 0)
            nAtomList << it->toInt ();
      } catch (const SxException &e)  {
         e.print ();
         SX_QUIT;
      }
   }
   nAtoms = nAtomList;
   // one line (type of coordinates)
   fp.topic ("type of coordinates");
   line = fp.getItem ();
   char letter = line.toUpper()(0);
   bool direct = true;
   switch (letter)  {
   case 'D' : direct = true;
              break;
   case 'C' : direct = false;
              break;
   case 'K' : direct = false;
              break;
   default  : // unknown
              fp.where ();
              cout << ": Failed to understand coordinate type: '"
                   << line << "'" << endl;
              SX_QUIT;
   }
   fp.nextLine ();

   structure.cell = cell;
   structure.startCreation ();
   Coord coord;
   int iAtom, iSpecies;
   fp.topic ("atomic coordinates");
   for (iSpecies = 0; iSpecies < nAtoms.getSize(); ++iSpecies)  {
      structure.newSpecies ();
      for (iAtom = 0; iAtom < nAtoms(iSpecies); ++iAtom)  {
         // --- atom coordinates
         fp >> coord(0) >> coord(1) >> coord(2);
         fp.nextLine ();

         if (direct) coord = structure.cell ^ coord; // Angstroem -> Bohr is already in the cell
         else coord = A2B * coord; // Angstroem -> Bohr
         structure.addAtom (iSpecies, coord);
      }
   }
   structure.endCreation ();
   structure.atomInfo->meta.update (SxAtomicStructure::Elements, chemNames);

   // --- mesh size
   fp.topic ("mesh");
   fp >> mesh(0) >> mesh(1) >> mesh(2);
   cout << "VASP mesh: " << mesh << endl;
   // --- read potential
   fp.topic ("potential");
   potential.resize (mesh.product ());
   SxVector3<int> x;
   for (x(2) = 0; x(2) < mesh(2); x(2)++)
      for (x(1) = 0; x(1) < mesh(1); x(1)++)
         for (x(0) = 0; x(0) < mesh(0); x(0)++)
            fp >> potential(mesh.getMeshIdx(x, SxMesh3D::Positive));
}
