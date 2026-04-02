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
#include <SxTextIO.h>
#include <SxIO.h>
#include <SxException.h>
#include <cstdarg>

void SxTextIO::open (const char *name, Mode modeIn)
{
   if (isOpen ())  {
#ifndef NDEBUG
      cout << "Warning: auto-closing file '" << fileName << "'." << endl;
#endif
      fclose (fp);
   }
   mode = modeIn;
   fileName = name;
   const char *fmode = NULL;
   if (mode == WriteOnly || mode == WriteSx)
      fmode = "w";
   else if (mode == ReadOnly)
      fmode = "r";
   else if (mode == Append || mode == AppendSx)
      fmode = "a";
   SX_CHECK (fmode);
   fp = sxfopen (name, fmode);

   if (mode == WriteSx)  {
#ifndef NDEBUG
      cout << "Warning: hard-coded format 'matrix' when writing '" << name
           << "' in sx mode." << endl;
#endif
      fprintf (fp, "format matrix;\n");
   } else if (mode == AppendSx)  {
      mode = WriteSx;
   } else if (mode == Append)  {
      mode = WriteOnly;
   }
}

void SxTextIO::printf (const char *fmt, ...)
{
   SX_CHECK (fp);
   SX_CHECK (mode == WriteOnly || mode == WriteSx);
   if (!fmt) return;
   va_list arg;
   va_start(arg, fmt);
   vfprintf (fp, fmt, arg);
   va_end (arg);
}

void SxTextIO::writeXYPlot (const SxVecRef<double> &yVec)
{
   ssize_t n = yVec.getSize();
   SxVector<double> xVec (n);
   SxVector<double>::Iterator x = xVec.begin();
   for (ssize_t i=0; i<n; i++)  {
      *x = (Real8)i;
      x++;
   }
   writeXYPlot (xVec, yVec);
}


void SxTextIO::writeXYPlot (const SxVecRef<double> &xVec,
                        const SxVecRef<double> &yVec)
{
   SX_CHECK (xVec.getSize() == yVec.getSize(),
             xVec.getSize(),   yVec.getSize());
   SX_CHECK (isOpen ());
   SX_CHECK (mode == WriteOnly, mode);
   SxVecRef<double>::Iterator x, y;
   const char* format;
   if (yVec.abs ().maxval () < 1e-6)
      format="%15.12f\t%15.12g\n";
   else
      format="%15.12f\t%15.12f\n";

   for (x = xVec.begin(), y=yVec.begin(); x != xVec.end(); x++, y++)
      fprintf (fp, format, *x, *y);
}

SxVector<double> SxTextIO::readXYPlot ()
{
   SX_CHECK (isOpen ());
   SX_CHECK (mode == ReadOnly, mode);
   SxList<double> values;
   double x,y;
   char c;
   while (!feof(fp)) {
      while (fscanf(fp," #%c",&c) == 1) {
         while (fgetc(fp) != '\n' && !feof(fp));
      }
      int nRead = fscanf(fp, "%lf %lf",&x,&y);
         if (nRead == 2) {
            values << x << y;
         } else break;
   }
   ssize_t dim = values.getSize();
   if (dim == 0 ) {
      sxprintf ("noData read.\n");
      SX_EXIT;
   }
   if (dim % 2 == 1) {
      sxprintf ("incomplete XYPlot.\n");
      SX_EXIT;
   }

   SxVector<double> result (dim/2,2);
   SxList<double>::ConstIterator it = values.begin ();
   for (int i = 0; i < result.getNCols (); i++)  {
      result(i,0) = *it++;
      result(i,1) = *it++;
   }
   return result;

}

void SxTextIO::writeXYPlot (const SxVecRef<double> &xVec,
                           const SxVecRef<double> &yVec,
                           int   xDecimals,
                           int   yDecimals)
{
   SX_CHECK (xVec.getSize() == yVec.getSize(),
             xVec.getSize(),   yVec.getSize());
   SX_CHECK (xDecimals >= 0 && xDecimals <= 30, xDecimals);
   if (yDecimals < 0) yDecimals = xDecimals;
   SX_CHECK (yDecimals >= 0 && yDecimals <= 30, yDecimals);
   SX_CHECK (isOpen ());
   SX_CHECK (mode == WriteOnly, mode);

   SxVecRef<double>::Iterator x, y;
   for (x = xVec.begin(), y=yVec.begin(); x != xVec.end(); ++x, ++y)  {
      fprintf (fp, "%.*f\t%.*f\n", xDecimals, *x, yDecimals, *y);
   }
}


void SxTextIO::writeXYPlot (const SxVecRef<SxComplex16> &yVec)
{
   ssize_t n = yVec.getSize();
   SxVector<double> xVec (n);
   SxVecRef<double>::Iterator x = xVec.begin();
   for (ssize_t i=0; i<n; i++)  {
      *x = (Real8)i;
      x++;
   }
   writeXYPlot (xVec, yVec);
}




void SxTextIO::writeXYPlot (const SxVecRef<double>      &xVec,
                           const SxVecRef<SxComplex16> &yVec)
{
   SX_CHECK (xVec.getSize() == yVec.getSize(),
             xVec.getSize(),   yVec.getSize());
   SX_CHECK (isOpen ());
   SX_CHECK (mode == WriteOnly, mode);
   fprintf (fp, "# x | y.re | y.im | y^2\n");
   SxVecRef<double>::Iterator      x;
   SxVecRef<SxComplex16>::Iterator y;
   for (x = xVec.begin(), y=yVec.begin(); x != xVec.end(); x++, y++)  {
      fprintf (fp, "%15.12f\t%15.12f\t%15.12f\t%15.12f\n",
               *x, (*y).re, (*y).im, (*y).absSqr ());
   }
}

void SxTextIO::writeNXYPlot (const SxVecRef<double> &yMat)
{
   SX_CHECK (isOpen ());
   SX_CHECK (mode == WriteOnly, mode);
   ssize_t nRows = yMat.getNRows ();
   ssize_t nCols = yMat.getNCols ();
   for (ssize_t iRow = 0; iRow < nRows; iRow++)  {
      fprintf (fp, "%i ", (int)iRow + 1);
      for (ssize_t iCol = 0; iCol < nCols; iCol++)  {
         fprintf (fp, "%.6f ", yMat(iRow, iCol));
      }
      fprintf (fp, "\n");
   }
}

//----------------------------------------------------------------------------
// Surface plot utility
//----------------------------------------------------------------------------
void SxTextIO::writeXYZPlot (const SxVecRef<PrecCoeffG> &M)
{
   SX_CHECK (isOpen ());
   SX_CHECK (mode == WriteOnly);
   ssize_t iG, iGp, nG = M.getNRows ();

   for (iG = 0; iG < nG; iG++)  {
      for (iGp= 0; iGp< nG; iGp++)   {
         fprintf (fp, "%i  %i  %g  %g \n",
                       (int)iG, (int)iGp, M(iG,iGp).re, M(iG,iGp).im );
      }
      fprintf (fp, "\n");
   }
}

void SxTextIO::writeVec (const SxString &varName,
                  const SxVecRef<SxComplex16> &val)
{
   if (mode == WriteOnly) {
#ifndef NDEBUG
      cout << "Warning: no varName='" << varName
           << "' needed when writing to ASCII file." << endl;
#endif
      writeMat (val);
      return;
   }
   SX_CHECK (mode == WriteSx, mode);
   ssize_t j;
   fprintf (fp, "\n");


   fprintf (fp,"%s = [\n", varName.ascii ());
      SxString line;
      line = SxString("");
      for (j = 0; j < val.getSize (); j++) {
         line += SxString (val(j).re, SxString("(%8.7f, "));
         line += SxString (val(j).im,SxString("%8.7f)"));
         if (j < val.getSize() - 1) line += SxString(",");
         line += SxString (val(j).im,SxString("\n"));
      }
   fprintf ( fp,"%s", line.ascii ());
   fprintf (fp, "];\n");
}

SxVector<double>
SxTextIO::readMat (int nElemRow, int nElemCol)
{
   SX_CHECK (mode == ReadOnly);
   SX_CHECK (nElemRow > 0, nElemRow);
   SX_CHECK (nElemCol > 0, nElemCol);
   SxVector<double> val(nElemRow, nElemCol);
   const int BUFLEN = 1024000;
   ssize_t i, j;
   char buffer[BUFLEN];
   char *bs = NULL;
   SxList<SxString> list;
   SxString tk;

   for (i = 0; i < nElemRow; i++) {
      bs = fgets (buffer, BUFLEN, fp);
      tk = SxString(bs);
      list = tk.tokenize (' ');
      if (list.last () == "\n") list.removeLast ();
      if (list.getSize () != nElemCol)  {
         cout << "Unexpected matrix size in " << fileName << endl;
         cout << "In row " << (i+1) << " of " << nElemRow << endl;
         cout << "Expected  " << nElemCol << " colums, but found "
              << list.getSize () << endl;
         SX_QUIT;
      }
      for (j = 0; j < nElemCol; j++)  {
         try {
            val(i, j)=list(j).toDouble ();
         } catch (SxException e)  {
            e.print ();
            SX_QUIT;
         }
      }
   }
   return val;
}

void SxTextIO::writeMat (const SxVecRef<double> &val)
{
   SX_CHECK (isOpen ());
   SX_CHECK (mode == WriteOnly, mode);
   for (ssize_t i = 0; i < val.getNRows (); i++) {
      for (ssize_t j = 0; j < val.getNCols (); j++)  {
         fprintf (fp, "%13.12e ", val(i, j));
      }
      fprintf (fp, "\n");
   }
}

void SxTextIO::writeMat (const SxString &varName, const SxVecRef<double> &val)
{
   if (mode == WriteOnly) {
#ifndef NDEBUG
      cout << "Warning: no varName='" << varName
           << "' needed when writing to ASCII file." << endl;
#endif
      writeMat (val);
      return;
   }
   SX_CHECK (isOpen ());
   SX_CHECK (mode == WriteSx, mode);
   ssize_t i, j;

   for (i = 0; i < val.getNCols (); i++) {
      fprintf (fp, "[%13.12e", val(i, 0));
      for (j = 1; j < val.getNCols (); j++) {
         fprintf (fp, ", %13.12e", val(i, j));
      }
      fprintf(fp, "]");
      if (i < (val.getNCols () - 1))
         fprintf(fp, ",");
      fprintf (fp, "\n");
   }
   fprintf (fp, "];\n");
}

SxList<SxList<SxVector3<double> > >
SxTextIO::loadStructureFHI98 () const
{
   const int BUFLEN = 1024000;
   int i, is, ia, ic;
   char buffer[BUFLEN];
   char *bs = NULL;
   int nSpecies, nAtoms;
   SxString tk;
   SxList<SxList<SxVector3<double> > >  returnValue;
   SxList<SxString> list;

   for (i = 0; i < 3; i++)  bs = fgets (buffer, 1024, fp);

   bs = fgets (buffer, 1024, fp);
   tk = SxString(bs);
   nSpecies = tk.toInt ();

   returnValue.resize (nSpecies);

   for (is=0; is < nSpecies; is++) {
      bs = fgets (buffer, 1024, fp);
      tk = SxString(bs);
      nAtoms = tk.toInt ();
      returnValue(is).resize (nAtoms);

      bs = fgets (buffer, 1024, fp);

      for (ia=0; ia < nAtoms; ia++)  {
         bs = fgets(buffer, 1024, fp);
         tk = SxString(bs);
         list = tk.tokenize (' ');
         for (ic = 0; ic < 3; ic++) {
            returnValue(is)(ia)(ic)=list(ic).toDouble ();
         }
      }
   }
   return returnValue;
}


SxList<SxString> SxTextIO::loadSpeciesFHI98 () const
{
   cout << "FHI98 support should be replaced by standard formats\n"
              "Use add-on converters instead" << endl;
   SX_EXIT;
   const int BUFLEN = 1024000;
   int i, is, ia;
   char buffer[BUFLEN];
   char *bs = NULL;
   int nSpecies, nAtoms;
   SxString tk;
   SxList<SxString>  returnValue;


   for (i = 0; i < 3; i++)  bs = fgets (buffer, 1024, fp);

   bs = fgets (buffer, 1024, fp);
   tk = SxString(bs);
   nSpecies = tk.toInt ();

   for (is=0; is < nSpecies; is++) {
      bs = fgets (buffer, 1024, fp);
      tk = SxString(bs);

      nAtoms = tk.toInt ();

      bs = fgets (buffer, 1024, fp);
      SxString help (bs);
      SxList<SxString> helpList = help.tokenize('\n');
      returnValue.append (helpList(0));

      for (ia=0; ia < nAtoms; ia++)  {
         bs = fgets(buffer, 1024, fp);
      }
   }
   return returnValue;
}

SxMatrix3<double> SxTextIO::loadCellFHI98 () const
{
   cout << "FHI98 support should be replaced by standard formats\n"
              "Use add-on converters instead" << endl;
   SX_EXIT;
   const int BUFLEN = 1024000;
   int i, ic;
   char buffer[BUFLEN];
   char *bs = NULL;
   SxString tk;
   SxMatrix3<double>   returnValue;


   for (i = 0; i < 3; i++)  {
         bs = fgets(buffer, 1024, fp);
         tk = SxString(bs);
         SxList<SxString> list;
         list = tk.tokenize (' ');
         for (ic = 0; ic < 3; ic++) {
            try {
               returnValue(i, ic)=list(ic).toDouble ();
            } catch (SxException e)  {
               e.print ();
               SX_QUIT;
            }
         }
   }
   return returnValue;
}
SxTextIO::ColumnWrite::ColumnWrite (SxTextIO &parent, int prec,
                                    bool smartFormat, int width)
   : file(parent), dx(0.)
{
   if (prec > 0)  {
      SX_CHECK (prec < 20);
      format = SxString ("%");
      if (width > 0)  {
         SX_CHECK (width < 20, width);
         format += width;
      }
      format += ".";
      format += prec;
      format += smartFormat ? "g" : "f";
   }
}

SxTextIO::ColumnWrite::~ColumnWrite ()
{
   SX_CHECK (file.mode == WriteOnly, file.mode);
   FILE *fp = file.getFp ();
   // determine format
   for (auto &vecIt : data)  {
      if (format.getSize () > 0) break;
      for (double val : vecIt)  {
         if (fabs(val) > 1e-6)  {
            format="%15.12f";
            break;
         }
      }
   }
   if (format.getSize () == 0) format="%15.12g";
   const char* fmt = format.ascii ();
   const char* dxfmt = (dx < 1e-4) ? "%.6g\t" : "%.6f\t";

   // write out data
   ssize_t N = data(0).getNRows ();
   for (ssize_t i = 0; i < N; ++i)  {
      if (dx > 0.)
         fprintf (fp, dxfmt, dx * double(i));
      bool tab = false;
      for (auto &vecIt : data)  {
         for (ssize_t j = 0; j < vecIt.getNCols (); ++j)  {
            if (tab) fprintf(fp, "\t");
            tab = true;
            fprintf (fp, fmt, vecIt(i, j));
         }
      }
      fprintf (fp, "\n");
   }
}


