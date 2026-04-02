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

#ifndef _SX_TEXT_IO_H_
#define _SX_TEXT_IO_H_

#include <SxIO2.h>
#include <stdio.h>
#include <SxString.h>
#include <SxVector.h>
//#include <SxMatrix3.h>
//#include <SxTypes.h>

/** \brief Read and write ASCII files

    \author Christoph Freysoldt, freysoldt@mpie.de
 */
class SX_EXPORT_IO2 SxTextIO
{
   protected:
      /// File pointer
      FILE *fp;

   public:
      /// Name of the file
      SxString fileName;

      enum Mode {
         /// Read ASCII file
         ReadOnly,
         /// Write ASCII file
         WriteOnly,
         /// Append to ASCII file
         Append,
         /// Write sx-style file
         WriteSx,
         /// Append to sx-style file
         AppendSx,
         /// Not set
         Unknown
      } mode;

      /// Empty constructor
      SxTextIO () : fp(NULL), mode (Unknown) {}

      /// Constructor that directly opens a file
      SxTextIO (const SxString &filename, Mode modeIn = WriteOnly)
         : fp(NULL)
      {
         open (filename.ascii (), modeIn);
      }

      SxTextIO (const char *filename, Mode modeIn = WriteOnly)
         : fp(NULL)
      {
         open (filename, modeIn);
      }

      /// Destructor (closes file, if open)
      ~SxTextIO () {
         if (fp) fclose (fp);
      }

      /// Open a file
      void open (const SxString &filename, Mode modeIn = WriteOnly)
      {
         open (filename.ascii (), modeIn);
      }

      /// Open a file
      void open (const char *name, Mode modeIn);

      /// Close the file
      void close ()
      {
         if (fp) fclose (fp);
         fp = NULL;
         mode = Unknown;
      }

      /// Check if file is open
      bool isOpen () const  { return fp; }

      /// Get fp for direct read/write
      FILE* getFp ()
      {
         SX_CHECK (isOpen ());
         return fp;
      }

      /// Direct printing
      inline void print (const char* text)
      {
         SX_CHECK (isOpen ());
         SX_CHECK (mode == WriteOnly || mode == WriteSx);
         fprintf (fp, "%s", text);
      }
      /// Direct printing
      inline void print (const SxString &text)
      {
         print (text.ascii ());
      }
      /// Direct formated printing
      void printf (const char *fmt, ...) __SXCHECK_FORMAT(2,3);
                                        // implicit this is arg 1

      void writeXYPlot (const SxVecRef<double> &yVec);
      /// write i * dx, yVec(i)
      template<StrideType Layout>
      void writeXYPlot (double dx, const SxVecRef<double, Layout> &yVec,
                        const char* printFinal = NULL);
      void writeXYPlot (const SxVecRef<SxComplex16> &yVec);
      void writeXYPlot (const SxVecRef<double> &xVec,
                        const SxVecRef<double> &yVec);
      void writeXYPlot (const SxVecRef<double>    &xVec,
                        const SxVecRef<SxComplex16> &yVec);
      SxVector<double> readXYPlot ();

      /**
        \brief XY-plot with specified precision
        \param xVec x-values
        \param yVec y-values
        \param xDecimals number of decimals for x values
        \param yDecimals number of decimals for y values
                         If it is not given, or negative, xDecimals is used

        \note The number of x-values and y-values must be the same.
        */
      void writeXYPlot (const SxVecRef<double> &xVec,
                        const SxVecRef<double> &yVec,
                        int   xDecimals,
                        int   yDecimals = -1);

      void writeNXYPlot (const SxVecRef<double> &yMat);

   protected:
      /** \brief Class to write out arbitrary number of columns

          This class collects all the columns and writes it in
          the destructor.
        */
      class ColumnWrite {
         private:
            SxList<SxVecRef<double, GeneralMatrix> > data;
            SxString format;
            friend class SxTextIO;
            ColumnWrite (SxTextIO &parent, int prec, bool smartFormat,
                         int width);
            SxTextIO &file;
            double dx;

            ColumnWrite (const ColumnWrite&) = delete;
            ColumnWrite (ColumnWrite &&) = default;
         public:
            template<StrideType Layout>
            ColumnWrite& operator<< (const SxVecRef<double, Layout> &in)
            {
#ifndef NDEBUG
               if (data.getSize () > 0)  {
                  SX_CHECK (data(0).getNRows () == in.getNRows (),
                            data(0).getNRows (), in.getNRows ());
               }
#endif
               data.append (const_cast<SxVecRef<double,Layout>&>(in));
               return *this;
            }
            ColumnWrite& operator<< (double dx_)  {
               SX_CHECK (dx_ > 0.);
               dx = dx_;
               return *this;
            }
            ~ColumnWrite ();
      };

   public:
      /** \brief Write arbitrary number of columns

          \note Stream-like use:
          file.writeMultiCol () << col1 << col2 << col3 << col4;

          The first column can be replaced by <dx * i> by appending the
          dx value

       */
      ColumnWrite writeMultiCol (int prec = -1, bool smartFormat = false,
                                 int width = -1)
      {
         return ColumnWrite (*this, prec, smartFormat, width);
      }


      /** \author Abdullah Alsharif
        */
      void writeXYZPlot (const SxVecRef<PrecCoeffG> &);

      /*\brief  I/O routines for fhi98 structure format, needed for
                compatiblity purposes

        \deprecate
       */
      SxList<SxList<SxVector3<double> > >  loadStructureFHI98 () const;

      /** \deprecate */
      SxList<SxString>  loadSpeciesFHI98 () const;

      /** \deprecate */
      SxMatrix3<double> loadCellFHI98 () const;

      /**
        @brief Writes vector to sx file.
        @param varName  name of variable
        @param val      vector to be written.
        */
      void writeVec (const SxString &varName,
                     const SxVecRef<SxComplex16> &val);
      /**
        @brief Writes matrix to sx file.
        @param varName  name of variable
        @param val      matrix to be written.
        */
      void writeMat (const SxString &varName, const SxVecRef<double> &val);

      /**
        @brief Writes matrix to ASCII file.
        @param val  matrix to be written.
        */
      void writeMat (const SxVecRef<double> &val);

      /** @brief Reads matrix from ASCII file
        @param nElemRow number of rows
        @param nElemCol number of columns
        */
      SxVector<double> readMat (int nRow, int nCol);

};

template<StrideType Layout>
void SxTextIO::writeXYPlot (double dx, const SxVecRef<double, Layout> &yVec,
                            const char* printFinal)
{
   SX_CHECK (isOpen ());
   SX_CHECK (mode == WriteOnly, mode);
   const char* format = "%15.12f\t%15.12g\n";
   for (auto it : yVec)
      if (fabs(it) >= 1e-6)  { format="%.12f\t% 15.12f\n"; break; }

   SxVecIt<double, Layout> y;
   int ix = 0;
   for (y = yVec.begin(); y != yVec.end (); ix++, y++)
      fprintf (fp, format, dx * ix, *y);
   if (printFinal)
      fprintf(fp, "%s\n", printFinal);
}

#endif /* _SX_TEXT_IO_H_ */
