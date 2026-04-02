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

#ifndef _SX_BIN_IO_H_
#define _SX_BIN_IO_H_

#include <SxIO2.h>
#include <stdio.h>
#include <SxString.h>
#include <SxException.h>
#include <SxVector.h>
//#include <SxTypes.h>
#include <SxMatrix3.h>
#ifdef _USE_H5
#  include <hdf5.h>
#endif /* _USE_H5 */

/**
  @ingroup Communication
  */
class SX_EXPORT_IO2 SxBinIO
{
   public:
      SxString filename;

      int ncId;
      FILE *fp;
      bool isOpen;
      mutable int err;
      /** The mode determines which kind of file is created, as well
          as the allowed operations.
        */
      enum Mode { /// Read netcdf binary file
                  BINARY_READ_ONLY,
                  /// Write netcdf binary file (classic format)
                  BINARY_WRITE_ONLY,
                  /// Write netcdf binary file (large file format, 64bit)
                  BINARY_WRITE_LARGE,
                  /// Write parallel NetCDF4 file using MPI-IO
                  BINARY_WRITE_PARALLEL,
                  //SX_WRITE_ONLY, SX_APPEND_ONLY,
                  /// Open as a scratch file (allow reading and writing)
                  SCRATCH_READ_WRITE,
                  /// Uninitialized
                  UNKNOWN };
      enum NcMode { WRITE_HEADER=0, WRITE_DATA=1 };
      enum Mode mode;
      mutable int  ncMode;   // 0 - write header, 1 - write data


      SxBinIO ();
      SxBinIO (const SxBinIO &);
      SxBinIO (const SxString &filename, Mode mode);
      ~SxBinIO ();
      void open (const SxString &filename, Mode mode);

      void close ();

      void ncError (int) const;
      void ncError () const;
      /// Raise error from reading/writing variable or dimension 'name'
      void ncError (int, const SxString &name) const;

      SxBinIO &operator= (const SxBinIO &);

      /** \brief Creates a unique file name for temprary files

          For creation of temprary files you have to ensure to use a
          filename, that is unique. This function returns such a
          filename.
          \param tmpDir  location of temprary directory
        */
      static SxString createTempName (const SxString &tmpDir="/tmp");

      int  getSize (const SxString &name) const;
      void setMode (int) const;
      void addDimension (const SxString &dimName, int size) const;
      int  getDimension (const SxString &dimName) const;
      bool contains (const SxString &varName) const;
      bool containsDim (const SxString &dimName) const;

      /// \name Scalar values
      //@{
      /// Write a string
      void write (const SxString &varName, const SxString &val) const;
      /// Write a double to a netCDF file
      void write (const SxString &varName, double val) const;
      /// Read a double from a netCDF file
      void read (const SxString &varName, double *val) const;
      //@}

      void read  (const SxString &varName, SxString *val) const;
      // --- SxList<int>
      void read  (const SxString &varName, SxList<int> *val, int nElem) const;

      /// \name SxVector3
      //@{
      /// write SxVector3<int>
      void write (const SxString &varName,
                  const SxVector3<int> &val,
                  const SxString &dimName) const;
      /// read SxVector3<int>
      void read  (const SxString &varName,
                  SxVector3<int> *val) const;
      /// write SxVector3<double>
      void write (const SxString &varName,
                  const SxVector3<double> &val,
                  const SxString &dimName) const;
      /// read SxVector3<double>
      void read  (const SxString &varName,
                  SxVector3<double> *val) const;
      //@}

      /// \name SxMatrix3
      //@{
      /// Write SxMatrix3 to netCDF
      void write (const SxString &varName,
                  const SxMatrix3<double> &val,
                  const SxString &dimName) const;
      /// Read SxMatrix3 from netCDF
      void read  (const SxString &varName,
                  SxMatrix3<double> *val) const;
      /** @param dim1Name dimension of array, i.e. val.getSize ()
          @param dim2Name dimension of matrix, i.e. 3
          @par Example:
        \code
   SxArray<SxMatrix3<double> > array;

   ...

   io.addDimension ("arraySize", array.getSize ());
   io.addDimension ("xyz", 3);
   io.write ("array", array, "arraySize", "xyz");
        \endcode
          @brief Writes an array of 3x3 matrices
        */
      void write (const SxString &varName,
                  const SxArray<SxMatrix3<double> > &val,
                  const SxString &dim1Name,
                  const SxString &dim2Name) const;
      /** The array must have the correct size.
          @par Example:
        \code
   int size = io.getDimension ("arraySize");
   SxArray<SxMatrix3<double> > array;
   array.resize (size);
   io.read ("array", &array);
        \endcode
          @brief Reads in an array of 3x3 matrices
        */
      void read  (const SxString &varName,
                  SxArray<SxMatrix3<double> > *val) const;
      //@}
      /// \name SxVector
      //@{
      /**
        @brief Writes vector to netCDF file.
        @param varName    netCDF variable name
        @param val        vector to be written
        @param dimName    netCDF name of the dimension
        */
      template<class T>
      inline void write (const SxString &varName,
                         const SxVecRef<T> &val,
                         const SxString &dimName) const
      {
         writeVec (varName, val, dimName);
      }

      /**
        @brief Reads vector from netCDF file
        @param varName  netCDF variable name
        @param val      pointer to a nElem vector.
                        It must have the correct size.
        @param nElem    number of elements in the vector
        */
      template<class T>
      inline void read (const SxString &varName,
                        SxVecRef<T> *val,
                        int nElem) const
      {
         readVec (varName, val, nElem);
      }

      void writeVec (const SxString &varName,
                  const SxVecRef<int> &val,
                  const SxString &dimName,
                  int offset=0, int nelem=0, int localOffset=0) const;
      void readVec (const SxString &varName,
                        SxVecRef<int> *val,
                        int nElem, int offset=0) const;


      void writeVec (const SxString &varName,
                  const SxVecRef<double> &val,
                  const SxString &dimName,
                  int offset=0, int nelem=0, int localOffset=0) const;
      void readVec  (const SxString &varName,
                     SxVecRef<double> *val,
                     int nElem, int offset=0) const;


      void writeVec (const SxString &varName,
                     const SxVecRef<SxComplex16> &val,
                     const SxString &dimName, int offset=0) const;
      void readVec  (const SxString &varName,
                     SxVecRef<SxComplex16> *val,
                     int nElem, int offset=0) const;

      // --------------------------------------------------------------------
      /**
        @param varName    netCDF variable name
        @param val        matrix to be written. The matrix may be a submatrix
                          of a larger netCDF matrix
        @param dimRowName row    (1st index) name of the netCDF matrix
        @param dimColName column (2nd index) name of the netCDF matrix
        @param offsetRow
        @param offsetCol  offsets in a larger matrix
        @brief Writes matrix to netCDF file.
        */
      template<class T>
      inline void write (const SxString &varName,
                         const SxVecRef<T> &val,
                         const SxString &dimRowName,
                         const SxString &dimColName) const
      {
         writeMat (varName, val, dimRowName, dimColName);
      }

      /**
        @param val pointer to a nElem1 x nElem2 matrix.
                   It must have the correct size.
        @param nElemRow number of rows
        @param nElemCol number of columns
        @brief Reads matrix from netCDF file
        */
      template<class T>
      inline void read (const SxString &varName,
                        SxVecRef<T> *val,
                        int nElemRow, int nElemCol) const
      {
         readMat (varName, val, nElemRow, nElemCol);
      }


      /**
        @param varName    netCDF variable name
        @param val        matrix to be written. The matrix may be a submatrix
                          of a larger netCDF matrix
        @param dimRowName row    (1st index) name of the netCDF matrix
        @param dimColName column (2nd index) name of the netCDF matrix
        @param offsetRow
        @param offsetCol  offsets in a larger matrix
        @brief Writes matrix to netCDF file.
        */
      void writeMat (const SxString &varName,
                     const SxVecRef<int> &val,
                     const SxString &dimRowName, const SxString &dimColName,
                     int offsetRow = 0, int offsetCol = 0) const;
      /**
        @param val pointer to a nElem1 x nElem2 matrix.
                   It must have the correct size.
        @param nElemRow number of rows
        @param nElemCol number of columns
        @param offsetRow
        @param offsetCol offsets in a larger matrix
        @brief Reads matrix from netCDF file
        */
      void readMat (const SxString &varName,
                 SxVecRef<int> *val,
                 int nElemRow, int nElemCol,
                 int offsetRow = 0, int offsetCol = 0) const;

      /**
        @param varName    netCDF variable name
        @param val        matrix to be written. The matrix may be a submatrix
                          of a larger netCDF matrix
        @param dim1Name   row    (1st index) name of the netCDF matrix
        @param dim2Name   column (2nd index) name of the netCDF matrix
        @param offsetRow
        @param offsetCol offsets in a larger matrix
        @brief Writes matrix to netCDF file.
        */
      void writeMat (const SxString &varName,
                  const SxVecRef<double> &val,
                  const SxString &dim1Name, const SxString &dim2Name,
                  int offsetRow = 0, int offsetCol = 0) const;
      /**
        @param varName     netCDF variable name
        @param val         matrix containing the rows to be written
        @param dim1Name    row    (1st index) name of the netCDF matrix
        @param dim2Name    column (2nd index) name of the netCDF matrix
        @param rowOffset   starting row
        @param rowNelem    number of rows to write
        @param localOffset starting row within val matrix
        @brief Writes some rows of a matrix in the netCDF file.
        */
      void writeRows (const SxString &varName,
                     const SxVecRef<double> &val,
                     const SxString &dim1Name, const SxString &dim2Name,
                     int rowOffset, int rowNelem, int localOffset) const;
      /**
        @param val pointer to a nElem1 x nElem2 matrix.
                   It must have the correct size.
        @param nElemRow number of rows
        @param nElemCol number of columns
        @param offsetRow
        @param offsetCol offsets in a larger matrix
        @brief Reads matrix from netCDF file
        */
      void readMat (const SxString &varName,
                    SxVecRef<double> *val,
                    int nElemRow, int nElemCol,
                    int offsetRow = 0, int offsetCol = 0) const;


      /**
        @param varName    netCDF variable name
        @param val        matrix to be written. The matrix may be a submatrix
                          of a larger netCDF matrix
        @param dim1Name   row    (1st index) name of the netCDF matrix
        @param dim2Name   column (2nd index) name of the netCDF matrix
        @param offsetRow
        @param offsetCol offsets in a larger matrix
        @brief Writes matrix to netCDF file.
        */
      void writeMat (const SxString &varName,
                  const SxVecRef<SxComplex16> &val,
                  const SxString &dim1Name, const SxString &dim2Name,
                  int offsetRow = 0, int offsetCol = 0) const;

      /**
        @param val pointer to a nElem1 x nElem2 matrix.
                   It must have the correct size.
        @param nElemRow number of rows
        @param nElemCol number of columns
        @param offsetRow
        @param offsetCol offsets in a larger matrix
        @brief Reads matrix from netCDF file
        */
      void readMat (const SxString &varName,
                 SxVecRef<SxComplex16> *val,
                 int nElemRow, int nElemCol,
                 int offsetRow = 0, int offsetCol = 0) const;
      //@}

      /// \name meshes
      //@{
      /// Write 3d mesh
      void writeMesh (const SxVecRef<double>  &mesh,
                      const SxMatrix3<double> &cell,
                      const SxVector3<int>    &dim);
      /// Write 3d meshes
      void writeMesh (const SxArray<SxVector<double> > &meshes,
                      const SxMatrix3<double>          &cell,
                      const SxVector3<int>             &dim) const;
      /// Read 3d meshes
      SxArray<SxVector<double> >  readMesh (SxMatrix3<double> *cellPtr=NULL,
                                            SxVector3<int>    *dimPtr =NULL)
                                            const;
      //@}


      /** \brief Add a double matrix variable
        \note This is meant for cases where the data is not available
              as a matrix, but will be written with the writeRow routine
        */
      void addDoubleVar (const SxString &varName,
                         const SxString &dim1Name,
                         const SxString &dim2Name);

      /** \brief Write a row in a double matrix variable
        @param varName name of the variable
        @param val     data to be written
        @param ir      row index
        @param offset  optional column offset

        \note This is meant for cases where the data is not available
              as a matrix.
        \note Note that netCDF handles rows more efficiently than columns
              since we use the C interface.
        */
      void writeRow (const SxString &varName,
                     const SxVecRef<double> &val,
                     int ir,
                     int offset = 0);

      /** \brief Read a row from a double matrix variable
        @param varName name of the variable
        @param val     data to be read
        @param ir      row index
        @param offset  optional column offset

        \note This is meant for cases where the data is not stored
              as a matrix.
        \note Note that netCDF handles rows more efficiently than columns
              since we use the C interface.
        */
      void readRow (const SxString &varName,
                    SxVecRef<double> *val,
                    int ir,
                    int offset = 0) const;

};
#endif // _SX_BIN_IO_H_

