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

#include <SxConfig.h>
#include <SxBinIO.h>
#include <stdio.h>
#include <stdlib.h>
#include <SxError.h>
#include <unistd.h>

#ifdef USE_PARALLEL_NETCDF4
#include <SxLoopMPI.h>
#include <mpi.h>
#include <netcdf_par.h>
#endif

#include <netcdf.h>
#undef  NC_MAX_VARS
#define NC_MAX_VARS 10
#undef  MAX_NC_VARS
#define MAX_NC_VARS 10

/** Check if mode is one of the netcdf modes */
#define IS_NETCDF_MODE(mode)  (   mode == BINARY_READ_ONLY \
                               || mode == BINARY_WRITE_ONLY \
                               || mode == BINARY_WRITE_LARGE \
                               || mode == BINARY_WRITE_PARALLEL)
/** Check if mode is one of the netcdf write modes */
#define IS_NETCDF_WRITE(mode) (   mode == BINARY_WRITE_ONLY \
                               || mode == BINARY_WRITE_LARGE \
                               || mode == BINARY_WRITE_PARALLEL)

#ifdef CYGWIN
   extern "C" int mkstemp (char *);
#endif 

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
SxBinIO::SxBinIO ()
   : ncId(-1), fp(NULL), isOpen(false), err(0), mode(UNKNOWN)
{
   ncMode = WRITE_HEADER;
}


SxBinIO::SxBinIO (const SxString &filename_, Mode mode_)
   : ncId(-1), fp(NULL), isOpen(false), err(0), mode(UNKNOWN)
{
   open (filename_, mode_);
}


SxBinIO::SxBinIO (const SxBinIO &in)
   : ncId(-1), fp(NULL), isOpen(false), err(0), mode(UNKNOWN)
{
   // Using this routine may be a mistake, in most cases due to
   // foo (SxBinIO io) instead of foo(SxBinIO &io)
   // in general, it is a bad idea to have two filehandles for the
   // same file
   cout << "Warning: file handle copy constructor SxBinIO(const SxBinIO&) called.\n";
   open (in.filename, in.mode);
}

SxBinIO::~SxBinIO ()
{
   try {
      close ();
   } catch (SxException e) {
      e.print ();
   }
}


SxBinIO &SxBinIO::operator= (const SxBinIO &in)
{
   if (this == &in)  return *this;

   close ();
   // Using this routine may be a mistake, in most cases due to
   // io = SxBinIO (...)  instead of io.open (...)
   // in general, it is a bad idea to have two filehandles for the
   // same file
   cout << "Warning: file handle copied via SxBinIO::operator=" << endl;
   open (in.filename, in.mode);

   return *this;
}

SxString SxBinIO::createTempName (const SxString &tmpDir)
{
   SxString res = tmpDir + "/tmpsxXXXXXX";  // template
#  ifdef HAVE_MKSTEMP
      char buffer[10240];
      memcpy (buffer, res.ascii(), res.getSize()+1);
      int fd = mkstemp (buffer);
      if (fd < 0)  {
         sxprintf ("Can not create temporary filename %s.\n", res.ascii());
         SX_EXIT;
      }
      res = buffer;
#  else
      cout << "Missing mkstemp. Needs workaround." << endl;
      SX_EXIT;  // not yet implemented
#  endif /* WIN32 */
// --- on non-BSD systems use:
// SxString res;
// res.resize (L_tmpnam);   or tempnam (uses TMPDIR evironment variable)
// tmpnam (res.ascii());

   // --- check, whether temporary filename is valid
   FILE *fpTmp = fopen (res.ascii(), "w");
   if ( !fpTmp )  {
      sxprintf ("Can not write to temporary directory %s.\n", tmpDir.ascii());
      SX_EXIT;
   }
   fclose (fpTmp);
   unlink (res.ascii());

   return res;
}


void SxBinIO::open (const SxString &filename_, Mode mode_)
{
   if (isOpen)  close ();

   bool verbose = false;

   filename = filename_;
   mode     = mode_;

   switch (mode)  {
      case BINARY_READ_ONLY  :     err = nc_open (filename.ascii(),
                                              NC_NOWRITE,
                                              &ncId);

                               if (err != NC_NOERR)  ncError (err);
                               isOpen = true;
                               ncMode = WRITE_HEADER;
                               break;

#ifndef USE_NETCDF4
      // use conventional NetCDF3 output
      case BINARY_WRITE_ONLY : err = nc_create (filename.ascii(),
                                                NC_CLOBBER,
                                                &ncId);
                               if (err != NC_NOERR)  ncError (err);
                               isOpen = true;
                               ncMode = WRITE_HEADER;
                               if (verbose)
                                  cout << "### SxBinIO::open() : BINARY_WRITE_ONLY" << endl;
                               break;
      case BINARY_WRITE_LARGE:
#  ifdef NC_64BIT_OFFSET
                               err = nc_create (filename.ascii(),
                                                NC_CLOBBER | NC_64BIT_OFFSET,
                                                &ncId);
#  else
    /* This is to support netcdf versions < 3.6 for some transition
       time. If you have problems with this part, update your netcdf
       library to 3.6 (see numlib/packages/ )
     */
                               err = nc_create (filename.ascii(),
                                                NC_CLOBBER,
                                                &ncId);
                               cout << endl   
 << "Warning: large netcdf files are not supported by netcdf library." << endl
 << "         The written file may be corrupted. Please update the"    << endl
 << "         library and recompile the S/PHI/ngX package if you need" << endl
 << "         large file (>2GB) support."                              << endl;
#  endif
                               if (err != NC_NOERR)  ncError (err);
                               isOpen = true;
                               ncMode = WRITE_HEADER;
                               if (verbose)
                                  cout << "### SxBinIO::open() : BINARY_WRITE_LARGE" << endl;
                               break;
#else
      // use NetCDF4 output
      case BINARY_WRITE_ONLY :
      case BINARY_WRITE_LARGE:
                                 err = nc_create (filename.ascii(),
                                       NC_NETCDF4,
                                       &ncId);
                                 if (err != NC_NOERR)  ncError (err);
                                 isOpen = true;
                                 ncMode = WRITE_HEADER;
                                 if (verbose)
                                    cout << "### SxBinIO::open() : NetCDF4-serial" << endl;
                                 break;
#endif

      case BINARY_WRITE_PARALLEL:
#ifdef USE_PARALLEL_NETCDF4
         /* NetCDF4/HDF5 and MPIIO */
                               err = nc_create_par (filename.ascii(),
                                        NC_NETCDF4 | NC_MPIIO,
                                        SxLoopMPI::MpiCommWorld(),
                                        SxLoopMPI::MpiInfoNull(),
                                        &ncId);
                               if (err != NC_NOERR)  ncError (err);
                               isOpen = true;
                               ncMode = WRITE_HEADER;
#else
                               cout << "Error: Sx was built without support for parallel NetCDF4 IO." << endl;
                               SX_EXIT;
#endif
                               if (verbose)
                                  cout << "### SxBinIO::open() : NetCDF4-parallel" << endl;
                               break;
                               
                               
      case SCRATCH_READ_WRITE: fp = fopen (filename.ascii(), "w+");
                               isOpen = (fp != NULL);
                               break;
      default                : cout << "Unknown mode in SxBinIO::open." << endl;
                               SX_EXIT;
                               break;
   }
}

void SxBinIO::close ()
{
   if (isOpen)  {
      if (IS_NETCDF_MODE(mode))  {
         err = nc_close (ncId);
         if (err != NC_NOERR)  ncError (err);
      }  else  {
         fclose (fp);
      }
   }
   fp     = NULL;
   ncId   = -1;
   isOpen = false;
   mode   = UNKNOWN;
}


void SxBinIO::ncError (int errNo) const
{
   if (mode == BINARY_READ_ONLY)
      throw SxException (( SxString("I/O Error: ") + nc_strerror(errNo)
                        + " occured during reading S/PHI/nX binary file "
                        + filename
                        ).ascii(),
                        __FILE__, __LINE__);
   else
      throw SxException (( SxString("I/O Error: ") + nc_strerror(errNo)
                        + " occured during writing S/PHI/nX binary file "
                        + filename
                        ).ascii(),
                        __FILE__, __LINE__);
}

void SxBinIO::ncError () const
{
   if (mode == BINARY_READ_ONLY)  {
      throw SxException (( SxString ("I/O Error: ") +  "Corrupt input file! "
                          "Use 'ncdump' to inspect the file:\n"
                        + filename).ascii(),
                        __FILE__, __LINE__);
   } else {
      SX_EXIT;
   }
}

void SxBinIO::ncError (int errNo, const SxString &name) const
{
   throw SxException (( SxString("I/O Error: ") + nc_strerror(errNo)
                     + " occured upon "
                     + (mode == BINARY_READ_ONLY ? "read" : "writ")
                     + "ing '" + name + "' from S/PHI/nX binary file "
                     + filename
                     ).ascii(),
                     __FILE__, __LINE__);
}


void SxBinIO::addDimension (const SxString &dimName, int size) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   int dimId=0;

   // does dimension exist already and is dimension of correct size?
   if (nc_inq_dimid (ncId, dimName.ascii(), &dimId) == NC_NOERR)  {
      size_t n = 0;
      err = nc_inq_dimlen  (ncId, dimId, &n);
      if (err != NC_NOERR)  ncError (err, dimName);
      if ((size_t)size != n)  {
         throw SxException (( SxString("I/O Error: ")
                           + "Overwriting of dimensions "
                           + dimName + " is not allowed.").ascii(),
                           __FILE__, __LINE__);
      }

      // dimension exists already and has the correct size
      return;
   }
   err = nc_def_dim (ncId, dimName.ascii(), size, &dimId);
   if (err != NC_NOERR)  ncError (err, dimName);
}


int SxBinIO::getDimension (const SxString &dimName) const
{
   int dimId = 0;
   err = nc_inq_dimid (ncId, dimName.ascii(), &dimId);
   if (err != NC_NOERR)  ncError (err, dimName);

   size_t size = 0;
   err = nc_inq_dimlen  (ncId, dimId, &size);
   if (err != NC_NOERR)  ncError (err, dimName);

   return (int)size;
}


bool SxBinIO::contains (const SxString &varName) const
{
   int varId = 0;
   // variable?
   err = nc_inq_varid (ncId, varName.ascii(), &varId);
   if (err == NC_NOERR)  return true;
   // global attribute?
   err = nc_inq_attid (ncId, NC_GLOBAL, varName.ascii (), &varId);
   if (err != NC_NOERR)  return false;
   else                  return true;
}


bool SxBinIO::containsDim (const SxString &dimName) const
{
   int dimId = 0;
   err = nc_inq_dimid (ncId, dimName.ascii(), &dimId);
   if (err != NC_NOERR)  return false;
   else                  return true;
}


void SxBinIO::setMode (int  m) const
{
   SX_CHECK (ncMode == WRITE_HEADER, ncMode);
   ncMode = m;
   if (m == WRITE_DATA)  {
      err = nc_enddef (ncId);
      if (err != NC_NOERR)  ncError (err);
   }
}


int SxBinIO::getSize (const SxString &name) const
{
   size_t size;

   // --- retrieve dimension
   int dimId=0;
   err = nc_inq_dimid (ncId, name.ascii(), &dimId);
   if (err != NC_NOERR)  ncError (err, name);

   err = nc_inq_dimlen (ncId, dimId, &size);
   if (err != NC_NOERR)  ncError (err, name);

   return int(size);
}


void SxBinIO::read (const SxString &varName, SxString *val) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   char buffer[10240];
   size_t len = 0;
   err = nc_inq_attlen (ncId, NC_GLOBAL, varName.ascii (), &len);
   char *buf = (len < 10239) ? buffer : new char[len+1];

   err = nc_get_att_text (ncId, NC_GLOBAL, varName.ascii(), buf);
   if (err != NC_NOERR)    ncError (err, varName);
   buf[len] = 0; // trailing 0

   *val = SxString(buffer);
   if (buf != buffer) delete [] buf;
}



template<class T>
static bool isValid (const SxVecRef<T> &vec)
{
   SX_CHECK (vec.getSize () > 0, vec.getSize ());
   for (int i = 0; i < vec.getSize (); i++)  {
      if (!::isValid (vec(i))) {
         cout << "ERROR: vector is corrupt!" << endl;
         return false;
      }
   }
   return true;
}

void SxBinIO::read (const SxString &varName, SxList<int> *val, int nElem) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   int id = 0;
   size_t start[] = {0};
   size_t count[] = {(size_t)nElem};

   err = nc_inq_varid (ncId, varName.ascii(), &id);
   if (err != NC_NOERR)  ncError (err, varName);

   SxVector<int> buffer (nElem);

   err = nc_get_vara_int (ncId, id, start, count, buffer.elements);
   if (err != NC_NOERR)    ncError (err, varName);
   //if (!isValid (buffer))  ncError (); // int is always valid

   val->removeAll ();
   for (ssize_t i=0; i < nElem; i++)  *val << buffer(i);
}



void SxBinIO::write (const SxString &varName, const SxString &val) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   if (ncMode == WRITE_HEADER)  {
      err = nc_put_att_text (ncId, NC_GLOBAL, varName.ascii(),
                             val.getSize(), val.ascii());
      if (err != NC_NOERR)  ncError (err, varName);

   }
}


void SxBinIO::write (const SxString &varName, double val) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   if (ncMode == WRITE_HEADER)  {
      err = nc_put_att_double (ncId, NC_GLOBAL, varName.ascii(),
                               NC_DOUBLE, 1, &val);
      if (err != NC_NOERR)  ncError (err, varName);
   }
}

void SxBinIO::read (const SxString &varName, double *val) const
{
   SX_CHECK (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   size_t len;

   err = nc_inq_attlen (ncId, NC_GLOBAL, varName.ascii (), &len);
   if (err != NC_NOERR)  ncError (err, varName);
   if (len != 1)
      throw SxException (( "Type missfit for " + varName
                         + ": Wanted number, found vector (" + double(len)
                         + " elements).").ascii());

   err = nc_get_att_double (ncId, NC_GLOBAL, varName.ascii (), val);
   if (err != NC_NOERR)  ncError (err, varName);
}

void SxBinIO::read (const SxString &varName, SxVector3<int> *val) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   int id = 0;
   size_t start[] = {0};
   size_t count[] = {3};

   err = nc_inq_varid (ncId, varName.ascii(), &id);
   if (err != NC_NOERR)  ncError (err, varName);

   err = nc_get_vara_int (ncId, id, start, count, val->v);
   if (err != NC_NOERR)  ncError (err, varName);

}


void SxBinIO::write (const SxString &varName,
                  const SxVector3<int> &val,
                  const SxString &dimName) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   int id = 0;

   // --- retrieve dimension
   int dimId=0;
   err = nc_inq_dimid (ncId, dimName.ascii(), &dimId);
   if (err != NC_NOERR)  ncError (err, dimName);

   if (ncMode == WRITE_HEADER)  {

      int dims[] = {dimId};

      err = nc_def_var (ncId, varName.ascii(), NC_INT, 1, dims, &id);
      if (err != NC_NOERR)  ncError (err, varName);

   }  else  {

      err = nc_inq_varid (ncId, varName.ascii(), &id);
      if (err != NC_NOERR)  ncError (err, varName);

      size_t start[] = {0};
      size_t count[] = {3};

      err = nc_put_vara_int (ncId, id, start, count,
                            (const int *)val.v);
      if (err != NC_NOERR)  ncError (err, varName);

   }
}



void SxBinIO::write (const SxString &varName,
                  const SxVector3<double> &val,
                  const SxString &dimName) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   int id = 0;

   // --- retrieve dimension
   int dimId=0;
   err = nc_inq_dimid (ncId, dimName.ascii(), &dimId);
   if (err != NC_NOERR)  ncError (err, dimName);

   if (ncMode == WRITE_HEADER)  {

      int dims[] = {dimId};

      err = nc_def_var (ncId, varName.ascii(), NC_DOUBLE, 1, dims, &id);
      if (err != NC_NOERR)  ncError (err, varName);

   }  else  {

      err = nc_inq_varid (ncId, varName.ascii(), &id);
      if (err != NC_NOERR)  ncError (err, varName);

      size_t start[] = {0};
      size_t count[] = {3};

      err = nc_put_vara_double (ncId, id, start, count,
                               (const double *)val.v);
      if (err != NC_NOERR)  ncError (err, varName);

   }
}


void SxBinIO::read (const SxString &varName, SxVector3<double> *val) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   int id = 0;
   size_t start[] = {0};
   size_t count[] = {3};

   err = nc_inq_varid (ncId, varName.ascii (), &id);
   if (err != NC_NOERR)  ncError (err, varName);

   err = nc_get_vara_double (ncId, id, start, count, val->v);
   if (err != NC_NOERR)  ncError (err, varName);
}


void SxBinIO::write (const SxString &varName,
                  const SxMatrix3<double> &val,
                  const SxString &dimName) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   int id = 0;

   // --- retrieve dimension
   int dimId=0;
   err = nc_inq_dimid (ncId, dimName.ascii(), &dimId);
   if (err != NC_NOERR)  ncError (err, dimName);

   if (ncMode == WRITE_HEADER)  {

      int dims[] = {dimId, dimId};

      err = nc_def_var (ncId, varName.ascii(), NC_DOUBLE, 2, dims, &id);
      if (err != NC_NOERR)  ncError (err, varName);

   }  else  {

      err = nc_inq_varid (ncId, varName.ascii(), &id);
      if (err != NC_NOERR)  ncError (err, varName);

      size_t start[] = {0, 0};
      size_t count[] = {3, 3};

      err = nc_put_vara_double (ncId, id, start, count,
                               (const double *)val.m);
      if (err != NC_NOERR)  ncError (err, varName);

   }
}


void SxBinIO::read (const SxString &varName, SxMatrix3<double> *val) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   int id = 0;
   size_t start[] = {0, 0};
   size_t count[] = {3, 3};

   err = nc_inq_varid (ncId, varName.ascii(), &id);
   if (err != NC_NOERR)  ncError (err, varName);

   err = nc_get_vara_double (ncId, id, start, count, (double *)&val->m);
   if (err != NC_NOERR)  ncError (err, varName);
}

void SxBinIO::write (const SxString &varName,
                  const SxArray<SxMatrix3<double> > &val,
                  const SxString &dim1Name,
                  const SxString &dim2Name) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   int id = 0;

   // --- retrieve dimension
   int dim1Id=0, dim2Id = 0;
   err = nc_inq_dimid (ncId, dim1Name.ascii (), &dim1Id);
   if (err != NC_NOERR)  ncError (err, dim1Name);
   err = nc_inq_dimid (ncId, dim2Name.ascii (), &dim2Id);
   if (err != NC_NOERR)  ncError (err, dim2Name);

   if (ncMode == WRITE_HEADER)  {

      int dims[] = {dim1Id, dim2Id, dim2Id};

      err = nc_def_var (ncId, varName.ascii (), NC_DOUBLE, 3, dims, &id);
      if (err != NC_NOERR)  ncError (err, varName);

   }  else  {

      err = nc_inq_varid (ncId, varName.ascii (), &id);
      if (err != NC_NOERR)  ncError (err, varName);

      size_t start[] = {0, 0, 0};
      size_t count[] = {1, 3, 3};

      // --- write one SxMatrix3 at a time
      for (ssize_t i = 0; i < val.getSize (); i++)  {
         start[0] = (size_t)i;
         err = nc_put_vara_double (ncId, id, start, count,
                                   (const double *)val(i).m);
         if (err != NC_NOERR)  ncError (err, varName);
      }

   }
}


void SxBinIO::read (const SxString &varName,
                 SxArray<SxMatrix3<double> > *val) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   int id = 0;
   size_t start[] = {0, 0, 0};
   size_t count[] = {1, 3, 3};

   err = nc_inq_varid (ncId, varName.ascii (), &id);
   if (err != NC_NOERR)  ncError (err, varName);

   for (ssize_t i = 0; i < val->getSize (); i++)  {
      start[0] = (size_t)i;
      err = nc_get_vara_double (ncId, id, start, count, (double *)&(*val)(i).m);
      if (err != NC_NOERR)  ncError (err, varName);
   }
}


// TODO   Check support for NetCDF4 carefully.
void SxBinIO::writeVec (const SxString &varName,
                  const SxVecRef<int> &val,
                  const SxString &dimName,
                  int offset, int nelem, int localOffset) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);
   SX_CHECK ((offset >= 0), offset);           // negative offsets are used
   SX_CHECK ((localOffset >= 0), localOffset); // during debugging

   int id = 0;

   // --- retrieve dimension
   int dimId=0;
   err = nc_inq_dimid (ncId, dimName.ascii(), &dimId);
   if (err != NC_NOERR)  ncError (err, dimName);

   if (ncMode == WRITE_HEADER)  {

      int dims[] = {dimId};

      err = nc_def_var (ncId, varName.ascii(), NC_INT, 1, dims, &id);
      if (err != NC_NOERR && !offset)  ncError (err, varName);

   }  else  {

      err = nc_inq_varid (ncId, varName.ascii(), &id);
      if (err != NC_NOERR)  ncError (err, varName);

      size_t start[] = {(size_t)offset};
      size_t count[] = {size_t(nelem == 0 ? val.getSize() : nelem)};

      int * rowPtr = val.elements;
      for (int i = 0; i < localOffset; i++) rowPtr++;

      err = nc_put_vara_int (ncId, id, start, count, rowPtr);
      if (err != NC_NOERR)  ncError (err, varName);

   }
}


void SxBinIO::readVec (const SxString &varName,
                 SxVecRef<int> *val, int nElem, int offset) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   int id = 0;
   size_t start[] = {(size_t)offset};
   size_t count[] = {(size_t)nElem};

   err = nc_inq_varid (ncId, varName.ascii(), &id);
   if (err != NC_NOERR)  ncError (err, varName);

   // val->resize (nElem);
   SX_CHECK (val->getSize () >= nElem, val->getSize (), nElem);

   err = nc_get_vara_int (ncId, id, start, count, val->elements);
   if (err != NC_NOERR)  ncError (err, varName);
   // if (!isValid (*val))  ncError (); // int is always valid
}



// TODO   check support for NetCDF4
void SxBinIO::writeVec (const SxString &varName,
                  const SxVecRef<double> &val,
                  const SxString &dimName,
                  int offset, int nelem, int localOffset) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);
   SX_CHECK ((offset >= 0), offset);
   SX_CHECK ((localOffset >= 0), localOffset);

   int id = 0;

   // --- retrieve dimension
   int dimId=0;
   err = nc_inq_dimid (ncId, dimName.ascii(), &dimId);
   if (err != NC_NOERR)  ncError (err, dimName);

   if (ncMode == WRITE_HEADER)
   {
      int dims[] = {dimId};

      err = nc_def_var (ncId, varName.ascii(), NC_DOUBLE, 1, dims, &id);
      if (err != NC_NOERR && !offset)  ncError (err, varName);
   }
   else
   {
      err = nc_inq_varid (ncId, varName.ascii(), &id);
      if (err != NC_NOERR)  ncError (err, varName);

      size_t start[] = {(size_t)offset};
      size_t count[] = {size_t(nelem == 0 ? val.getSize() : nelem)};

      double * rowPtr = val.elements;
      for (int i = 0; i < localOffset; i++) rowPtr++;

      err = nc_put_vara_double (ncId, id, start, count, rowPtr);
      if (err != NC_NOERR)  ncError (err, varName);
   }
}


void SxBinIO::readVec (const SxString &varName,
                 SxVecRef<double> *val, int nElem, int offset) const
{
   SX_CHECK     (val);
   SX_CHECK (mode == BINARY_READ_ONLY, mode);

   int id = 0;
   size_t start[] = {(size_t)offset};
   size_t count[] = {(size_t)nElem};

   err = nc_inq_varid (ncId, varName.ascii(), &id);
   if (err != NC_NOERR)  ncError (err, varName);

   // val->resize (nElem);
   SX_CHECK (val->getSize () >= nElem, val->getSize (), nElem);

   err = nc_get_vara_double (ncId, id, start, count, val->elements);
   if (err != NC_NOERR)  ncError (err, varName);
   if (!isValid (*val))  ncError ();
}


// TODO   Check support for NetCDF4 parallel IO
void SxBinIO::writeVec (const SxString &varName,
                  const SxVecRef<SxComplex16> &val,
                  const SxString &dimName, int offset) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode));
   int idRe=0, idIm=0;

   // --- retrieve dimension
   int dimId=0;
   err = nc_inq_dimid (ncId, dimName.ascii(), &dimId);
   if (err != NC_NOERR)  ncError (err, dimName);

   if (ncMode == WRITE_HEADER)  {

      int dims[] = {dimId};

      err = nc_def_var (ncId, (varName+".re").ascii(),
            NC_DOUBLE, 1, dims, &idRe);
      if (err != NC_NOERR && !offset)  ncError (err, varName);
      err = nc_def_var (ncId, (varName+".im").ascii(),
            NC_DOUBLE, 1, dims, &idIm);
      if (err != NC_NOERR && !offset)  ncError (err, varName);

   }  else  {

      err = nc_inq_varid (ncId, (varName+".re").ascii(), &idRe);
      if (err != NC_NOERR)  ncError (err, varName);
      err = nc_inq_varid (ncId, (varName+".im").ascii(), &idIm);
      if (err != NC_NOERR)  ncError (err, varName);

      size_t start[] = {(size_t)offset};
      size_t count[] = {(size_t)val.getSize()};

      err = nc_put_vara_double (ncId, idRe, start, count,
                                SxVector<double>(val.real()).elements);
      if (err != NC_NOERR)  ncError (err, varName);
      err = nc_put_vara_double (ncId, idIm, start, count,
                                SxVector<double>(val.imag()).elements);
      if (err != NC_NOERR)  ncError (err, varName);

   }
}


void SxBinIO::readVec (const SxString &varName,
                 SxVecRef<SxComplex16> *val, int nElem, int offset) const
{
   SX_CHECK      (val);
   SX_CHECK (val->getSize() == nElem, val->getSize(), nElem);
   SX_CHECK  (mode == BINARY_READ_ONLY, mode);

   int idRe=0, idIm=0;
   size_t start[] = {(size_t)offset};
   size_t count[] = {(size_t)nElem};

   err = nc_inq_varid (ncId, (varName+".re").ascii(), &idRe);
   if (err != NC_NOERR)  ncError (err, varName);
   err = nc_inq_varid (ncId, (varName+".im").ascii(), &idIm);
   if (err != NC_NOERR)  ncError (err, varName);

// val->resize (nElem);

   SxVector<double> vecRe(nElem), vecIm(nElem);

   err = nc_get_vara_double (ncId, idRe, start, count, vecRe.elements);
   if (err != NC_NOERR)  ncError (err, varName);
   err = nc_get_vara_double (ncId, idIm, start, count, vecIm.elements);
   if (err != NC_NOERR)  ncError (err, varName);

   SxVectorReIm<false>::real (*val) = vecRe;
   SxVectorReIm<false>::imag (*val) = vecIm;
   if (!isValid (*val))  ncError();
}

// --- SxVector: matrices
void SxBinIO::writeMat (const SxString &varName,
                  const SxVecRef<double> &val,
                  const SxString &dim1Name, const SxString &dim2Name,
                  int offsetRow, int offsetCol) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   int id = 0;

   // --- retrieve dimension
   int dim1Id=0, dim2Id=0;
   err = nc_inq_dimid (ncId, dim1Name.ascii(), &dim1Id);
   if (err != NC_NOERR)  ncError (err, dim1Name);
   err = nc_inq_dimid (ncId, dim2Name.ascii(), &dim2Id);
   if (err != NC_NOERR)  ncError (err, dim2Name);

   if (ncMode == WRITE_HEADER)  {
      int dims[] = {dim1Id, dim2Id};

      err = nc_def_var (ncId, varName.ascii(), NC_DOUBLE, 2, dims, &id);
      if (err != NC_NOERR)  ncError (err, varName);

   }  else  {

      err = nc_inq_varid (ncId, varName.ascii(), &id);
      if (err != NC_NOERR)  ncError (err, varName);

      size_t size1, size2;
      err = nc_inq_dimlen (ncId, dim1Id, &size1);
      if (err != NC_NOERR)  ncError (err, dim1Name);
      SX_CHECK (size_t(offsetRow + val.getNRows ()) <= size1,
                offsetRow, val.getNRows (), size1);
      err = nc_inq_dimlen (ncId, dim2Id, &size2);
      if (err != NC_NOERR)  ncError (err, dim2Name);
      SX_CHECK (size_t(offsetCol + val.getNCols ()) <= size2,
                offsetCol, val.getNCols (), size2);

      size_t start[] = {(size_t)offsetRow, (size_t)offsetCol};
      size_t count[] = {(size_t)val.getNRows (), 1};

      // write matrix columnwise, because matrix is stored columnwise
      for (ssize_t col = 0; col < val.getNCols (); col++)  {
         err = nc_put_vara_double (ncId, id, start, count,
               val.colRef(col).elements);
         if (err != NC_NOERR)  ncError (err, varName);
         start[1]++; // next column
      }
   }
}

void SxBinIO::readMat (const SxString &varName,
                 SxVecRef<double> *val,
                 int nElemRow, int nElemCol,
                 int offsetRow, int offsetCol) const
{
   SX_CHECK      (val);
   SX_CHECK (val->getNRows () == nElemRow, val->getNRows (), nElemRow);
   SX_CHECK (val->getNCols () == nElemCol, val->getNCols (), nElemCol);
   SX_CHECK  (mode == BINARY_READ_ONLY, mode);


   int id = 0;
   size_t start[] = {(size_t)offsetRow, (size_t)offsetCol};
   size_t count[] = {(size_t)nElemRow, 1}; // read 1 column at a time

   err = nc_inq_varid (ncId, varName.ascii (), &id);
   if (err != NC_NOERR)  ncError (err, varName);

   // read matrix columnwise, because matrix is stored columnwise
   for (ssize_t col = 0; col < nElemCol; col++)  {
      err = nc_get_vara_double (ncId, id, start, count,
                                val->colRef(col).elements);
      if (err != NC_NOERR)  ncError (err, varName);
      start[1]++; // next column
   }
}

void SxBinIO::writeMat (const SxString &varName,
                  const SxVecRef<SxComplex16> &val,
                  const SxString &dim1Name, const SxString &dim2Name,
                  int offsetRow, int offsetCol) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   int idRe=0, idIm=0;

   // --- retrieve dimensions
   int dimRowId=0, dimColId=0;
   err = nc_inq_dimid (ncId, dim1Name.ascii(), &dimRowId);
   if (err != NC_NOERR)  ncError (err, dim1Name);
   err = nc_inq_dimid (ncId, dim2Name.ascii(), &dimColId);
   if (err != NC_NOERR)  ncError (err, dim2Name);

   if (ncMode == WRITE_HEADER)  {
      int dims[] = {dimRowId, dimColId};

      err = nc_def_var (ncId, (varName+".re").ascii(),
            NC_DOUBLE, 2, dims, &idRe);
      if (err != NC_NOERR)  ncError (err, varName);
      err = nc_def_var (ncId, (varName+".im").ascii(),
            NC_DOUBLE, 2, dims, &idIm);
      if (err != NC_NOERR)  ncError (err, varName);

   }  else  {

      err = nc_inq_varid (ncId, (varName+".re").ascii(), &idRe);
      if (err != NC_NOERR)  ncError (err, varName);
      err = nc_inq_varid (ncId, (varName+".im").ascii(), &idIm);
      if (err != NC_NOERR)  ncError (err, varName);

      size_t sizeRow, sizeCol;
      err = nc_inq_dimlen (ncId, dimRowId, &sizeRow);
      if (err != NC_NOERR)  ncError (err, dim1Name);
      SX_CHECK (size_t(offsetRow + val.getNRows ()) <= sizeRow,
                offsetRow, val.getNRows (), sizeRow);
      err = nc_inq_dimlen (ncId, dimColId, &sizeCol);
      if (err != NC_NOERR)  ncError (err, dim2Name);
      SX_CHECK (size_t(offsetCol + val.getNCols ()) <= sizeCol,
                offsetCol, val.getNCols (), sizeCol);

      size_t start[] = {(size_t)offsetRow, (size_t)offsetCol};
      size_t count[] = {(size_t)val.getNRows (), 1};

      // --- write matrix columnwise, because matrices are stored columnwise
      for (ssize_t col = 0; col < val.getNCols (); col++)  {
         err = nc_put_vara_double (ncId, idRe, start, count,
               SxVector<double>(val.colRef(col).real()).elements);
         if (err != NC_NOERR)  ncError (err, varName);
         err = nc_put_vara_double (ncId, idIm, start, count,
               SxVector<double>(val.colRef(col).imag()).elements);
         if (err != NC_NOERR)  ncError (err, varName);
         start[1]++;  // next column
      }
   }
}

void SxBinIO::readMat (const SxString &varName,
                 SxVecRef<SxComplex16> *val,
                 int nElemRow, int nElemCol,
                 int offsetRow, int offsetCol) const
{
   SX_CHECK      (val);
   SX_CHECK (val->getNRows () == nElemRow, val->getNRows (), nElemRow);
   SX_CHECK (val->getNCols () == nElemCol, val->getNCols (), nElemCol);

   SX_CHECK  (mode == BINARY_READ_ONLY, mode);

   int idRe=0, idIm=0;
   size_t start[] = {(size_t)offsetRow, (size_t)offsetCol};
   size_t count[] = {(size_t)nElemRow, 1};  // read 1 column at a time

   err = nc_inq_varid (ncId, (varName+".re").ascii(), &idRe);
   if (err != NC_NOERR)  ncError (err, varName);
   err = nc_inq_varid (ncId, (varName+".im").ascii(), &idIm);
   if (err != NC_NOERR)  ncError (err, varName);

   SxVector<double>    colRe(nElemRow), colIm(nElemRow);

   // --- read matrices columnwise, because they are stored columnwise
   for (ssize_t col = 0; col < nElemCol; col++)  {
      err = nc_get_vara_double (ncId, idRe, start, count, colRe.elements);
      if (err != NC_NOERR)  ncError (err, varName);
      err = nc_get_vara_double (ncId, idIm, start, count, colIm.elements);
      if (err != NC_NOERR)  ncError (err, varName);

      SxVecRef<SxComplex16> colRef = val->colRef(col);
      SxVectorReIm<false>::real(colRef) = colRe;
      SxVectorReIm<false>::imag(colRef) = colIm;
      start[1]++;  // next column
   }

   if (!isValid (*val))  ncError();
}


// ---
void SxBinIO::writeMat (const SxString &varName,
                  const SxVecRef<int> &val,
                  const SxString &dimRowName, const SxString &dimColName,
                  int offsetRow, int offsetCol) const
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);

   int id = 0;

   // --- retrieve dimension
   int dimRowId=0, dimColId=0;
   err = nc_inq_dimid (ncId, dimRowName.ascii (), &dimRowId);
   if (err != NC_NOERR)  ncError (err, dimRowName);
   err = nc_inq_dimid (ncId, dimColName.ascii (), &dimColId);
   if (err != NC_NOERR)  ncError (err, dimColName);

   if (ncMode == WRITE_HEADER)  {

      int dims[] = {dimRowId, dimColId};

      err = nc_def_var (ncId, varName.ascii (), NC_INT, 2, dims, &id);
      if (err != NC_NOERR)  ncError (err, varName);

   }  else  {

      err = nc_inq_varid (ncId, varName.ascii (), &id);
      if (err != NC_NOERR)  ncError (err, varName);

      size_t sizeRow, sizeCol;
      err = nc_inq_dimlen (ncId, dimRowId, &sizeRow);
      if (err != NC_NOERR)  ncError (err, dimRowName);
      SX_CHECK (size_t(offsetRow + val.getNRows ()) <= sizeRow,
                offsetRow, val.getNRows (), sizeRow);
      err = nc_inq_dimlen (ncId, dimColId, &sizeCol);
      if (err != NC_NOERR)  ncError (err, dimColName);
      SX_CHECK (size_t(offsetCol + val.getNCols ()) <= sizeCol,
                offsetCol, val.getNCols (), sizeCol);

      size_t start[] = {(size_t)offsetRow, (size_t)offsetCol};
      size_t count[] = {(size_t)val.getNRows (), 1}; // 1 column at a time

      // write matrix columnwise, because matrix is stored columnwise
      for (ssize_t col = 0; col < val.getNCols (); col++)  {
         err = nc_put_vara_int (ncId, id, start, count,
                                val.colRef(col).elements);
         if (err != NC_NOERR)  ncError (err, varName);
         start[1]++; // next column
      }

   }
}

// TODO   Check support for parallel IO
void SxBinIO::writeRows (const SxString &varName,
                  const SxVecRef<double> &val,
                  const SxString &dim1Name, const SxString &dim2Name,
                  int rowOffset, int rowNelem, int localOffset) const
{
   SX_CHECK ((rowOffset >= 0), rowOffset);
   SX_CHECK ((localOffset >= 0), localOffset);

   int id = 0;


   if ( IS_NETCDF_WRITE(mode) ) {
      // --- retrieve dimension
      int dim1Id=0, dim2Id=0;
      err = nc_inq_dimid (ncId, dim1Name.ascii(), &dim1Id);
      if (err != NC_NOERR)  ncError (err, dim1Name);
      err = nc_inq_dimid (ncId, dim2Name.ascii(), &dim2Id);
      if (err != NC_NOERR)  ncError (err, dim2Name);

      if (ncMode == WRITE_HEADER)  {

         int dims[] = {dim1Id, dim2Id};

         err = nc_def_var (ncId, varName.ascii(), NC_DOUBLE, 2, dims, &id);
         if (err != NC_NOERR)  ncError (err, varName);

      }  else  {

         err = nc_inq_varid (ncId, varName.ascii(), &id);
         if (err != NC_NOERR)  ncError (err, varName);

         size_t size1, size2;
         err = nc_inq_dimlen (ncId, dim1Id, &size1);
         if (err != NC_NOERR)  ncError (err, dim1Name);
         err = nc_inq_dimlen (ncId, dim2Id, &size2);
         if (err != NC_NOERR)  ncError (err, dim2Name);

         // rowOffset and rowNelem are default parameters initialized to zero
         size_t start[] = {(size_t)rowOffset, 0};
         size_t count[] = {size_t(rowNelem == 0 ? val.getNRows () : rowNelem), 1};
         SX_CHECK (ssize_t(localOffset + count[0]) <= val.getNRows (),
                   localOffset, count[0], val.getNRows ());
         SX_CHECK (rowOffset + count[0] <= size1, rowOffset, count[0], size1);

         // write matrix columnwise, because matrix is stored columnwise
         for (ssize_t col = 0; col < val.getNCols (); col++)
         {
            double * rowPtr = val.colRef(col).elements;
            rowPtr += localOffset;

            err = nc_put_vara_double (ncId, id, start, count, rowPtr );
            if (err != NC_NOERR)  ncError (err, varName);

            start[1]++; // next column
         }
      }
   }

}


void SxBinIO::readMat (const SxString &varName,
                 SxVecRef<int> *val,
                 int nElemRow, int nElemCol,
                 int offsetRow, int offsetCol) const
{
   SX_CHECK      (val);
   SX_CHECK (val->getNRows () == nElemRow, val->getNRows (), nElemRow);
   SX_CHECK (val->getNCols () == nElemCol, val->getNCols (), nElemCol);
   SX_CHECK  (mode == BINARY_READ_ONLY, mode);


   int id = 0;
   size_t start[] = {(size_t)offsetRow, (size_t)offsetCol};
   size_t count[] = {(size_t)nElemRow, 1}; // read 1 column at a time

   err = nc_inq_varid (ncId, varName.ascii (), &id);
   if (err != NC_NOERR)  ncError (err, varName);

   // read matrix columnwise, because matrix is stored columnwise
   for (ssize_t col = 0; col < nElemCol; col++)  {
      err = nc_get_vara_int (ncId, id, start, count, val->colRef(col).elements);
      if (err != NC_NOERR)  ncError (err, varName);
      start[1]++; // next column
   }
}



void SxBinIO::writeMesh (const SxVecRef<double>  &mesh,
                      const SxMatrix3<double> &cell,
                      const SxVector3<int>    &dim)
{
   SxArray<SxVector<double> > tmpArray (1);
   tmpArray(0) = mesh;
   writeMesh (tmpArray, cell, dim);
}

void SxBinIO::writeMesh (const SxArray<SxVector<double> > &meshes,
                      const SxMatrix3<double>          &cell,
                      const SxVector3<int>             &dim) const
{
   int nMeshes = int(meshes.getSize());
   int nElem = int(meshes(0).getSize ());
   SX_CHECK (nElem == dim.product (), nElem, dim.product ());

   // --- write header, create dimensions
   addDimension ("nMeshes",  nMeshes);
   addDimension ("xyz",      3);
   addDimension ("meshSize", nElem);

   // --- write data
   writeVec ("dim", SxVector<int>(SxList<int> () << dim(0) << dim(1) << dim(2)),
          "xyz");
   write ("cell", cell, "xyz");


   SxString varName;

   for (int i=0; i < nMeshes; i++)  {
      varName = "mesh";
      if (nMeshes > 1)  varName += SxString("-") + i;
      writeVec (varName, meshes(i), "meshSize");
   }
}

SxArray<SxVector<double> > SxBinIO::readMesh (SxMatrix3<double> *cellPtr,
                                             SxVector3<int>    *dimPtr) const
{
   int nMeshes = getDimension ("nMeshes");
   int nElem   = getDimension ("meshSize");
   SxVector3<int> dim;  read("dim", &dim);

   SxArray<SxVector<double> >  meshes (nMeshes);
   SxString file;
   for (int i=0; i < nMeshes; i++)  {
      file = "mesh";
      if (nMeshes > 1)  file += SxString("-") + i;

      meshes(i).resize (nElem);
      readVec (file, &meshes(i), nElem);
      if (!isValid (meshes(i)))  ncError();
   }

   if (cellPtr)  read ("cell", cellPtr);
   if (dimPtr)   read ("dim",  dimPtr);

   return meshes;
}


void SxBinIO::addDoubleVar (const SxString &varName,
                         const SxString &dim1Name,
                         const SxString &dim2Name)
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);
   SX_CHECK (ncMode == WRITE_HEADER);

   // --- retrieve dimension
   int dim1Id=0, dim2Id = 0;
   err = nc_inq_dimid (ncId, dim1Name.ascii (), &dim1Id);
   if (err != NC_NOERR)  ncError (err, dim1Name);
   err = nc_inq_dimid (ncId, dim2Name.ascii (), &dim2Id);
   if (err != NC_NOERR)  ncError (err, dim2Name);

   int id = 0;
   int dims[] = {dim1Id, dim2Id};

   err = nc_def_var (ncId, varName.ascii (), NC_DOUBLE, 2, dims, &id);
   if (err != NC_NOERR)  ncError (err, varName);

}

void SxBinIO::writeRow (const SxString &varName,
                     const SxVecRef<double> &val,
                     int ir,
                     int offset)
{
   SX_CHECK (IS_NETCDF_WRITE(mode), mode);
   SX_CHECK (ncMode == WRITE_DATA);

   int id = 0;
   err = nc_inq_varid (ncId, varName.ascii (), &id);
   if (err != NC_NOERR)  ncError (err, varName);

   size_t start[] = {(size_t)ir, (size_t)offset};
   size_t count[] = {1, (size_t)val.getSize ()};

   // --- write column
   err = nc_put_vara_double (ncId, id, start, count, val.elements);
   if (err != NC_NOERR)  ncError (err, varName);
}

void SxBinIO::readRow (const SxString &varName,
                    SxVecRef<double> *val,
                    int ir,
                    int offset) const
{
   SX_CHECK  (mode == BINARY_READ_ONLY, mode);
   SX_CHECK (val);

   int id = 0;
   err = nc_inq_varid (ncId, varName.ascii (), &id);
   if (err != NC_NOERR)  ncError (err, varName);

   size_t start[] = {(size_t)ir, (size_t)offset};
   size_t count[] = {1, (size_t)val->getSize ()};

   // --- write column
   err = nc_get_vara_double (ncId, id, start, count, val->elements);
   if (err != NC_NOERR)  ncError (err, varName);
}
