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

#ifndef _SX_PARALLEL_INIT_H_
#define _SX_PARALLEL_INIT_H_

#include <SxIO2.h>
#include <SxParallelHierarchy.h>

/** \brief Write out the current parallel hierarchy in sx format */
void SX_EXPORT_IO2 writeParallelHierarchy(const SxString &filename);

/** Setup parallel hierarchy from sx file */
void SX_EXPORT_IO2 buildHierarchy(const SxString &hierarchyFile,
                                  SxParallelHierarchy &ph);
#endif // _SX_PARALLEL_INIT_H_
