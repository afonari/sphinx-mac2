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

#ifndef _SX_ADD_UTIL_H_
#define _SX_ADD_UTIL_H_


#ifdef WIN32
#  if defined(_EXPORT_sxaddutil)
#     define SX_EXPORT_ADD_UTIL __declspec(dllexport)
#  else
#     define SX_EXPORT_ADD_UTIL __declspec(dllimport)
#  endif
#else
#  define SX_EXPORT_ADD_UTIL
#endif

#endif /* _SX_ADD_UTIL_H_ */
