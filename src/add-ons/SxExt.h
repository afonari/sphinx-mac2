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

#ifndef _SX_EXT_H_
#define _SX_EXT_H_


#ifdef WIN32
#  if defined(_EXPORT_sxext)
#     define SX_EXPORT_EXT __declspec(dllexport)
#  else
#     define SX_EXPORT_EXT __declspec(dllimport)
#  endif
#else
#  define SX_EXPORT_EXT
#endif

#endif /* _SX_EXT_H_ */
