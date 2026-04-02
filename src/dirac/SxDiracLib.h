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

#ifndef SX_EXPORT_DIRAC
#ifdef WIN32
#  if defined(_EXPORT_sxdirac)
#     define SX_EXPORT_DIRAC __declspec(dllexport)
#  else
#     define SX_EXPORT_DIRAC __declspec(dllimport)
#  endif
#else
#  define SX_EXPORT_DIRAC
#endif
#endif

#define NONE_M char(0x80)
