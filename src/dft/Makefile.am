# ---------------------------------------------------------------------------
#
#      The ab-initio based multiscale library
#
#                  S / P H I / n X
#
#      Copyright:  Max-Planck-Institut fuer Eisenforschung GmbH
#                  40237 Duesseldorf, Germany
#
#      Contact:    https://sxrepo.mpie.de
#      Authors:    see sphinx/AUTHORS
#      License:    see sphinx/LICENSE
#
# ---------------------------------------------------------------------------

NAME = "sxdft"
MODE = "libOnly"
VERSION = "1:0:0"

SX_DEPLIBS  = $(SXLGPL) $(SXDIRAC) $(SXGEOM) $(SXIO2) $(SXIO) $(SXFS) $(SXIPC)
SX_DEPLIBS += $(SXMATH) $(SXUTIL) $(SXNONSTD) $(SXMPI)
NUM_DEPLIBS = $(FFT) $(BLAS)

EXTRA_DIST=d3Data.h

#SUBDIRS =
#if COND_BUILD_EXAMPLES
#   SUBDIRS += examples
#endif


# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist

AM_CPPFLAGS += @CXX_PEDANTIC@
