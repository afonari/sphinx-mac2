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

NAME = "sxdirac"
MODE = "libOnly"
VERSION = "1:0:0"

#SUBDIRS =
#if COND_BUILD_EXAMPLES
#   SUBDIRS += examples
#endif


SX_DEPLIBS  = $(SXIO2) $(SXGEOM) $(SXIO) $(SXIPC) $(SXMATH) $(SXFS) $(SXUTIL) $(SXNONSTD)
SX_DEPLIBS += $(SXMPI)
NUM_DEPLIBS = $(NETCDF)

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist

AM_CPPFLAGS += @CXX_PEDANTIC@

