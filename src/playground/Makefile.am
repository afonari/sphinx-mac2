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

NAME = "sxplay"
MODE = "withApps"
VERSION = "1:0:0"

SX_DEPLIBS  = $(SXEXT) $(SXADDUTIL) $(SXSTRUCT) $(SXDFT) $(SXLGPL) $(SXCLASSIC) $(SXEXX)
SX_DEPLIBS += $(SXDIRAC) $(SXGEOM) $(SXIO2) $(SXIO) $(SXIPC) $(SXFS)
SX_DEPLIBS += $(SXMATH) $(SXUTIL) $(SXNONSTD)
SX_DEPLIBS += $(SXMPI)
NUM_DEPLIBS = $(FFT) $(BLAS) $(NETCDF) $(REGEX) $(SHORTCUT)

# ----------------------------------------------------------------------------
ACLOCAL_AMFLAGS = -I system/m4

include $(top_srcdir)/system/make.am.filelist
include $(top_srcdir)/system/make.am.apps

AM_CPPFLAGS += @CXX_PEDANTIC@

