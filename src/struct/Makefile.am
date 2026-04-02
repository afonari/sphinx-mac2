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

NAME = "sxstruct"
MODE = "libOnly"
VERSION = "1:0:0"

SX_DEPLIBS  = $(SXDFT) $(SXLGPL) $(SXCLASSIC) $(SXDIRAC) $(SXGEOM) $(SXIO2)
SX_DEPLIBS += $(SXIO) $(SXIPC) $(SXFS) $(SXMATH)
SX_DEPLIBS += $(SXUTIL) $(SXNONSTD)
SX_DEPLIBS += $(SXMPI)
NUM_DEPLIBS =

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist

AM_CPPFLAGS += @CXX_PEDANTIC@
