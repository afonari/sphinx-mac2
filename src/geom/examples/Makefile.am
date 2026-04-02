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

NAME = "geom_examples"
MODE = "withBin"

SX_DEPLIBS  = $(SXGEOM) $(SXDIRAC) $(SXIPC) $(SXIO) $(SXFS) $(SXIPC) $(SXMATH)
SX_DEPLIBS += $(SXUTIL) $(SXNONSTD)
NUM_DEPLIBS = $(BLAS)

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist
include $(top_srcdir)/system/make.am.apps

$(top_builddir)/geom/libsxgeom.la:
	$(MAKE) -C .. libsxgeom.la
