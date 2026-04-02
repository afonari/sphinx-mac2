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

NAME = "io2_examples"
MODE = "withBin"

SX_DEPLIBS  = $(SXIO2) $(SXIO) $(SXIPC) $(SXMATH) $(SXFS)
SX_DEPLIBS += $(SXUTIL) $(SXNONSTD)
NUM_DEPLIBS = $(NETCDF)

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist
include $(top_srcdir)/system/make.am.apps


$(top_builddir)/io2/libsxio2.la:
	$(MAKE) -C .. libsxio2.la
EXTRA_DIST=parser.std parser.sx
