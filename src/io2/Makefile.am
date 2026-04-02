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

NAME = "sxio2"
MODE = "libOnly"
VERSION = "1:0:0"

#SUBDIRS =
#if COND_BUILD_EXAMPLES
#   SUBDIRS += examples
#endif


SX_DEPLIBS  = $(SXIO) $(SXIPC) $(SXMATH) $(SXFS) $(SXUTIL) $(SXNONSTD)
SX_DEPLIBS += $(SXMPI)
NUM_DEPLIBS = $(NETCDF)

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist

AM_CPPFLAGS += @CXX_PEDANTIC@

# --- netcdf.h uses "long long" which is not allowed in ISO C++ (-pedantic)
AM_CXXFLAGS += -Wno-long-long

# --- LEX and YACC support
AUX_SRC       = SxParser.lex.cpp
AUX_SRC      += SxParser.tab.cpp SxParser.tab.hpp
AM_YFLAGS     = -pSxParser_ -d
BUILT_SOURCES = $(AUX_SRC)
EXTRA_DIST    = SxParser.lpp SxParser.ypp rapidxml-1.13
CLEANFILES   += $(AUX_SRC) SxParser.output


%.lex.lo: AM_CXXFLAGS += @LEX_CXXFLAGS@ -Wno-unused
%.tab.lo: AM_CXXFLAGS += @LEX_CXXFLAGS@

%.lex.cpp: %.lpp %.tab.cpp
	$(LEX) -PSxParser_ $(LEXFLAGS) -o $@ $<

%.tab.cpp %.tab.hpp: %.ypp
	$(YACC) -pSxParser_ -d $< -o $@
