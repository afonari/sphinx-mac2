# ----------------------------------------------------------------------------
#
#           The general purpose cross platform C/C++ framework
#
#                       S x A c c e l e r a t e
#
#           Home:       https://www.sxlib.de
#           License:    Apache 2
#           Authors:    see src/AUTHORS
#
# ----------------------------------------------------------------------------

NAME      = "sxio"
MODE      = "libOnly"
VERSION   = "1:0:0"
NAME_HASH = 2090740616  # sxloghash



SX_DEPLIBS  = $(SXIPC) $(SXFS) $(SXUTIL) $(SXNONSTD)
NUM_DEPLIBS = $(REGEX)

#SUBDIRS =
#if COND_BUILD_EXAMPLES
#   SUBDIRS += examples
#endif

# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist

AM_CXXFLAGS += @CXX_PEDANTIC@

