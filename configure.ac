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


AC_INIT([S/PHI/nX], [m4_normalize(m4_esyscmd([cat ./VERSION]))], [https://sxrepo.mpie.de], [sphinx])

ac_default_prefix=`pwd`

AM_INIT_AUTOMAKE([foreign])

AC_CONFIG_FILES([Makefile])

# set prefix to current directory if unset

SX_TITLE([S/PHI/nX])
sx_accelerate=`pwd`/sxaccelerate

SX_CONFIG_SUBDIR([sxaccelerate],[--with-sxmath])
SX_CONFIG_SUBDIR([src],[--with-sxaccelerate=$sx_accelerate/src])

AC_OUTPUT

