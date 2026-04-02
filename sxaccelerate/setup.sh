#!/bin/sh

OS=`uname`

case $OS in
   MINGW*)
      # SxAccelerate is typically nested. Support standalone building
      test -L sxaccelerate || ln -sf . sxaccelerate
      ;;
   *)
      autoreconf -vif
      ;;
esac

if test $? -ne 0; then
   echo "ERROR: Could not generate configure environment."
   exit 1
fi
