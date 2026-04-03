# Mac

arch -x86_64 ./configure --without-tools --disable-debug --enable-mkl --enable-mklfft --with-mklpath="$INTEL_MKL"/include:"$INTEL_MKL"/lib --disable-openmp --disable-mpi --disable-numlibschecks CXXFLAGS="-arch x86_64" CFLAGS="-arch x86_64" LDFLAGS="-arch x86_64" --prefix=$PWD/install_dir

make -j8 && make install

make -C src/tools sxdefectalign
cp src/tools/.libs/sxdefectalign install_dir/bin/

# To list required libs to copy
otool -L install_dir/bin/sxdefectalign

install_dir/bin/sxdefectalign:
    /Users/fonari/src/sphinx-mac2/install_dir/lib/libsxext.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac2/install_dir/lib/libsxaddutil.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac2/install_dir/lib/libsxexx.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac2/install_dir/lib/libsxstruct.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac2/install_dir/lib/libsxclassic.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac2/install_dir/lib/libsxdft.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac2/install_dir/lib/libsxlgpl.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac2/install_dir/lib/libsxdirac.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac2/install_dir/lib/libsxgeom.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac2/install_dir/lib/libsxio2.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac2/install_dir/lib/libsxmath.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac2/install_dir/lib/libsxmpi.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac2/install_dir/lib/libsxio.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac2/install_dir/lib/libsxipc.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac2/install_dir/lib/libsxfs.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    /Users/fonari/src/sphinx-mac2/install_dir/lib/libsxutil.1.dylib (compatibility version 2.0.0, current version 2.0.0)
    @rpath/libmkl_intel_lp64.dylib (compatibility version 0.0.0, current version 0.0.0)
    @rpath/libmkl_core.dylib (compatibility version 0.0.0, current version 0.0.0)
    @rpath/libmkl_sequential.dylib (compatibility version 0.0.0, current version 0.0.0)
    /usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1356.0.0)
    /System/Library/Frameworks/CoreFoundation.framework/Versions/A/CoreFoundation (compatibility version 150.0.0, current version 4201.0.0)
    /usr/lib/libc++.1.dylib (compatibility version 1.0.0, current version 2000.67.0)
