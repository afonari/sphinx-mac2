All instructions are for after `source build_env` in the mmshare.

# Mac

Arch for rosseta (not sure)
```
arch -x86_64 ./configure --without-tools --disable-debug --enable-mkl --enable-mklfft --with-mklpath="$INTEL_MKL"/include:"$INTEL_MKL"/lib --disable-openmp --disable-mpi --disable-numlibschecks --prefix=$PWD/install_dir CXXFLAGS="-arch x86_64" CFLAGS="-arch x86_64" LDFLAGS="-arch x86_64"
```

```
make -j8 && make install
make -C src/tools sxdefectalign
cp src/tools/.libs/sxdefectalign install_dir/bin/
```

## To list required libs to copy
```
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
```

# Linux

```
unset PERL5LIB
./configure ac_cv_lib_gfortran_rand=no ac_cv_lib_gfortran__gfortran_copy_string=no --without-tools --disable-debug --enable-mkl --enable-mklfft --with-mklpath=$INTEL_MKL --disable-openmp --disable-mpi --disable-numlibschecks --prefix=$PWD/install_dir
```

```
make -j8 && make install
make -C src/tools sxdefectalign
cp src/tools/.libs/sxdefectalign install_dir/bin/
```

## To list required libs to copy
```
ldd install_dir/bin/sxdefectalign 
	linux-vdso.so.1 (0x00007ffce4de3000)
	libsxext.so.1 => /home/fonari/QE-builds/sphinx/sphinx-mac2/install_dir/lib/libsxext.so.1 (0x00007fd750de8000)
	libsxaddutil.so.1 => /home/fonari/QE-builds/sphinx/sphinx-mac2/install_dir/lib/libsxaddutil.so.1 (0x00007fd750d64000)
	libsxstruct.so.1 => /home/fonari/QE-builds/sphinx/sphinx-mac2/install_dir/lib/libsxstruct.so.1 (0x00007fd750bb8000)
	libsxdft.so.1 => /home/fonari/QE-builds/sphinx/sphinx-mac2/install_dir/lib/libsxdft.so.1 (0x00007fd7506ca000)
	libsxlgpl.so.1 => /home/fonari/QE-builds/sphinx/sphinx-mac2/install_dir/lib/libsxlgpl.so.1 (0x00007fd750609000)
	libsxexx.so.1 => /home/fonari/QE-builds/sphinx/sphinx-mac2/install_dir/lib/libsxexx.so.1 (0x00007fd75055e000)
	libsxclassic.so.1 => /home/fonari/QE-builds/sphinx/sphinx-mac2/install_dir/lib/libsxclassic.so.1 (0x00007fd750516000)
	libsxdirac.so.1 => /home/fonari/QE-builds/sphinx/sphinx-mac2/install_dir/lib/libsxdirac.so.1 (0x00007fd750464000)
	libsxgeom.so.1 => /home/fonari/QE-builds/sphinx/sphinx-mac2/install_dir/lib/libsxgeom.so.1 (0x00007fd7503af000)
	libsxio2.so.1 => /home/fonari/QE-builds/sphinx/sphinx-mac2/install_dir/lib/libsxio2.so.1 (0x00007fd75035c000)
	libsxio.so.1 => /home/fonari/QE-builds/sphinx/sphinx-mac2/install_dir/lib/libsxio.so.1 (0x00007fd750341000)
	libsxipc.so.1 => /home/fonari/QE-builds/sphinx/sphinx-mac2/install_dir/lib/libsxipc.so.1 (0x00007fd750303000)
	libsxfs.so.1 => /home/fonari/QE-builds/sphinx/sphinx-mac2/install_dir/lib/libsxfs.so.1 (0x00007fd7502b4000)
	libsxmath.so.1 => /home/fonari/QE-builds/sphinx/sphinx-mac2/install_dir/lib/libsxmath.so.1 (0x00007fd7501d0000)
	libsxutil.so.1 => /home/fonari/QE-builds/sphinx/sphinx-mac2/install_dir/lib/libsxutil.so.1 (0x00007fd75014c000)
	libsxmpi.so.1 => /home/fonari/QE-builds/sphinx/sphinx-mac2/install_dir/lib/libsxmpi.so.1 (0x00007fd750132000)
	libmkl_intel_lp64.so => /scr/fonari/2026-2/schrodinger_buildenv_packages/.pixi/envs/schrodinger/mkl-2017.1.132/lib/intel64/libmkl_intel_lp64.so (0x00007fd74f745000)
	libmkl_core.so => /scr/fonari/2026-2/schrodinger_buildenv_packages/.pixi/envs/schrodinger/mkl-2017.1.132/lib/intel64/libmkl_core.so (0x00007fd74dc9e000)
	libmkl_sequential.so => /scr/fonari/2026-2/schrodinger_buildenv_packages/.pixi/envs/schrodinger/mkl-2017.1.132/lib/intel64/libmkl_sequential.so (0x00007fd74cf2f000)
	libpcre2-8.so.0 => /lib64/libpcre2-8.so.0 (0x00007fd74ccab000)
	libdl.so.2 => /lib64/libdl.so.2 (0x00007fd74caa7000)
	libpthread.so.0 => /lib64/libpthread.so.0 (0x00007fd74c887000)
	libstdc++.so.6 => /scr/fonari/2026-2/schrodinger_buildenv_packages/.pixi/envs/schrodinger/gcc-11.3.0/lib64/libstdc++.so.6 (0x00007fd74c66a000)
	libm.so.6 => /lib64/libm.so.6 (0x00007fd74c2e8000)
	libgcc_s.so.1 => /scr/fonari/2026-2/schrodinger_buildenv_packages/.pixi/envs/schrodinger/gcc-11.3.0/lib64/libgcc_s.so.1 (0x00007fd74c2ce000)
	libc.so.6 => /lib64/libc.so.6 (0x00007fd74bef8000)
	/lib64/ld-linux-x86-64.so.2 (0x00007fd750d21000)
```

