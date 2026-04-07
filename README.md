All instructions are for after `source build_env` in the mmshare.

# Mac

Arch for rosseta (not sure)
```
arch -x86_64 ./configure --without-tools --disable-debug --enable-mkl --enable-mklfft --with-mklpath="$INTEL_MKL"/include:"$INTEL_MKL"/lib --disable-openmp --disable-mpi --disable-numlibschecks --prefix=$PWD/install_dir CXXFLAGS="-arch x86_64" CFLAGS="-arch x86_64" LDFLAGS="-arch x86_64"
```

```
make AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=: install && make AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=: -C src/tools sxdefectalign && cp src/tools/.libs/sxdefectalign install_dir/bin/
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
./configure ac_cv_lib_gfortran_rand=no ac_cv_lib_gfortran__gfortran_copy_string=no --without-tools --disable-debug --enable-mkl --enable-mklfft --with-mklpath=$INTEL_MKL --disable-openmp --disable-mpi --disable-numlibschecks --prefix=$PWD/install_dir
```

```
make -j8 AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=: install && make AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=: -C src/tools sxdefectalign && cp src/tools/.libs/sxdefectalign install_dir/bin/
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

# Windows:

```
./configure --without-tools --disable-debug --enable-mkl --enable-mklfft --with-mklpath="$INTEL_MKL"/include:"$INTEL_MKL"/lib/intel64_win --disable-openmp --disable-mpi --disable-numlibschecks --prefix=$PWD/install_dir LIBS="-ldbghelp"
```

```
make -j8 AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=: install && make AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=: -C src/tools sxdefectalign.exe && cp src/tools/.libs/sxdefectalign.exe install_dir/bin/
```

```
$ ldd  /c/source/sphinx-mac2/install_dir/bin/sxdefectalign.exe
        ntdll.dll => /c/Windows/SYSTEM32/ntdll.dll (0x7ff9edac0000)
        KERNEL32.DLL => /c/Windows/System32/KERNEL32.DLL (0x7ff9ec4c0000)
        KERNELBASE.dll => /c/Windows/System32/KERNELBASE.dll (0x7ff9eb260000)
        msvcrt.dll => /c/Windows/System32/msvcrt.dll (0x7ff9ec410000)
        libsxdft-1.dll => /c/source/sphinx-mac2/install_dir/bin/libsxdft-1.dll (0x7ff9b71d0000)
        libsxaddutil-1.dll => /c/source/sphinx-mac2/install_dir/bin/libsxaddutil-1.dll (0x7ff9cfb60000)
        libsxdirac-1.dll => /c/source/sphinx-mac2/install_dir/bin/libsxdirac-1.dll (0x7ff9cfa80000)
        libsxio2-1.dll => /c/source/sphinx-mac2/install_dir/bin/libsxio2-1.dll (0x7ff9cf910000)
        libsxgeom-1.dll => /c/source/sphinx-mac2/install_dir/bin/libsxgeom-1.dll (0x7ff9c1e50000)
        libgcc_s_seh-1.dll => /mingw64/bin/libgcc_s_seh-1.dll (0x7ff9e3130000)
        libsxmath-1.dll => /c/source/sphinx-mac2/install_dir/bin/libsxmath-1.dll (0x7ff9c1d50000)
        libsxutil-1.dll => /c/source/sphinx-mac2/install_dir/bin/libsxutil-1.dll (0x7ff9c1c70000)
        libsxlgpl-1.dll => /c/source/sphinx-mac2/install_dir/bin/libsxlgpl-1.dll (0x7ff9beb30000)
        libsxfs-1.dll => /c/source/sphinx-mac2/install_dir/bin/libsxfs-1.dll (0x7ff9cf880000)
        libsxio-1.dll => /c/source/sphinx-mac2/install_dir/bin/libsxio-1.dll (0x7ff9d9160000)
        libsxmpi-1.dll => /c/source/sphinx-mac2/install_dir/bin/libsxmpi-1.dll (0x7ff9e28e0000)
        libstdc++-6.dll => /mingw64/bin/libstdc++-6.dll (0x7ff9be780000)
        libsxstruct-1.dll => /c/source/sphinx-mac2/install_dir/bin/libsxstruct-1.dll (0x7ff9be5a0000)
        libstdc++-6.dll => /mingw64/bin/libstdc++-6.dll (0x22dcfc50000)
        libsxnonstd-1.dll => /c/source/sphinx-mac2/install_dir/bin/libsxnonstd-1.dll (0x7ff9d32a0000)
        libwinpthread-1.dll => /mingw64/bin/libwinpthread-1.dll (0x7ff9d2ce0000)
        dbghelp.dll => /c/Windows/SYSTEM32/dbghelp.dll (0x7ff9e5e10000)
        ucrtbase.dll => /c/Windows/System32/ucrtbase.dll (0x7ff9eb030000)
        mkl_rt.dll => /c/builds/2026-2/schrodinger_buildenv_packages/.pixi/envs/schrodinger/mkl-2017.1.040/redist/intel64_win/mkl/mkl_rt.dll (0x7ff9b5160000)
        libsxipc-1.dll => /c/source/sphinx-mac2/install_dir/bin/libsxipc-1.dll (0x7ff9cd4d0000)
        ADVAPI32.dll => /c/Windows/System32/ADVAPI32.dll (0x7ff9ed9c0000)
        sechost.dll => /c/Windows/System32/sechost.dll (0x7ff9ec220000)
        bcrypt.dll => /c/Windows/System32/bcrypt.dll (0x7ff9eb710000)
        libsxclassic-1.dll => /c/source/sphinx-mac2/install_dir/bin/libsxclassic-1.dll (0x7ff9c76c0000)
        RPCRT4.dll => /c/Windows/System32/RPCRT4.dll (0x7ff9ec2e0000)
        USER32.dll => /c/Windows/System32/USER32.dll (0x7ff9eba30000)
        win32u.dll => /c/Windows/System32/win32u.dll (0x7ff9eb740000)
        dbgcore.DLL => /c/Windows/SYSTEM32/dbgcore.DLL (0x7ff9e5bb0000)
        GDI32.dll => /c/Windows/System32/GDI32.dll (0x7ff9ec120000)
        gdi32full.dll => /c/Windows/System32/gdi32full.dll (0x7ff9eb140000)
        msvcp_win.dll => /c/Windows/System32/msvcp_win.dll (0x7ff9eb950000)
```

# Notes

* `AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=:` needed to disable calls to autoconf when make is called
* architecture is set to sse4.2 (might not matter for Mac)

# Test

See STU 34409. Example:

```
sxdefectalign --charge 3 --eps 12.9 --vdef v_elec_relaxed_vga_222_chgm3.cub --vref v_elec_pristine_222.cub --ecut 60 --pos 0.5,0.5,0.5 --relative --average 5.2 --qe
Reading mesh+cell...
mesh=135 x 135 x 135
Reading atoms...
New species Z=31
New species Z=33
Reading potential...
cell defect = [a1={21.2161,0,0},a2={0,21.2161,0},a3={0,0,21.2161}]
Reading mesh+cell...
mesh=135 x 135 x 135
Reading atoms...
New species Z=31
New species Z=33
Reading potential...
cell bulk = [a1={21.2161,0,0},a2={0,21.2161,0},a3={0,0,21.2161}]
Excess Electrons = 3 located at {10.608,10.608,10.608}
Atomic structure specified.
ng=74892
V average: -0.0268551 eV
Averaging (33.0881 points)
Averaging (33.0881 points)
Averaging (33.0881 points)
vAlign=0 eV
+-----------------------------------------------------------------------------
=== Intermediate results (unscreened) ===
Isolated energy       : 3.59048
Periodic energy       : 2.98868
Difference (Hartree)  : -0.601801
Difference (eV)       : -16.3758
+-----------------------------------------------------------------------------
Calculation performed with epsilon = 12.9
+-----------------------------------------------------------------------------
Defect correction (eV): 1.26944 (incl. screening & alignment)
```
