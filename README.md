# GOTM

A comprehensive description of GOTM including compilation instructions are given at the official [GOTM homepage](http://gotm.net).

This copy of GOTM code is used as a driver to setup simulations to compare different ocean surface boundary layer schemes with Langmuir turbulence. It is built on [GOTM](https://github.com/gotm-model/code) version 5.0.0. New capabilities include the following.

* Stokes drift support. The full profile (grid cell averaged), surface values and the penetration depth of Stokes drift are calculated from wave spectrum (e.g., wave buoy data) or partitioned surface Stokes drift (e.g., WAVEWATCH III output).
* [CVMix](https://github.com/CVMix/CVMix-src) support.
* Langmuir mixing parameterization and Langmuir turbulence enhanced entrainment in KPP via CVMix ([Li et al., 2016](https://doi.org/10.1016%2Fj.ocemod.2015.07.020); [Li and Fox-Kemper, 2017](https://doi.org/10.1175%2FJPO-D-17-0085.1)).
* Langmuir turbulence parameterization in tropical cyclone conditions via CVMix ([Reichl et al., 2016](https://doi.org/10.1175/JPO-D-15-0106.1)).
* ePBL and ePBL-LT.
* OSMOSIS scheme.
* Second moment closure models of Langmuir turbulence ([Harcourt 2013](https://doi.org/10.1175%2FJPO-D-12-0105.1), [2015](https://doi.org/10.1175%2FJPO-D-14-0046.1)).
* UCLA ROMS KPP ([McWilliams et al., 2009](https://doi.org/10.1175%2F2009JPO4130.1)).

# Install

Refer to the [software requirements](http://gotm.net/portfolio/software/) and [instructions for Linux/Mac](http://gotm.net/software/linux/) for a general guide of installing GOTM.

## CMake

[CMake](https://cmake.org) is used in this version of GOMT to configure the code. The options and rules are listed in `${GOTM_ROOT}/src/CMakeLists.txt`. Options other than the default can be set by adding flags when compile GOTM with `cmake`, e.g., use `-DGOTM_USE_CVMix=true` to compile with CVMix.

## NetCDF

[NetCDF](https://www.unidata.ucar.edu/software/netcdf/) is required to compile GOTM by default. CMake uses `nf-config` to determine the appropriate path for `NetCDF_LIBRARIES` and `NetCDF_INCLUDE_DIRS`. One can check if `nf-config` is working correctly by typing
```sh
nf-config --all
```
in the terminal. See `${GOTM_ROOT}/src/cmake/Modules/FindNetCDF.cmake` for more rules CMake uses to find the NetCDF library.

**A note for Mac users**: If NetCDF is installed using [Homebrew](https://brew.sh) (at least for NetCDF v4.5.0), `nf-config` is not working correctly and gives an error `nf-config not yet implemented for cmake builds`. [MacPort](https://www.macports.org) version seems fine. A workaround is to use your own `nf-config` instead of the Homebrew version. An example (NetCDF v4.5.0 installed in `/usr/local` with Fortran compiler `gfortran`) is given in `${GOTM_ROOT}/scripts/nf-config`.

## CVMix

To use [CVMix](https://github.com/CVMix/CVMix-src) in GOTM, CVMix needs to be compiled separately. See [here](https://github.com/CVMix/CVMix-src) for more details on how to compile CVMix. Set the environment variable `CVMIX_ROOT` to the directory where precompiled CVMix is located. CMake uses it to locate the CVMix library. The rules are set in the file `${GOTM_ROOT}/src/cmake/Modules/FindCVMix.cmake`.

To compile GOTM with CVMix, add the flag `-DGOTM_USE_CVMix=true` for `cmake`, i.e.,
```sh
cmake ${srcdir} -DGOTM_USE_CVMix=true
make install
```
where `${srcdir}` is the directory of the GOTM source code.

Also make sure `#define KPP_CVMIX` is set in `${GOTM_ROOT}/include/cppdefs.h`. This allows compiling GOTM without CVMix by `#undef KPP_CVMIX`.
