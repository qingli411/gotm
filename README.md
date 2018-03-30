# GOTM

A comprehensive description of GOTM including compilation instructions are given at the official [GOTM homepage](http://gotm.net).

This copy of GOTM code is used to setup simulations to compare different ocean surface boundary layer schemes. It is built on [GOTM](https://github.com/gotm-model/code) version 5.0.0. New capabilities include the following.

* KPP using [CVMix](https://github.com/CVMix/CVMix-src)

* Langmuir mixing parameterization and Langmuir turbulence enhanced entrainment in KPP (Li et al., 2016; Li and Fox-Kemper, 2017), both in GOTM version of KPP and via CVMix. _To be tested._

* OSMOSIS scheme _To be included._

* Second moment closure models of Langmuir turbulence (Harcourt 2013, 2015) _To be included._

* Stokes drift calculated from input wave spectrum from wave buoy, or partitioned surface Stokes drift from wave simulations.

# Install

Refer to the [software requirements](http://gotm.net/portfolio/software/) and [instructions for Linux/Mac](http://gotm.net/software/linux/) for a general guide of installing GOTM.

## CMake

[CMake](https://cmake.org) is used in this version of GOMT to configure the code. The options and rules are listed in `${GOTM_ROOT}/src/CMakeLists.txt`.

## NetCDF

[NetCDF](https://www.unidata.ucar.edu/software/netcdf/) is required to compile GOTM by default. CMake uses `nf-config` to determine the correct path for `NetCDF_LIBRARIES` and `NetCDF_INCLUDE_DIRS`. One can check if `nf-config` is working correctly by typing
```sh
nf-config --all
```
in the terminal. See `${GOTM_ROOT}/src/cmake/Modules/FindNetCDF.cmake` for more rules CMake uses to find the NetCDF library.

A note for Mac users: if NetCDF is installed using [Homebrew](https://brew.sh) (at least for NetCDF v4.5.0), `nf-config` is not working correctly and gives an error `nf-config not yet implemented for cmake builds`. [MacPort](https://www.macports.org) version seems fine. A fix is to use your own `nf-config` instead of the Homebrew version. An example (NetCDF v4.5.0 installed in `/usr/local` with Fortran compiler `gfortran`) is given in `${GOTM_ROOT}/scripts/nf-config`.

## CVMix

To use [CVMix](https://github.com/CVMix/CVMix-src) in GOTM, CVMix needs to be compiled separately. See [here](https://github.com/CVMix/CVMix-src) for more details on how to compile CVMix.

To compile GOTM with CVMix, add the flag `-DGOTM_USE_CVMix=true` for `cmake`, i.e.,
```sh
cmake ${srcdir} -DGOTM_USE_CVMix=true
make install
```
where `${srcdir}` is the directory of the GOTM source code.

The file `${GOTM_ROOT}/src/cmake/Modules/FindCVMix.cmake` sets the rules used by CMake to find the CVMix library. By default it assumes the compiled CVMix is located in either `${HOME}/CVMix-src` or `${HOME}/local/CVMix-src`. This file should be modified accordingly if CVMix is located in a different directory. Here is how it works. It tries to determine the root directory of CVMix (`CVMix_PREFIX`) by looking for `src/cvmix_driver.F90` in the above directories. Then it determines `CVMix_LIBRARIES` by looking for `libcvmix` in `${CVMix_PREFIX}/lib` and `CVMix_INCLUDE_DIRS` by looking for `cvmix_kpp.mod` in `${CVMix_PREFIX}/include`.

Also make sure `#define KPP_CVMIX` is set in `${GOTM_ROOT}/include/cppdefs.h`
