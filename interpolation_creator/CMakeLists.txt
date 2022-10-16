# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

cmake_minimum_required(VERSION 3.10)
set(version_major 0)
set(version_minor 1)
set(version_patch 1)
set(version ${version_major}.${version_minor}.${version_patch})
project(real2GAME VERSION ${version})
enable_language(Fortran)
add_executable(
interpolation_creator
src/control.f90
../mo_shared.f90)
add_definitions(-Wall)
find_package(OpenMP)
SET(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -O2 -ftree-vectorize -funroll-loops -Wall -Wpedantic -I/usr/include \
-I/usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15/ -I/usr/lib/aarch64-linux-gnu/fortran/gfortran-mod-15")
# switching on debugging options
if(DEBUGGING)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbacktrace -fdump-core")
endif()
target_link_libraries(interpolation_creator m netcdff eccodes_f90)
enable_testing()
