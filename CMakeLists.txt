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
real2game
interpolator/control.F90
mo_shared.F90)
add_definitions(-Wall)
find_package(OpenMP)
# aggressive optimization
if (FAST)
  SET(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -Ofast -ftree-vectorize -funroll-loops -Wall -Wpedantic -I/usr/include")
# default optimization
else()
  SET(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -O2 -ftree-vectorize -funroll-loops -Wall -Wpedantic -I/usr/include")
endif()
# switching on debugging options
if(DEBUGGING)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbacktrace -fdump-core")
endif()
target_link_libraries(real2game netcdff)
enable_testing()
