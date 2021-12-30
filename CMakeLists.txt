# This source file is part of GAME-DA, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME-DA

cmake_minimum_required(VERSION 3.10)
set(version_major 0)
set(version_minor 1)
set(version_patch 1)
set(version ${version_major}.${version_minor}.${version_patch})
project(game-da VERSION ${version})
enable_language(Fortran)
add_executable(
game-da src/coordinator.c
src/oi.c
src/lin_algebra_dry.c
src/lin_algebra_moist.c
src/lin_algebra_wind.c
src/lin_algebra.f90)
find_package(OpenMP)
SET(CMAKE_C_FLAGS "${OpenMP_C_FLAGS} -O2 -Wall")
SET(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -O2 -Wall -Wno-c-binding-type -Wno-unused-dummy-argument")
target_link_libraries(game-da m netcdf geos95 atmostracers)
enable_testing()
