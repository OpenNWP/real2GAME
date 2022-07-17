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
src/coordinator.c
src/coordinator.f90)
add_definitions(-Wall)
find_package(OpenMP)
SET(CMAKE_C_FLAGS "${OpenMP_C_FLAGS} -O2 -Wall")
SET(CMAKE_Fortran_FLAGS "-O2 ${OpenMP_Fortran_FLAGS} -Wall -Wno-c-binding-type -I/usr/include -L/usr/lib/x86_64-linux-gnu -lnetcdff")
target_link_libraries(interpolation_creator m eccodes netcdf)
enable_testing()
