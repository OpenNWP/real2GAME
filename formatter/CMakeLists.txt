# This source file is part of real2GAME, which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/real2GAME

cmake_minimum_required(VERSION 3.10)
set(version_major 0)
set(version_minor 1)
set(version_patch 1)
set(version ${version_major}.${version_minor}.${version_patch})
project(real2GAME VERSION ${version})
add_executable(
formatter
src/control.f90)
add_definitions(-Wall)
find_package(OpenMP)
SET(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -O2 -Wall -Wpedantic -Wno-c-binding-type -I/usr/include -L/usr/lib/x86_64-linux-gnu")
target_link_libraries(formatter m netcdff eccodes)
enable_testing()
