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
src/formatter.c)
add_definitions(-Wall)
find_package(OpenMP)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O2 -Wall -Wpedantic")
target_link_libraries(formatter m netcdf eccodes)
enable_testing()
