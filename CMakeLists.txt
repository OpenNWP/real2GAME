cmake_minimum_required(VERSION 3.10)
project(ndvar VERSION 0.0.2)
set(CMAKE_C_COMPILER mpicc)
add_executable(
ndvar src/coordinator.c
src/oi.c
src/lin_algebra_dry.c
src/lin_algebra_moist.c)
add_definitions(-Wall)
find_package(OpenMP)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O1")
target_link_libraries(ndvar m netcdf geos95 atmostracers)
install(TARGETS ndvar DESTINATION bin)
enable_testing()
