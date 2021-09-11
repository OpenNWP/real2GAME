cmake_minimum_required(VERSION 3.10)
project(game-da VERSION 0.0.2)
add_executable(
game-da src/coordinator.c
src/oi.c
src/lin_algebra_dry.c
src/lin_algebra_moist.c)
add_definitions(-Wall)
find_package(OpenMP)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
target_link_libraries(game-da m netcdf geos95 atmostracers)
install(TARGETS game-da DESTINATION bin)
enable_testing()
