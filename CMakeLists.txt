cmake_minimum_required(VERSION 3.10)

project(matrix)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

SET(CMAKE_CXX_FLAGS "-O3")
SET(CMAKE_C_FLAGS "-O3")

add_executable(run src/run.cpp)
add_executable(benchmark src/benchmark.cpp)
