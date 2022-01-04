cmake_minimum_required(VERSION 3.10)

project(matrix)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(run src/run.cpp)
add_executable(benchmark src/benchmark.cpp)
