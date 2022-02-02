cmake_minimum_required(VERSION 3.10)

project(num++)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

SET(CMAKE_CXX_FLAGS "-O3")
SET(CMAKE_C_FLAGS "-O3")

add_executable(run src/run.cpp)
add_executable(benchmark src/benchmark.cpp)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip
)
FetchContent_MakeAvailable(googletest)


include(GoogleTest)
enable_testing()
add_executable(matrix_test tests/matrix_test.cpp)
target_link_libraries(matrix_test gtest_main)

gtest_discover_tests(matrix_test)
