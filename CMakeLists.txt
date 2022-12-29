cmake_minimum_required(VERSION 3.10)

project(num++)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

SET(CMAKE_CXX_FLAGS "-Ofast -ffast-math -march=native")
SET(CMAKE_C_FLAGS "-Ofast")

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

foreach(test_name matrix linear linalg)
  add_executable(${test_name}_test tests/${test_name}_test.cpp)
  target_link_libraries(${test_name}_test gtest_main)
  gtest_discover_tests(${test_name}_test)
endforeach()
