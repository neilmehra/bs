cmake_minimum_required(VERSION 3.14)
project(bs LANGUAGES CXX)
enable_testing()

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1) # to help clangd detect sources
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

add_library(bs INTERFACE)
target_include_directories(bs INTERFACE include)

# external libs
include(FetchContent)
set(BENCHMARK_ENABLE_TESTING off)    # to suppress benchmark internal tests
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
)

FetchContent_Declare(googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG main)

FetchContent_MakeAvailable(googletest googlebenchmark)

FetchContent_Declare(
  tbb
  GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
  GIT_TAG v2022.0.0
  GIT_SHALLOW ON
)
FetchContent_MakeAvailable(tbb)

# Link gtest
file(GLOB_RECURSE TEST_SOURCES "tests/*.cpp")
foreach(TEST_SOURCE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_SOURCE})
    target_link_libraries(${TEST_NAME} bs gtest gtest_main TBB::tbb)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()

# Link googlebenchmark
file(GLOB_RECURSE BENCH_SOURCES "benches/*.cpp")
foreach(BENCH_SOURCE ${BENCH_SOURCES})
    get_filename_component(BENCH_NAME ${BENCH_SOURCE} NAME_WE)
    add_executable(${BENCH_NAME} ${BENCH_SOURCE})
    target_link_libraries(${BENCH_NAME} bs benchmark::benchmark TBB::tbb)
    add_test(NAME ${BENCH_NAME} COMMAND ${BENCH_NAME})
endforeach()
