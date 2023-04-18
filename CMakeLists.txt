cmake_minimum_required(VERSION 3.16)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release;Profile" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-0g -ggdb" CACHE STRING "")
set(CMAKE_C_FLAGS_DEBUG "-0g -ggdb" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native" CACHE STRING "")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native" CACHE STRING "")
set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g -fno-omit-frame-pointer" CACHE STRING "")
set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_RELEASE} -g -fno-omit-frame-pointer" CACHE STRING "")

project(
  multiqueue_experiments
  VERSION 3.0
  LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(utils)

disable_intree_builds()

# Read the system's L1 cache-line size and page size
execute_process(
  COMMAND getconf LEVEL1_DCACHE_LINESIZE
  OUTPUT_VARIABLE DEFAULT_L1_CACHE_LINESIZE
  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(
  COMMAND getconf PAGESIZE
  OUTPUT_VARIABLE DEFAULT_PAGESIZE
  OUTPUT_STRIP_TRAILING_WHITESPACE)

find_package(PAPI)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)
find_package(TBB)

option("MULTIQUEUE_EXP_BUILD_TESTS" "Build unit tests for experiment tools" OFF)
if(PAPI_FOUND)
  option("MULTIQUEUE_EXP_PERFORMANCE_COUNTERS" "Enable performance counters in benchmarks" ON)
else()
  message(STATUS "PAPI library not found, performance counters disabled")
endif()
option("MULTIQUEUE_EXP_PGO_GENERATE"
       "Enable instrumentation for profile-guided optimizations" OFF)
option("MULTIQUEUE_EXP_PGO_USE" "Enable profile-guided optimizations" OFF)
option("MULTIQUEUE_EXP_COVERAGE" "Enable flags for coverage tests" OFF)
option("MULTIQUEUE_EXP_ASAN"
       "Instrucment benchmarks with the address sanitizer" OFF)
option("MULTIQUEUE_EXP_UBSAN"
       "Instrument benchmarks with the undefined behaviour sanitizer" OFF)
option("MULTIQUEUE_EXP_TSAN" "Instrument benchmarks with the thread sanitizer"
       OFF)
set(L1_CACHE_LINESIZE
    ${DEFAULT_L1_CACHE_LINESIZE}
    CACHE STRING "Specify the assumed L1 cache linesize (bytes)")
set(PAGESIZE
    ${DEFAULT_PAGESIZE}
    CACHE STRING "Specify the assumed pagesize (bytes)")

if(NOT "${L1_CACHE_LINESIZE}" MATCHES "^[0-9]+")
  message(SEND_ERROR "Invalid L1 cache linesize")
endif()
if(NOT "${PAGESIZE}" MATCHES "^[0-9]+")
  message(SEND_ERROR "Invalid pagesize")
endif()

if(MULTIQUEUE_EXP_ASAN AND MULTIQUEUE_EXP_TSAN)
  message(
    SEND_ERROR "Cannot enable both address sanitizer and thread sanitizer")
endif()

string(
  CONCAT WARNING_FLAGS
         "-Wall;"
         "-Wextra;"
         "-Wcast-align;"
         "-Wconversion;"
         "-Wnoexcept;"
         "-Wnon-virtual-dtor;"
         "-Wnull-dereference;"
         "-Woverloaded-virtual;"
         "-Wpedantic;"
         "-Wsuggest-attribute=const;"
         "-Wsuggest-attribute=noreturn;"
         "-Wshadow;"
         "-Wsign-conversion;"
         "-Wsuggest-override;"
         "-Wundef;"
         "-Wunreachable-code;"
         "-Wunused;")

set(EXTERNAL_DIR "${CMAKE_CURRENT_LIST_DIR}/third_party")

add_subdirectory(multiqueue)
add_subdirectory(src)
add_subdirectory(tools)
add_subdirectory(wrapper)

if(MULTIQUEUE_EXP_BUILD_TESTS AND TARGET Catch2::Catch2WithMain)
  include(CTest)
  add_subdirectory(tests)
endif()
