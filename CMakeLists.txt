cmake_minimum_required(VERSION 3.16)

project(
  multiqueue_experiments
  VERSION 3.0
  LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(utils)
include(CMakeDependentOption)

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

if(NOT PAPI_FOUND)
  message(STATUS "PAPI library not found, performance counters not available")
endif()

option("EXP_BUILD_TESTS" "Build unit tests for experiment tools" OFF)
set(L1_CACHE_LINESIZE
    ${DEFAULT_L1_CACHE_LINESIZE}
    CACHE STRING "The assumed L1 cache linesize (bytes)")
set(PAGESIZE
    ${DEFAULT_PAGESIZE}
    CACHE STRING "The assumed pagesize (bytes)")
set(NUM_NUMA_NODES
    1
    CACHE STRING "The number of NUMA nodes")
set(CORES_PER_NUMA_NODE
    256
    CACHE STRING "The number of cores per NUMA node")
option(EXP_TUNING_TARGETS "Generate tuning targets" OFF)
if(EXP_TUNING_TARGETS)
  set(EXP_TUNING_POP_PQS
      "1;2;4;8"
      CACHE
        STRING
        "Generate tuning targets for the specified list of number of pop pqs")
  set(EXP_TUNING_BUFFER_SIZES
      "0;4;16;64;256;1024"
      CACHE STRING
            "Generate tuning targets for the specified list of buffer sizes")
  set(EXP_TUNING_HEAP_ARITIES
      "2;4;8;16"
      CACHE STRING
            "Generate tuning targets for the specified list of heap arities")
endif()

if(NOT "${L1_CACHE_LINESIZE}" MATCHES "^[0-9]+")
  message(SEND_ERROR "Invalid L1 cache linesize")
endif()
if(NOT "${PAGESIZE}" MATCHES "^[0-9]+")
  message(SEND_ERROR "Invalid pagesize")
endif()

add_subdirectory(external/multiqueue)
add_subdirectory(benchmarks)
add_subdirectory(tools)

if(EXP_BUILD_TESTS AND TARGET Catch2::Catch2WithMain)
  include(CTest)
  add_subdirectory(tests)
endif()
