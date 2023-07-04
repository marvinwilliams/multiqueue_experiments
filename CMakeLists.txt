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
find_package(TBB)

option("MULTIQUEUE_EXP_BUILD_TESTS" "Build unit tests for experiment tools" OFF)
if(PAPI_FOUND)
  option("MULTIQUEUE_EXP_PERFORMANCE_COUNTERS"
         "Enable performance counters in benchmarks" ON)
else()
  message(STATUS "PAPI library not found, performance counters disabled")
endif()
option("MULTIQUEUE_EXP_PGO" "Use profile-guided optimizations" OFF)
cmake_dependent_option(
  "MULTIQUEUE_EXP_PGO_USE" "Use the generated profile for optimizations" OFF
  "MULTIQUEUE_EXP_PGO" OFF)
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
set(NUM_NUMA_NODES
    1
    CACHE STRING "Specify the number of NUMA nodes")
set(CORES_PER_NUMA_NODE
    256
    CACHE STRING "Specify the number of cores per NUMA node")
option(MULTIQUEUE_EXP_COMPARE_STRICT "Repeat comparison if key is stale" ON)
option(MULTIQUEUE_EXP_COUNT_STATS "Count various stats" OFF)
option(MULTIQUEUE_EXP_NO_SCAN_ON_FAILED_POP
       "Don't scan queues for elements if popping failes" OFF)

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
         "-Wshadow;"
         "-Wsign-conversion;"
         "-Wsuggest-override;"
         "-Wundef;"
         "-Wunreachable-code;"
         "-Wunused;")

add_subdirectory(multiqueue)
# add_subdirectory(third_party/smq/Galois-2.2.1)
add_subdirectory(wrapper)
add_subdirectory(benchmarks)
add_subdirectory(tools)

if(MULTIQUEUE_EXP_BUILD_TESTS AND TARGET Catch2::Catch2WithMain)
  include(CTest)
  add_subdirectory(tests)
endif()
