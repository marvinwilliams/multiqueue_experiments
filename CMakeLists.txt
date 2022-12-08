cmake_minimum_required(VERSION 3.16)

project(multiqueue_experiments VERSION 3.0 LANGUAGES CXX C)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(utils)
disable_intree_builds()

# Read the system's L1 cache-line size and page size
execute_process(COMMAND getconf LEVEL1_DCACHE_LINESIZE OUTPUT_VARIABLE DEFAULT_L1_CACHE_LINESIZE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND getconf PAGESIZE OUTPUT_VARIABLE DEFAULT_PAGESIZE OUTPUT_STRIP_TRAILING_WHITESPACE)

option("MULTIQUEUE_EXP_BUILD_TESTS" "Build unit tests for experiment tools" OFF)
option("MULTIQUEUE_EXP_PC" "Enable performance counters in benchmarks" OFF)
option("MULTIQUEUE_EXP_PROFILING" "Add profiling compiler flags" OFF)
option("MULTIQUEUE_EXP_ASAN" "Instrucment benchmarks with the address sanitizer" OFF)
option("MULTIQUEUE_EXP_UBSAN" "Instrument benchmarks with the undefined behaviour sanitizer" OFF)
option("MULTIQUEUE_EXP_TSAN" "Instrument benchmarks with the thread sanitizer" OFF)
set(L1_CACHE_LINESIZE ${DEFAULT_L1_CACHE_LINESIZE} CACHE STRING "Specify the assumed L1 cache linesize (bytes)")
set(PAGESIZE ${DEFAULT_PAGESIZE} CACHE STRING "Specify the assumed pagesize (bytes)")

if(NOT "${L1_CACHE_LINESIZE}" MATCHES "^[0-9]+")
  message(SEND_ERROR "Invalid L1 cache linesize")
endif()
if(NOT "${PAGESIZE}" MATCHES "^[0-9]+")
  message(SEND_ERROR "Invalid pagesize")
endif()

if(MULTIQUEUE_EXP_ASAN AND MULTIQUEUE_EXP_TSAN)
  message(SEND_ERROR "Cannot enable both address sanitizer and thread sanitizer")
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
  "-Wunused;"
)

set(EXTERNAL_DIR "${CMAKE_CURRENT_LIST_DIR}/third_party")

add_subdirectory(multiqueue)
add_subdirectory(multififo)
add_subdirectory(src)
add_subdirectory(tools)
add_subdirectory(wrapper)

if(MULTIQUEUE_EXP_BUILD_TESTS AND TARGET Catch2::Catch2WithMain)
  include(CTest)
  add_subdirectory(tests)
endif()
