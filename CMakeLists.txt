cmake_minimum_required(VERSION 3.16)

include(cmake/utils.cmake)
discourage_intree_builds()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Read the system's L1 cache-line size and page size
execute_process(COMMAND getconf LEVEL1_DCACHE_LINESIZE OUTPUT_VARIABLE DEFAULT_L1_CACHE_LINESIZE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND getconf PAGESIZE OUTPUT_VARIABLE DEFAULT_PAGESIZE OUTPUT_STRIP_TRAILING_WHITESPACE)

option("experiments_BUILD_TESTS" "Build unit tests for experiment tools" OFF)
option("experiments_PC" "Enable performance counters in benchmarks" OFF)
option("experiments_PRESERVE_FP" "Do not omit the frame pointer" OFF)
option("experiments_ASAN" "Instrucment benchmarks with the address sanitizer" OFF)
option("experiments_UBSAN" "Instrument benchmarks with the undefined behaviour sanitizer" OFF)
option("experiments_TSAN" "Instrument benchmarks with the thread sanitizer" OFF)
option("experiments_EXACT_TERMINATION" "Use exact termination for sssp and knapsack benchmarks" ON)
option("experiments_COUNT_STATS" "Count and print stats for the sssp and knapsack benchmarks" ON)
set(L1_CACHE_LINESIZE ${DEFAULT_L1_CACHE_LINESIZE} CACHE STRING "Specify the assumed L1 cache linesize (bytes)")
set(PAGESIZE ${DEFAULT_PAGESIZE} CACHE STRING "Specify the assumed pagesize (bytes)")

if(NOT "${L1_CACHE_LINESIZE}" MATCHES "^[0-9]+")
  message(SEND_ERROR "Invalid L1 cache linesize")
endif()

if(NOT "${PAGESIZE}" MATCHES "^[0-9]+")
  message(SEND_ERROR "Invalid pagesize")
endif()

if(experiments_ASAN AND experiments_TSAN)
  message(FATAL_ERROR "Cannot enable both address sanitizer and thread sanitizer")
endif()

project(multiqueue_experiments VERSION 2.1.0 LANGUAGES CXX C)

# Find the pthread library
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)

find_package(TBB)

set(EXTERNAL_DIR "${CMAKE_CURRENT_LIST_DIR}/external")

add_subdirectory(utils)
add_subdirectory(tools)
add_subdirectory(wrapper)

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

add_library(benchmark_base INTERFACE)
target_include_directories(benchmark_base INTERFACE "${CMAKE_SOURCE_DIR}/external" "utils")
target_compile_options(benchmark_base INTERFACE $<$<CONFIG:Release>:-march=native>)
if(CMAKE_COMPILER_IS_GNUCXX)
  target_link_options(benchmark_base INTERFACE $<$<CONFIG:Release>:-flto>)
endif()
target_compile_options(benchmark_base INTERFACE ${WARNING_FLAGS})
target_compile_definitions(benchmark_base INTERFACE L1_CACHE_LINESIZE=${L1_CACHE_LINESIZE} PAGESIZE=${PAGESIZE})

if(experiments_PC)
  target_compile_definitions(benchmark_base INTERFACE USE_PAPI)
endif()
if(experiments_PRESERVE_FP)
  target_compile_options(benchmark_base INTERFACE -fno-omit-frame-pointer)
endif()
if(experiments_ASAN)
  target_compile_options(benchmark_base INTERFACE -fno-omit-frame-pointer)
  target_compile_options(benchmark_base INTERFACE -fsanitize=address)
  target_link_options(benchmark_base INTERFACE -fsanitize=address)
endif()
if(experiments_UBSAN)
  target_compile_options(benchmark_base INTERFACE -fno-omit-frame-pointer)
  target_compile_options(benchmark_base INTERFACE -fsanitize=undefined)
  target_link_options(benchmark_base INTERFACE -fsanitize=undefined)
endif()
if(experiments_TSAN)
  target_compile_options(benchmark_base INTERFACE -fno-omit-frame-pointer)
  target_compile_options(benchmark_base INTERFACE -fsanitize=thread)
  target_link_options(benchmark_base INTERFACE -fsanitize=thread)
endif()

add_library(quality INTERFACE)
target_sources(quality INTERFACE stress_test.cpp)
target_compile_features(quality INTERFACE cxx_std_17)
target_link_libraries(quality INTERFACE thread_coordination benchmark_base)
target_compile_definitions(quality INTERFACE QUALITY_MODE)

add_library(throughput INTERFACE)
target_sources(throughput INTERFACE stress_test.cpp)
target_compile_features(throughput INTERFACE cxx_std_17)
target_link_libraries(throughput INTERFACE thread_coordination benchmark_base)
target_compile_definitions(throughput INTERFACE THROUGHPUT_MODE)

add_library(sssp INTERFACE)
target_sources(sssp INTERFACE shortest_path.cpp)
target_compile_features(sssp INTERFACE cxx_std_17)
target_link_libraries(sssp INTERFACE thread_coordination benchmark_base)
if(experiments_EXACT_TERMINATION)
  target_compile_definitions(sssp INTERFACE EXACT_TERMINATION)
endif()
if(experiments_COUNT_STATS)
  target_compile_definitions(sssp INTERFACE COUNT_STATS)
endif()

add_library(knapsack_packed INTERFACE)
target_sources(knapsack_packed INTERFACE knapsack.cpp)
target_compile_features(knapsack_packed INTERFACE cxx_std_17)
target_link_libraries(knapsack_packed INTERFACE thread_coordination benchmark_base)
target_compile_definitions(knapsack_packed INTERFACE PACKED_VALUE)
if(experiments_EXACT_TERMINATION)
  target_compile_definitions(knapsack_packed INTERFACE EXACT_TERMINATION)
endif()
if(experiments_COUNT_STATS)
  target_compile_definitions(knapsack_packed INTERFACE COUNT_STATS)
endif()

add_library(knapsack_heap INTERFACE)
target_sources(knapsack_heap INTERFACE knapsack.cpp)
target_compile_features(knapsack_heap INTERFACE cxx_std_17)
target_link_libraries(knapsack_heap INTERFACE thread_coordination benchmark_base)
target_compile_definitions(knapsack_heap INTERFACE HEAP_VALUE)
if(experiments_EXACT_TERMINATION)
  target_compile_definitions(knapsack_heap INTERFACE EXACT_TERMINATION)
endif()
if(experiments_COUNT_STATS)
  target_compile_definitions(knapsack_heap INTERFACE COUNT_STATS)
endif()

add_library(knapsack_explicit INTERFACE)
target_sources(knapsack_explicit INTERFACE knapsack.cpp)
target_compile_features(knapsack_explicit INTERFACE cxx_std_17)
target_link_libraries(knapsack_explicit INTERFACE thread_coordination benchmark_base)
target_compile_definitions(knapsack_explicit INTERFACE EXPLICIT_VALUE)
if(experiments_EXACT_TERMINATION)
  target_compile_definitions(knapsack_explicit INTERFACE EXACT_TERMINATION)
endif()
if(experiments_COUNT_STATS)
  target_compile_definitions(knapsack_explicit INTERFACE COUNT_STATS)
endif()

function(add_benchmarks_for name)
  add_executable(quality_${name})
  target_link_libraries(quality_${name} PRIVATE ${name} quality)

  add_executable(throughput_${name})
  target_link_libraries(throughput_${name} PRIVATE ${name} throughput)

  add_executable(sssp_${name})
  target_link_libraries(sssp_${name} PRIVATE ${name} sssp)

  add_executable(knapsack_packed_${name})
  target_link_libraries(knapsack_packed_${name} PRIVATE ${name} knapsack_packed)

  add_executable(knapsack_heap_${name})
  target_link_libraries(knapsack_heap_${name} PRIVATE ${name} knapsack_heap tbbmalloc)

  add_executable(knapsack_explicit_${name})
  target_link_libraries(knapsack_explicit_${name} PRIVATE ${name} knapsack_explicit)

  if(BUILD_TESTING)
    add_test(NAME quality_${name}_small COMMAND /bin/bash -c "$<TARGET_FILE:quality_${name}> -p 10000 -n 10000 -j 1 | $<TARGET_FILE:evaluate_quality> -v")
    add_test(NAME quality_${name}_large COMMAND /bin/bash -c "$<TARGET_FILE:quality_${name}> -j 8 | $<TARGET_FILE:evaluate_quality> -v")
    add_test(NAME sssp_${name}_small COMMAND /bin/bash -c "$<TARGET_FILE:sssp_${name}> -j 1 -f data/NY_graph.gr -c data/NY_solution.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME sssp_${name}_large COMMAND /bin/bash -c "$<TARGET_FILE:sssp_${name}> -j 8 -f data/CAL_graph.gr -c data/CAL_solution.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME knapsack_packed_${name}_small COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_packed_${name}> -j 1 -f  data/knapsack_small.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME knapsack_packed_${name}_large COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_packed_${name}> -j 8 -f data/knapsack_large.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME knapsack_heap_${name}_small COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_external_${name}> -j 1 -f  data/knapsack_small.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME knapsack_heap_${name}_large COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_external_${name}> -j 8 -f data/knapsack_large.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME knapsack_explicit_${name}_small COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_explicit_${name}> -j 1 -f  data/knapsack_small.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME knapsack_explicit_${name}_large COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_explicit_${name}> -j 8 -f data/knapsack_large.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
endfunction()

add_subdirectory(multiqueue_simple)

add_library(mq_random INTERFACE)
target_compile_definitions(mq_random INTERFACE PQ_MQ_RANDOM MULTIQUEUE_DEFAULT_SELECTION_STRATEGY=selection_strategy::Random)
target_link_libraries(mq_random INTERFACE multiqueue::multiqueue)

add_library(mq_sticky INTERFACE)
target_compile_definitions(mq_sticky INTERFACE PQ_MQ_STICKY MULTIQUEUE_DEFAULT_SELECTION_STRATEGY=selection_strategy::Sticky)
target_link_libraries(mq_sticky INTERFACE multiqueue::multiqueue)

add_library(mq_swapping INTERFACE)
target_compile_definitions(mq_swapping INTERFACE PQ_MQ_SWAPPING MULTIQUEUE_DEFAULT_SELECTION_STRATEGY=selection_strategy::Swapping)
target_link_libraries(mq_swapping INTERFACE multiqueue::multiqueue)

add_library(mq_perm INTERFACE)
target_compile_definitions(mq_perm INTERFACE PQ_MQ_PERM MULTIQUEUE_DEFAULT_SELECTION_STRATEGY=selection_strategy::Permuting)
target_link_libraries(mq_perm INTERFACE multiqueue::multiqueue)

add_benchmarks_for(mq_random)
add_benchmarks_for(mq_sticky)
add_benchmarks_for(mq_swapping)
add_benchmarks_for(mq_perm)
add_benchmarks_for(capq)
add_benchmarks_for(klsm256)
add_benchmarks_for(klsm1024)
add_benchmarks_for(klsm4096)
add_benchmarks_for(linden)
add_benchmarks_for(spraylist)

if(experiments_BUILD_TESTS)
    include(CTest)
  if(NOT TARGET Catch2::Catch2WithMain)
    add_subdirectory("multiqueue/external/Catch2" EXCLUDE_FROM_ALL)
  endif()
  add_subdirectory(tests)
endif()
