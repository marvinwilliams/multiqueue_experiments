cmake_minimum_required(VERSION 3.16)
message(STATUS "Hostname" ${HOSTNAME})
message(STATUS "Host" ${HOST})
message(STATUS "Hostname" $ENV{HOSTNAME})
message(STATUS "Host" $ENV{HOST})

include(cmake/utils.cmake)
discourage_intree_builds()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Read the system's L1 cache-line size and page size
execute_process(COMMAND getconf LEVEL1_DCACHE_LINESIZE OUTPUT_VARIABLE DEFAULT_L1_CACHE_LINESIZE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND getconf PAGESIZE OUTPUT_VARIABLE DEFAULT_PAGESIZE OUTPUT_STRIP_TRAILING_WHITESPACE)

option("MULTIQUEUE_EXP_BUILD_TESTS" "Build unit tests for experiment tools" OFF)
option("MULTIQUEUE_EXP_PC" "Enable performance counters in benchmarks" OFF)
option("MULTIQUEUE_EXP_PRESERVE_FP" "Do not omit the frame pointer" OFF)
option("MULTIQUEUE_EXP_ASAN" "Instrucment benchmarks with the address sanitizer" OFF)
option("MULTIQUEUE_EXP_UBSAN" "Instrument benchmarks with the undefined behaviour sanitizer" OFF)
option("MULTIQUEUE_EXP_TSAN" "Instrument benchmarks with the thread sanitizer" OFF)
option("MULTIQUEUE_EXP_EXACT_TERMINATION" "Use exact termination for sssp and knapsack benchmarks" ON)
option("MULTIQUEUE_EXP_COUNT_STATS" "Count and print stats for the sssp and knapsack benchmarks" ON)
set(L1_CACHE_LINESIZE ${DEFAULT_L1_CACHE_LINESIZE} CACHE STRING "Specify the assumed L1 cache linesize (bytes)")
set(PAGESIZE ${DEFAULT_PAGESIZE} CACHE STRING "Specify the assumed pagesize (bytes)")

if(NOT "${L1_CACHE_LINESIZE}" MATCHES "^[0-9]+")
  message(SEND_ERROR "Invalid L1 cache linesize")
endif()

if(NOT "${PAGESIZE}" MATCHES "^[0-9]+")
  message(SEND_ERROR "Invalid pagesize")
endif()

if(MULTIQUEUE_EXP_ASAN AND MULTIQUEUE_EXP_TSAN)
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
target_include_directories(benchmark_base INTERFACE "${CMAKE_SOURCE_DIR}/external" "utils" "")
target_include_directories(benchmark_base INTERFACE "${CMAKE_SOURCE_DIR}/external" "utils")
target_compile_options(benchmark_base INTERFACE $<$<CONFIG:Release>:-march=native>)
if(CMAKE_COMPILER_IS_GNUCXX)
  target_link_options(benchmark_base INTERFACE $<$<CONFIG:Release>:-flto>)
endif()
target_compile_options(benchmark_base INTERFACE ${WARNING_FLAGS})
target_compile_definitions(benchmark_base INTERFACE L1_CACHE_LINESIZE=${L1_CACHE_LINESIZE} PAGESIZE=${PAGESIZE})

if(MULTIQUEUE_EXP_PC)
  target_compile_definitions(benchmark_base INTERFACE USE_PAPI)
endif()
if(MULTIQUEUE_EXP_PRESERVE_FP)
  target_compile_options(benchmark_base INTERFACE -fno-omit-frame-pointer)
endif()
if(MULTIQUEUE_EXP_ASAN)
  target_compile_options(benchmark_base INTERFACE -fno-omit-frame-pointer)
  target_compile_options(benchmark_base INTERFACE -fsanitize=address)
  target_link_options(benchmark_base INTERFACE -fsanitize=address)
endif()
if(MULTIQUEUE_EXP_UBSAN)
  target_compile_options(benchmark_base INTERFACE -fno-omit-frame-pointer)
  target_compile_options(benchmark_base INTERFACE -fsanitize=undefined)
  target_link_options(benchmark_base INTERFACE -fsanitize=undefined)
endif()
if(MULTIQUEUE_EXP_TSAN)
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
if(MULTIQUEUE_EXP_EXACT_TERMINATION)
  target_compile_definitions(sssp INTERFACE EXACT_TERMINATION)
endif()
if(MULTIQUEUE_EXP_COUNT_STATS)
  target_compile_definitions(sssp INTERFACE COUNT_STATS)
endif()

add_library(knapsack_packed INTERFACE)
target_sources(knapsack_packed INTERFACE knapsack.cpp)
target_compile_features(knapsack_packed INTERFACE cxx_std_17)
target_link_libraries(knapsack_packed INTERFACE thread_coordination benchmark_base)
target_compile_definitions(knapsack_packed INTERFACE PACKED_VALUE)
if(MULTIQUEUE_EXP_EXACT_TERMINATION)
  target_compile_definitions(knapsack_packed INTERFACE EXACT_TERMINATION)
endif()
if(MULTIQUEUE_EXP_COUNT_STATS)
  target_compile_definitions(knapsack_packed INTERFACE COUNT_STATS)
endif()

add_library(knapsack_heap INTERFACE)
target_sources(knapsack_heap INTERFACE knapsack.cpp)
target_compile_features(knapsack_heap INTERFACE cxx_std_17)
target_link_libraries(knapsack_heap INTERFACE thread_coordination benchmark_base)
target_compile_definitions(knapsack_heap INTERFACE HEAP_VALUE)
if(MULTIQUEUE_EXP_EXACT_TERMINATION)
  target_compile_definitions(knapsack_heap INTERFACE EXACT_TERMINATION)
endif()
if(MULTIQUEUE_EXP_COUNT_STATS)
  target_compile_definitions(knapsack_heap INTERFACE COUNT_STATS)
endif()

add_library(knapsack_explicit INTERFACE)
target_sources(knapsack_explicit INTERFACE knapsack.cpp)
target_compile_features(knapsack_explicit INTERFACE cxx_std_17)
target_link_libraries(knapsack_explicit INTERFACE thread_coordination benchmark_base)
target_compile_definitions(knapsack_explicit INTERFACE EXPLICIT_VALUE)
if(MULTIQUEUE_EXP_EXACT_TERMINATION)
  target_compile_definitions(knapsack_explicit INTERFACE EXACT_TERMINATION)
endif()
if(MULTIQUEUE_EXP_COUNT_STATS)
  target_compile_definitions(knapsack_explicit INTERFACE COUNT_STATS)
endif()

add_custom_target(knapsack_all)

function(add_throughput targets)
  foreach(target ${targets})
    add_executable(throughput_${target})
    target_link_libraries(throughput_${target} PRIVATE ${target} throughput)
  endforeach()
endfunction()

function(add_quality targets)
  foreach(target ${targets})
    add_executable(quality_${target})
    target_link_libraries(quality_${target} PRIVATE ${target} quality)
    if(BUILD_TESTING)
      add_test(NAME quality_${target}_small COMMAND /bin/bash -c "$<TARGET_FILE:quality_${target}> -p 10000 -n 10000 -j 1 | $<TARGET_FILE:evaluate_quality> -v")
      add_test(NAME quality_${target}_large COMMAND /bin/bash -c "$<TARGET_FILE:quality_${target}> -j 8 | $<TARGET_FILE:evaluate_quality> -v")
    endif()
  endforeach()
endfunction()

function(add_sssp targets)
  foreach(target ${targets})
    add_executable(sssp_${target})
    target_link_libraries(sssp_${target} PRIVATE ${target} sssp)
    if(BUILD_TESTING)
      add_test(NAME sssp_${target}_small COMMAND /bin/bash -c "$<TARGET_FILE:sssp_${target}> -j 1 -f data/NY_graph.gr -c data/NY_solution.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
      add_test(NAME sssp_${target}_large COMMAND /bin/bash -c "$<TARGET_FILE:sssp_${target}> -j 8 -f data/CAL_graph.gr -c data/CAL_solution.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endif()
  endforeach()
endfunction()

function(add_knapsack_packed targets)
  foreach(target ${targets})
    add_executable(knapsack_packed_${target})
    target_link_libraries(knapsack_packed_${target} PRIVATE ${target} knapsack_packed)
    add_dependencies(knapsack_all knapsack_packed_${target})
    if(BUILD_TESTING)
      add_test(NAME knapsack_packed_${target}_small COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_packed_${target}> -j 1 -f  data/knapsack_small.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
      add_test(NAME knapsack_packed_${target}_large COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_packed_${target}> -j 8 -f data/knapsack_large.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endif()
  endforeach()
endfunction()

function(add_knapsack_heap targets)
  foreach(target ${targets})
    add_executable(knapsack_heap_${target})
    target_link_libraries(knapsack_heap_${target} PRIVATE ${target} knapsack_heap tbbmalloc)
    add_dependencies(knapsack_all knapsack_heap_${target})
    if(BUILD_TESTING)
      add_test(NAME knapsack_heap_${target}_small COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_external_${target}> -j 1 -f  data/knapsack_small.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
      add_test(NAME knapsack_heap_${target}_large COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_external_${target}> -j 8 -f data/knapsack_large.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endif()
  endforeach()
endfunction()

function(add_knapsack_explicit targets)
  foreach(target ${targets})
    add_executable(knapsack_explicit_${target})
    target_link_libraries(knapsack_explicit_${target} PRIVATE ${target} knapsack_explicit)
    add_dependencies(knapsack_all knapsack_explicit_${target})
    if(BUILD_TESTING)
      add_test(NAME knapsack_explicit_${target}_small COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_explicit_${target}> -j 1 -f  data/knapsack_small.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
      add_test(NAME knapsack_explicit_${target}_large COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_explicit_${target}> -j 8 -f data/knapsack_large.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endif()
  endforeach()
endfunction()

add_subdirectory(multiqueue)
add_subdirectory(multififo)
add_subdirectory(mt-kahypar)

add_library(mq_random INTERFACE)
target_compile_definitions(mq_random INTERFACE PQ_MQ MULTIQUEUE_DEFAULT_SELECTION_STRATEGY=selection_strategy::Random)
target_link_libraries(mq_random INTERFACE multiqueue::multiqueue)

add_library(mq_sticky INTERFACE)
target_compile_definitions(mq_sticky INTERFACE PQ_MQ MQ_HAS_STICKINESS MULTIQUEUE_DEFAULT_SELECTION_STRATEGY=selection_strategy::Sticky)
target_link_libraries(mq_sticky INTERFACE multiqueue::multiqueue)

add_library(mq_swapping INTERFACE)
target_compile_definitions(mq_swapping INTERFACE PQ_MQ MQ_HAS_STICKINESS MULTIQUEUE_DEFAULT_SELECTION_STRATEGY=selection_strategy::Swapping)
target_link_libraries(mq_swapping INTERFACE multiqueue::multiqueue)

add_library(mq_perm INTERFACE)
target_compile_definitions(mq_perm INTERFACE PQ_MQ MQ_HAS_STICKINESS MULTIQUEUE_DEFAULT_SELECTION_STRATEGY=selection_strategy::Permuting)
target_link_libraries(mq_perm INTERFACE multiqueue::multiqueue)

add_library(mf_sticky INTERFACE)
target_compile_definitions(mf_sticky INTERFACE PQ_MF MQ_HAS_STICKINESS MULTIFIFO_DEFAULT_SELECTION_STRATEGY=selection_strategy::Sticky)
target_link_libraries(mf_sticky INTERFACE multififo::multififo)

add_library(mf_swapping INTERFACE)
target_compile_definitions(mf_swapping INTERFACE PQ_MF MQ_HAS_STICKINESS MULTIFIFO_DEFAULT_SELECTION_STRATEGY=selection_strategy::Swapping)
target_link_libraries(mf_swapping INTERFACE multififo::multififo)

add_throughput("mq_random;mq_sticky;mq_swapping;mq_perm;capq;klsm256;klsm1024;klsm4096;linden;spraylist;tbb_q;tbb_pq;mf_sticky;mf_swapping")
add_quality("mq_random;mq_sticky;mq_swapping;mq_perm;capq;klsm256;klsm1024;klsm4096;linden;spraylist;tbb_pq")
add_sssp("mq_random;mq_sticky;mq_swapping;mq_perm;capq;klsm256;klsm1024;klsm4096;linden;spraylist;tbb_q;tbb_pq;mf_sticky")
add_knapsack_packed("mq_random;mq_sticky;mq_swapping;mq_perm;capq;klsm256;klsm1024;klsm4096;linden;spraylist;tbb_q;tbb_pq;mf_sticky")
add_knapsack_heap("mq_random;mq_sticky;mq_swapping;mq_perm;capq;klsm256;klsm1024;klsm4096;linden;spraylist;tbb_q;tbb_pq;mf_sticky")
add_knapsack_explicit("mq_random;mq_sticky;mq_swapping;mq_perm;klsm256;klsm1024;klsm4096;tbb_q;tbb_pq;mf_sticky")

if(MULTIQUEUE_EXP_BUILD_TESTS)
    include(CTest)
  if(NOT TARGET Catch2::Catch2WithMain)
    add_subdirectory("multiqueue/external/Catch2" EXCLUDE_FROM_ALL)
  endif()
  add_subdirectory(tests)
endif()
