cmake_minimum_required(VERSION 3.16)

include(cmake/utils.cmake)
discourage_intree_builds()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Read the system's L1 cache-line size and page size
execute_process(COMMAND getconf LEVEL1_DCACHE_LINESIZE OUTPUT_VARIABLE DEFAULT_L1_CACHE_LINESIZE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND getconf PAGESIZE OUTPUT_VARIABLE DEFAULT_PAGESIZE OUTPUT_STRIP_TRAILING_WHITESPACE)

option("MULTIQUEUE_EXP_BUILD_TESTS" "Build unit tests for experiment tools" OFF)
option("MULTIQUEUE_EXP_PC" "Enable performance counters in benchmarks" OFF)
option("MULTIQUEUE_EXP_NULL" "Don't do actual work" OFF)
option("MULTIQUEUE_EXP_PROFILING" "Add profiling compiler flags" OFF)
option("MULTIQUEUE_EXP_ASAN" "Instrucment benchmarks with the address sanitizer" OFF)
option("MULTIQUEUE_EXP_UBSAN" "Instrument benchmarks with the undefined behaviour sanitizer" OFF)
option("MULTIQUEUE_EXP_TSAN" "Instrument benchmarks with the thread sanitizer" OFF)
option("MULTIQUEUE_EXP_COUNT_STATS" "Count and print stats for the sssp and knapsack benchmarks" ON)
option(MULTIQUEUE_EXP_PQ_STD "Use std::priority_queue as priority queue" OFF)
option(MULTIQUEUE_EXP_DISABLE_BUFFERING "Disable buffering" OFF)
set(MULTIQUEUE_EXP_DELETION_BUFFERSIZE "" CACHE STRING "Deletion buffer size")
set(MULTIQUEUE_EXP_INSERTION_BUFFERSIZE "" CACHE STRING "Insertion buffer size")
set(MULTIQUEUE_EXP_HEAP_ARITY "" CACHE STRING "Heap arity")
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

project(multiqueue_experiments VERSION 3.0 LANGUAGES CXX C)

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
target_compile_definitions(benchmark_base INTERFACE L1_CACHE_LINESIZE=${L1_CACHE_LINESIZE})
target_compile_definitions(benchmark_base INTERFACE PAGESIZE=${PAGESIZE})

if(MULTIQUEUE_EXP_PC)
  find_library(PAPI_LIB papi)
  find_path(PAPI_INCLUDE_DIR papi.h)
  target_compile_definitions(benchmark_base INTERFACE USE_PAPI)
  target_link_libraries(benchmark_base INTERFACE ${PAPI_LIB})
  target_include_directories(benchmark_base INTERFACE ${PAPI_INCLUDE_DIR})
endif()
if(MULTIQUEUE_EXP_NULL)
  target_compile_definitions(benchmark_base INTERFACE NULL_WORK)
endif()
if(MULTIQUEUE_EXP_PROFILING)
  target_compile_options(benchmark_base INTERFACE -g -fno-omit-frame-pointer)
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

if(MULTIQUEUE_EXP_PQ_STD)
  target_compile_definitions(benchmark_base INTERFACE MQ_PQ_STD)
endif()
if(MULTIQUEUE_EXP_DISABLE_BUFFERING)
  target_compile_definitions(benchmark_base INTERFACE MQ_DISABLE_BUFFERING)
endif()
if(NOT MULTIQUEUE_EXP_DELETION_BUFFERSIZE STREQUAL "")
  target_compile_definitions(benchmark_base INTERFACE MQ_DELETION_BUFFERSIZE=${MULTIQUEUE_EXP_DELETION_BUFFERSIZE})
endif()
if(NOT MULTIQUEUE_EXP_INSERTION_BUFFERSIZE STREQUAL "")
  target_compile_definitions(benchmark_base INTERFACE MQ_INSERTION_BUFFERSIZE=${MULTIQUEUE_EXP_INSERTION_BUFFERSIZE})
endif()
if(NOT MULTIQUEUE_EXP_HEAP_ARITY STREQUAL "")
  target_compile_definitions(benchmark_base INTERFACE MQ_HEAP_ARITY=${MULTIQUEUE_EXP_HEAP_ARITY})
endif()

add_library(quality INTERFACE)
target_sources(quality INTERFACE quality.cpp)
target_compile_features(quality INTERFACE cxx_std_17)
target_link_libraries(quality INTERFACE thread_coordination benchmark_base)

add_library(throughput INTERFACE)
target_sources(throughput INTERFACE throughput.cpp)
target_compile_features(throughput INTERFACE cxx_std_17)
target_link_libraries(throughput INTERFACE thread_coordination benchmark_base)

add_library(insert_delete INTERFACE)
target_sources(insert_delete INTERFACE insert_delete.cpp)
target_compile_features(insert_delete INTERFACE cxx_std_17)
target_link_libraries(insert_delete INTERFACE thread_coordination benchmark_base)

add_library(sssp INTERFACE)
target_sources(sssp INTERFACE shortest_path.cpp)
target_compile_features(sssp INTERFACE cxx_std_17)
target_link_libraries(sssp INTERFACE thread_coordination benchmark_base)
if(MULTIQUEUE_EXP_COUNT_STATS)
  target_compile_definitions(sssp INTERFACE COUNT_STATS)
endif()

add_library(knapsack_packed INTERFACE)
target_sources(knapsack_packed INTERFACE knapsack.cpp)
target_compile_features(knapsack_packed INTERFACE cxx_std_17)
target_link_libraries(knapsack_packed INTERFACE thread_coordination benchmark_base)
target_compile_definitions(knapsack_packed INTERFACE PACKED_VALUE)
target_compile_definitions(knapsack_packed INTERFACE HINT_MODE)
if(MULTIQUEUE_EXP_COUNT_STATS)
  target_compile_definitions(knapsack_packed INTERFACE COUNT_STATS)
endif()

add_library(knapsack_heap INTERFACE)
target_sources(knapsack_heap INTERFACE knapsack.cpp)
target_compile_features(knapsack_heap INTERFACE cxx_std_17)
target_link_libraries(knapsack_heap INTERFACE thread_coordination benchmark_base)
target_compile_definitions(knapsack_heap INTERFACE HEAP_VALUE)
target_compile_definitions(knapsack_heap INTERFACE HINT_MODE)
if(MULTIQUEUE_EXP_COUNT_STATS)
  target_compile_definitions(knapsack_heap INTERFACE COUNT_STATS)
endif()

add_library(knapsack_explicit INTERFACE)
target_sources(knapsack_explicit INTERFACE knapsack.cpp)
target_compile_features(knapsack_explicit INTERFACE cxx_std_17)
target_link_libraries(knapsack_explicit INTERFACE thread_coordination benchmark_base)
target_compile_definitions(knapsack_explicit INTERFACE EXPLICIT_VALUE)
target_compile_definitions(knapsack_explicit INTERFACE HINT_MODE)
if(MULTIQUEUE_EXP_COUNT_STATS)
  target_compile_definitions(knapsack_explicit INTERFACE COUNT_STATS)
endif()

add_custom_target(knapsack_all)

function(add_throughput targets)
  foreach(target ${targets})
    add_executable(throughput_${target})
    target_link_libraries(throughput_${target} PRIVATE ${target} throughput)

    add_executable(insert_delete_${target})
    target_link_libraries(insert_delete_${target} PRIVATE ${target} insert_delete)
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
      add_test(NAME knapsack_heap_${target}_small COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_heap_${target}> -j 1 -f  data/knapsack_small.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
      add_test(NAME knapsack_heap_${target}_large COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_heap_${target}> -j 8 -f data/knapsack_large.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
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

add_library(mq_none INTERFACE)
target_compile_definitions(mq_none INTERFACE PQ_MQ MQ_STICK_POLICY_NONE)
target_link_libraries(mq_none INTERFACE multiqueue::multiqueue)

add_library(mq_random INTERFACE)
target_compile_definitions(mq_random INTERFACE PQ_MQ MQ_STICK_POLICY_RANDOM)
target_link_libraries(mq_random INTERFACE multiqueue::multiqueue)

add_library(mq_random_strict INTERFACE)
target_compile_definitions(mq_random_strict INTERFACE PQ_MQ MQ_STICK_POLICY_RANDOM_STRICT)
target_link_libraries(mq_random_strict INTERFACE multiqueue::multiqueue)

add_library(mq_swapping INTERFACE)
target_compile_definitions(mq_swapping INTERFACE PQ_MQ MQ_STICK_POLICY_SWAPPING)
target_link_libraries(mq_swapping INTERFACE multiqueue::multiqueue)

add_library(mq_perm INTERFACE)
target_compile_definitions(mq_perm INTERFACE PQ_MQ MQ_STICK_POLICY_PERMUTATION)
target_link_libraries(mq_perm INTERFACE multiqueue::multiqueue)

add_library(mf INTERFACE)
target_compile_definitions(mf INTERFACE PQ_MF)
target_link_libraries(mf INTERFACE multififo::multififo)

add_throughput("mq_none;mq_random;mq_random_strict;mq_swapping;mq_perm;capq;klsm256;klsm1024;klsm4096;linden;spraylist;tbb_q;tbb_pq;mf")
add_quality("mq_none;mq_random;mq_random_strict;mq_swapping;mq_perm;capq;klsm256;klsm1024;klsm4096;linden;spraylist;tbb_pq")
add_sssp("mq_none;mq_random;mq_random_strict;mq_swapping;mq_perm;capq;klsm256;klsm1024;klsm4096;linden;spraylist;tbb_pq;tbb_q;mf")
add_knapsack_packed("mq_none;mq_random;mq_random_strict;mq_swapping;mq_perm;capq;klsm256;klsm1024;klsm4096;linden;spraylist;tbb_pq;tbb_q;mf")
add_knapsack_heap("mq_none;mq_random;mq_random_strict;mq_swapping;mq_perm;capq;klsm256;klsm1024;klsm4096;linden;spraylist;tbb_pq;tbb_q;mf")
add_knapsack_explicit("mq_none;mq_random;mq_random_strict;mq_swapping;mq_perm;klsm256;klsm1024;klsm4096;tbb_pq;tbb_q;mf")

if(MULTIQUEUE_EXP_BUILD_TESTS)
  include(CTest)
  if(NOT TARGET Catch2::Catch2WithMain)
    add_subdirectory("multiqueue/external/Catch2" EXCLUDE_FROM_ALL)
  endif()
  add_subdirectory(tests)
endif()

if(CMAKE_EXPORT_COMPILE_COMMANDS)
  set(CMAKE_CXX_STANDARD 17)
endif()
