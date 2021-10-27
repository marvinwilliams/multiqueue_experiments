cmake_minimum_required(VERSION 3.16)

project(multiqueue_experiments VERSION 0.1.0 LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

include(utils)
include(CMakeDependentOption)
include(CMakePackageConfigHelpers)
include(CTest)
include(GNUInstallDirs)

prevent_intree_builds()

option("experiments_BUILD_TESTS" "Build unit tests for multiqueue experiments" OFF)
option("experiments_PROFILE" "Enable profiling in the benchmarks" OFF)
option("experiments_FP" "Compile benchmarks with frame" OFF)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_DEBUG_POSTFIX _debug)

# Read the system's L1 cache-line size and page size
execute_process(COMMAND getconf LEVEL1_DCACHE_LINESIZE OUTPUT_VARIABLE L1_CACHE_LINESIZE OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT "${L1_CACHE_LINESIZE}" MATCHES "^[0-9]+")
  message(WARNING "Could not read L1 cache-line size, using default: 64")
  set(L1_CACHE_LINESIZE 64)
endif()

# Generate header containing definitions for L1 cache-line size and page size
configure_file(
  "system_config.hpp.in"
  "system_config.hpp"
)
add_library("experiment_base" INTERFACE)

string(
  CONCAT WARNING_FLAGS
  "-Wall;"
  "-Wextra;"
  "-Wpedantic;"
)

target_compile_options("experiment_base" INTERFACE ${WARNING_FLAGS})

add_subdirectory(multiqueue_simple)

# Find the pthread library
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)

find_package(TBB)

set(EXTERNAL_DIR "${CMAKE_CURRENT_LIST_DIR}/external")

if(experiments_BUILD_TESTS)
  add_subdirectory("${EXTERNAL_DIR}/Catch2" EXCLUDE_FROM_ALL)
endif()
list(APPEND CMAKE_MODULE_PATH "${EXTERNAL_DIR}/Catch2/extras")

function(generate_mq_target name compiler_define)
  add_library(${name} INTERFACE)
  target_compile_definitions(${name} INTERFACE ${compiler_define})
  target_link_libraries(${name} INTERFACE multiqueue::multiqueue)
  set_property(TARGET ${name} PROPERTY pq_name ${name})
endfunction()

generate_mq_target(mq_random PQ_MQ_RANDOM)
generate_mq_target(mq_sticky PQ_MQ_STICKY)

add_subdirectory(utils)
add_subdirectory(tools)
add_subdirectory(wrapper)

add_library(benchmark_base INTERFACE)
target_include_directories(benchmark_base INTERFACE "${CMAKE_SOURCE_DIR}/external")
target_link_libraries(benchmark_base INTERFACE experiment_base)
target_compile_options(benchmark_base INTERFACE $<$<CONFIG:Release>:-march=native>)
if(CMAKE_COMPILER_IS_GNUCXX)
  target_link_options(benchmark_base INTERFACE $<$<CONFIG:Release>:-flto>)
endif()

add_library(quality INTERFACE)
target_link_libraries(quality INTERFACE thread_coordination benchmark_base)
target_include_directories(quality INTERFACE "utils" "wrapper" "${CMAKE_BINARY_DIR}")
target_compile_features(quality INTERFACE cxx_std_17)
target_sources(quality INTERFACE stress_test.cpp)
target_compile_definitions(quality INTERFACE QUALITY)

add_library(throughput INTERFACE)
target_link_libraries(throughput INTERFACE thread_coordination benchmark_base)
target_include_directories(throughput INTERFACE "utils" "wrapper" "${CMAKE_BINARY_DIR}")
target_compile_features(throughput INTERFACE cxx_std_17)
target_sources(throughput INTERFACE stress_test.cpp)
target_compile_options(throughput INTERFACE -fno-omit-frame-pointer)
target_compile_definitions(throughput INTERFACE THROUGHPUT)

add_library(sssp INTERFACE)
target_link_libraries(sssp INTERFACE thread_coordination benchmark_base)
target_include_directories(sssp INTERFACE "utils" "wrapper" "${CMAKE_BINARY_DIR}")
target_compile_features(sssp INTERFACE cxx_std_17)
target_sources(sssp INTERFACE shortest_path.cpp)

add_library(knapsack INTERFACE)
target_link_libraries(knapsack INTERFACE thread_coordination benchmark_base)
target_include_directories(knapsack INTERFACE "utils" "wrapper" "${CMAKE_BINARY_DIR}")
target_compile_features(knapsack INTERFACE cxx_std_17)
target_sources(knapsack INTERFACE knapsack.cpp)

add_library(knapsack_external INTERFACE)
target_link_libraries(knapsack_external INTERFACE thread_coordination benchmark_base)
target_include_directories(knapsack_external INTERFACE "utils" "wrapper" "${CMAKE_BINARY_DIR}")
target_compile_features(knapsack_external INTERFACE cxx_std_17)
target_sources(knapsack_external INTERFACE knapsack_external.cpp)

function(add_benchmarks_for priority_queue_target)
  get_property(name TARGET ${priority_queue_target} PROPERTY pq_name)

  add_executable(quality_${name})
  set_target_properties(quality_${name} PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
  target_link_libraries(quality_${name} PRIVATE ${priority_queue_target} quality)

  add_executable(throughput_${name})
  set_target_properties(throughput_${name} PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
  target_link_libraries(throughput_${name} PRIVATE ${priority_queue_target} throughput)

  add_executable(sssp_${name})
  set_target_properties(sssp_${name} PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
  target_link_libraries(sssp_${name} PRIVATE ${priority_queue_target} sssp)

  add_executable(knapsack_${name})
  set_target_properties(knapsack_${name} PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
  target_link_libraries(knapsack_${name} PRIVATE ${priority_queue_target} knapsack)

  add_executable(knapsack_external_${name})
  set_target_properties(knapsack_external_${name} PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
  target_link_libraries(knapsack_external_${name} PRIVATE ${priority_queue_target} knapsack_external tbbmalloc)

  if(BUILD_TESTING)
    add_test(NAME quality_${name}_small COMMAND /bin/bash -c "$<TARGET_FILE:quality_${name}> -p 10000 -n 10000 -j 1 | $<TARGET_FILE:evaluate_quality> -v")
    add_test(NAME quality_${name}_large COMMAND /bin/bash -c "$<TARGET_FILE:quality_${name}> -j 8 | $<TARGET_FILE:evaluate_quality> -v")
    add_test(NAME sssp_${name}_small COMMAND /bin/bash -c "$<TARGET_FILE:sssp_${name}> -j 1 -f data/NY_graph.gr -c data/NY_solution.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME sssp_${name}_large COMMAND /bin/bash -c "$<TARGET_FILE:sssp_${name}> -j 8 -f data/CAL_graph.gr -c data/CAL_solution.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME knapsack_${name}_small COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_${name}> -j 1 -f  data/knapsack_small.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME knapsack_${name}_large COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_${name}> -j 8 -f data/knapsack_large" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME knapsack_external_${name}_small COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_external_${name}> -j 1 -f  data/knapsack_external_small.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME knapsack_external_${name}_large COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_external_${name}> -j 8 -f data/knapsack_external_large" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
endfunction()

function(add_mq_configuration target bufsize)
  get_property(name TARGET ${target} PROPERTY pq_name)
  set(target_name ${name}_bs_${bufsize})
  add_library(${target_name} INTERFACE)
  target_link_libraries(${target_name} INTERFACE ${target})
  target_compile_definitions(
    ${target_name}
    INTERFACE
    MQ_CONFIG_DELETION_BUFFER_SIZE=${bufsize}
    MQ_CONFIG_INSERTION_BUFFER_SIZE=${bufsize}
  )
  set_property(TARGET ${target_name} PROPERTY pq_name ${target_name})
  add_benchmarks_for(${target_name})
endfunction()

add_mq_configuration(mq_random 1)
add_mq_configuration(mq_random 2)
add_mq_configuration(mq_random 4)
add_mq_configuration(mq_random 8)
add_mq_configuration(mq_random 16)
add_mq_configuration(mq_random 64)

add_mq_configuration(mq_sticky 16)
add_mq_configuration(mq_sticky 64)
add_mq_configuration(mq_sticky 512)

add_benchmarks_for(capq_wrapper)
add_benchmarks_for(klsm256_wrapper)
add_benchmarks_for(klsm1024_wrapper)
add_benchmarks_for(klsm4096_wrapper)
add_benchmarks_for(linden_wrapper)
add_benchmarks_for(spraylist_wrapper)

# add_executable(distribution)
# target_link_libraries(distribution thread_coordination benchmark_base multiqueue::multiqueue)
# target_include_directories(distribution PRIVATE "utils" "${CMAKE_BINARY_DIR}")
# target_include_directories(distribution PRIVATE "${CMAKE_SOURCE_DIR}/external")
# target_compile_features(distribution PRIVATE cxx_std_17)
# target_sources(distribution PRIVATE distribution.cpp)

if(experiments_BUILD_TESTS)
  add_subdirectory(unit_tests)
endif()
