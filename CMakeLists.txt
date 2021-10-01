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

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

add_library("experiment_base" INTERFACE)

string(
  CONCAT WARNING_FLAGS
  "-Wall;"
  "-Wextra;"
  "-Wpedantic;"
)

target_compile_options(
  "experiment_base" INTERFACE ${WARNING_FLAGS}
)

add_subdirectory(multiqueue_simple)

# Find the pthread library
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)

find_package(TBB)

function(generate_mq_target name compiler_define)
  add_library(${name} INTERFACE)
  target_compile_definitions(${name} INTERFACE ${compiler_define})
  target_link_libraries(${name} INTERFACE multiqueue::multiqueue)
  set_property(TARGET ${name} PROPERTY pq_name ${name})
endfunction()

generate_mq_target(mq PQ_MQ_RANDOM)
generate_mq_target(nbmq PQ_MQ_STICKY)

set(EXTERNAL_DIR "${CMAKE_CURRENT_LIST_DIR}/external")

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
target_compile_definitions(throughput INTERFACE THROUGHPUT)

add_library(shortest_path INTERFACE)
target_link_libraries(shortest_path INTERFACE thread_coordination benchmark_base)
target_include_directories(shortest_path INTERFACE "utils" "wrapper" "${CMAKE_BINARY_DIR}")
target_compile_features(shortest_path INTERFACE cxx_std_17)
target_sources(shortest_path INTERFACE shortest_path.cpp)

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
  add_executable(${name}_quality)
  if(BUILD_TESTING)
    add_test(NAME ${name}_quality_small COMMAND /bin/bash -c "$<TARGET_FILE:${name}_quality> -n 0 -j 1 | $<TARGET_FILE:evaluate_quality> -v")
    add_test(NAME ${name}_quality_large COMMAND /bin/bash -c "$<TARGET_FILE:${name}_quality> -n 1000000 -j 8 | $<TARGET_FILE:evaluate_quality> -v")
  endif()
  target_link_libraries(${name}_quality PRIVATE ${priority_queue_target} quality)
  add_executable(${name}_throughput)
  target_link_libraries(${name}_throughput PRIVATE ${priority_queue_target} throughput)
  add_executable(${name}_shortest_path)
  target_link_libraries(${name}_shortest_path PRIVATE ${priority_queue_target} shortest_path)
  if(BUILD_TESTING)
    add_test(NAME ${name}_sssp_small COMMAND /bin/bash -c "$<TARGET_FILE:${name}_shortest_path> -j 1 -f data/NY_graph.gr -c data/NY_solution.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME ${name}_sssp_large COMMAND /bin/bash -c "$<TARGET_FILE:${name}_shortest_path> -j 8 -f data/CAL_graph.gr -c data/CAL_solution.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
  add_executable(${name}_knapsack)
  target_link_libraries(${name}_knapsack PRIVATE ${priority_queue_target} knapsack)
  if(BUILD_TESTING)
    add_test(NAME ${name}_knapsack_small COMMAND /bin/bash -c "$<TARGET_FILE:${name}_knapsack> -j 1 -f  data/knapsack_small.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME ${name}_knapsack_large COMMAND /bin/bash -c "$<TARGET_FILE:${name}_knapsack> -j 8 -f data/knapsack_large" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
  add_executable(${name}_knapsack_external)
  set_target_properties(${name}_knapsack_external PROPERTIES DEBUG_POSTFIX _debug)
  target_link_libraries(${name}_knapsack_external PRIVATE ${priority_queue_target} knapsack_external tbbmalloc)
  target_compile_definitions(${name}_knapsack_external PRIVATE KNAPSACK_DEBUG_PRINT)
endfunction()

function(add_mq_configuration c k ibs dbs)
  get_property(name TARGET mq PROPERTY pq_name)
  set(target_name ${name}_c_${c}_k_${k}_ibs_${ibs}_dbs_${dbs})
  add_library(${target_name} INTERFACE)
  target_link_libraries(${target_name} INTERFACE mq)
  target_compile_definitions(${target_name}
    INTERFACE
      MQ_CONFIG_C=${c}
      MQ_CONFIG_K=${k}
      MQ_CONFIG_DELETION_BUFFER_SIZE=${dbs}
      MQ_CONFIG_INSERTION_BUFFER_SIZE=${ibs}
      )
  set_property(TARGET ${target_name} PROPERTY pq_name ${target_name})
  add_benchmarks_for(${target_name})
endfunction()

function(add_nbmq_configuration c k)
  get_property(name TARGET nbmq PROPERTY pq_name)
  set(target_name ${name}_c_${c}_k_${k})
  add_library(${target_name} INTERFACE)
  target_link_libraries(${target_name} INTERFACE nbmq)
  target_compile_definitions(${target_name}
    INTERFACE
      MQ_CONFIG_C=${c}
      MQ_CONFIG_K=${k}
      )
  set_property(TARGET ${target_name} PROPERTY pq_name ${target_name})
  add_benchmarks_for(${target_name})
endfunction()

function(add_asmq_configuration c k ibs dbs)
  get_property(name TARGET asmq PROPERTY pq_name)
  set(target_name ${name}_c_${c}_k_${k}_ibs_${ibs}_dbs_${dbs})
  add_library(${target_name} INTERFACE)
  target_link_libraries(${target_name} INTERFACE asmq)
  target_compile_definitions(${target_name}
    INTERFACE
      MQ_CONFIG_C=${c}
      MQ_CONFIG_K=${k}
      MQ_CONFIG_DELETION_BUFFER_SIZE=${dbs}
      MQ_CONFIG_INSERTION_BUFFER_SIZE=${ibs}
      )
  set_property(TARGET ${target_name} PROPERTY pq_name ${target_name})
  add_benchmarks_for(${target_name})
endfunction()

add_nbmq_configuration(2 1)
add_nbmq_configuration(4 1)
add_nbmq_configuration(4 4)
add_nbmq_configuration(8 1)
add_nbmq_configuration(8 8)
add_nbmq_configuration(8 16)
add_nbmq_configuration(16 1)
add_nbmq_configuration(16 16)
add_nbmq_configuration(16 64)

add_mq_configuration(2 1 16 16)
add_mq_configuration(4 1 16 16)
add_mq_configuration(4 4 16 16)
add_mq_configuration(8 1 16 16)
add_mq_configuration(8 8 16 16)
add_mq_configuration(8 16 16 16)
add_mq_configuration(16 1 16 16)
add_mq_configuration(16 16 16 16)
add_mq_configuration(16 64 16 16)

add_asmq_configuration(3 4 16 16)
add_asmq_configuration(3 16 16 16)
add_asmq_configuration(3 64 16 16)

add_benchmarks_for(capq_wrapper)
add_benchmarks_for(klsm256_wrapper)
add_benchmarks_for(klsm1024_wrapper)
add_benchmarks_for(dlsm_wrapper)
add_benchmarks_for(linden_wrapper)
add_benchmarks_for(spraylist_wrapper)

add_executable(distribution)
target_link_libraries(distribution thread_coordination benchmark_base multiqueue::multiqueue)
target_include_directories(distribution PRIVATE "utils" "${CMAKE_BINARY_DIR}")
target_include_directories(distribution PRIVATE "${CMAKE_SOURCE_DIR}/external")
target_compile_features(distribution PRIVATE cxx_std_17)
target_sources(distribution PRIVATE distribution.cpp)

if (experiments_BUILD_TESTS)
  add_subdirectory(unit_tests)
endif()
