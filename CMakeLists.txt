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
option("MQ_NO_BUFFERING" "Disable buffering in the mq" OFF)
option("MQ_IMPLICIT_LOCK" "Use implicit locking" OFF)
set("MQ_DBUF_SIZE" "8" CACHE STRING "MQ deletion buffer size")
set("MQ_IBUF_SIZE" "8" CACHE STRING "MQ insertion buffer size")
set("MQ_HEAP_DEGREE" "8" CACHE STRING "MQ heap degree")
set(L1_CACHE_LINESIZE ${DEFAULT_L1_CACHE_LINESIZE} CACHE STRING "Specify the assumed L1 cache linesize (bytes)")
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

project(multiqueue_experiments VERSION 2.0.0 LANGUAGES CXX C)

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
  target_compile_options(benchmark_base INTERFACE -fsanitize=address)
  target_link_options(benchmark_base INTERFACE -fsanitize=address)
endif()
if(experiments_UBSAN)
  target_compile_options(benchmark_base INTERFACE -fsanitize=undefined)
  target_link_options(benchmark_base INTERFACE -fsanitize=undefined)
endif()
if(experiments_TSAN)
  target_compile_options(benchmark_base INTERFACE -fsanitize=thread)
  target_link_options(benchmark_base INTERFACE -fsanitize=thread)
endif()

add_library(quality INTERFACE)
target_link_libraries(quality INTERFACE thread_coordination benchmark_base)
target_compile_features(quality INTERFACE cxx_std_17)
target_sources(quality INTERFACE stress_test.cpp)
target_compile_definitions(quality INTERFACE QUALITY)

add_library(throughput INTERFACE)
target_link_libraries(throughput INTERFACE thread_coordination benchmark_base)
target_compile_features(throughput INTERFACE cxx_std_17)
target_sources(throughput INTERFACE stress_test.cpp)
target_compile_definitions(throughput INTERFACE THROUGHPUT)

add_library(sssp INTERFACE)
target_link_libraries(sssp INTERFACE thread_coordination benchmark_base)
target_compile_features(sssp INTERFACE cxx_std_17)
target_sources(sssp INTERFACE shortest_path.cpp)

add_library(knapsack INTERFACE)
target_link_libraries(knapsack INTERFACE thread_coordination benchmark_base)
target_compile_features(knapsack INTERFACE cxx_std_17)
target_sources(knapsack INTERFACE knapsack.cpp)

add_library(knapsack_external INTERFACE)
target_link_libraries(knapsack_external INTERFACE thread_coordination benchmark_base)
target_compile_features(knapsack_external INTERFACE cxx_std_17)
target_sources(knapsack_external INTERFACE knapsack_external.cpp)

function(add_benchmarks_for priority_queue_target)
  get_property(name TARGET ${priority_queue_target} PROPERTY pq_name)

  add_executable(quality_${name})
  target_link_libraries(quality_${name} PRIVATE ${priority_queue_target} quality)

  add_executable(throughput_${name})
  target_link_libraries(throughput_${name} PRIVATE ${priority_queue_target} throughput)

  add_executable(sssp_${name})
  target_link_libraries(sssp_${name} PRIVATE ${priority_queue_target} sssp)

  add_executable(knapsack_${name})
  target_link_libraries(knapsack_${name} PRIVATE ${priority_queue_target} knapsack)

  add_executable(knapsack_external_${name})
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

add_subdirectory(multiqueue_simple)

function(generate_mq_target name compiler_define)
  add_library(${name} INTERFACE)
  target_compile_definitions(${name} INTERFACE ${compiler_define})
  if(MQ_NO_BUFFERING)
    target_compile_definitions(${name} INTERFACE MQ_NO_BUFFERING)
  endif()
  if(MQ_DBUF_SIZE)
    target_compile_definitions(${name} INTERFACE MQ_DBUF_SIZE=${MQ_DBUF_SIZE})
  endif()
  if(MQ_IBUF_SIZE)
    target_compile_definitions(${name} INTERFACE MQ_IBUF_SIZE=${MQ_IBUF_SIZE})
  endif()
  if(MQ_HEAP_DEGREE)
    target_compile_definitions(${name} INTERFACE MQ_HEAP_DEGREE=${MQ_HEAP_DEGREE})
  endif()
  if(MQ_IMPLICIT_LOCK)
    target_compile_definitions(${name} INTERFACE MQ_IMPLICIT_LOCK)
  endif()
  target_link_libraries(${name} INTERFACE multiqueue::multiqueue)
  set_property(TARGET ${name} PROPERTY pq_name ${name})
endfunction()

generate_mq_target(mq_random PQ_MQ_RANDOM)
generate_mq_target(mq_sticky PQ_MQ_STICKY)
generate_mq_target(mq_swapping PQ_MQ_SWAPPING)
generate_mq_target(mq_perm PQ_MQ_PERM)

add_benchmarks_for(mq_random)
add_benchmarks_for(mq_sticky)
add_benchmarks_for(mq_swapping)
add_benchmarks_for(mq_perm)
add_benchmarks_for(capq_wrapper)
add_benchmarks_for(klsm256_wrapper)
add_benchmarks_for(klsm1024_wrapper)
add_benchmarks_for(klsm4096_wrapper)
add_benchmarks_for(linden_wrapper)
add_benchmarks_for(spraylist_wrapper)

# add_executable(distribution)
# target_link_libraries(distribution thread_coordination benchmark_base multiqueue::multiqueue)
# target_include_directories(distribution PRIVATE "${CMAKE_SOURCE_DIR}/external")
# target_compile_features(distribution PRIVATE cxx_std_17)
# target_sources(distribution PRIVATE distribution.cpp)

if(experiments_BUILD_TESTS)
  if(NOT multiqueue_BUILD_TESTS)
    include(CTest)
    add_subdirectory("multiqueue/external/Catch2" EXCLUDE_FROM_ALL)
  endif()
  add_subdirectory(tests)
endif()
