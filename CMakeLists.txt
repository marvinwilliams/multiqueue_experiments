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

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)
find_package(TBB)
find_package(PAPI)

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
add_subdirectory(common)
add_subdirectory(tools)
add_subdirectory(wrapper)

if(MULTIQUEUE_EXP_BUILD_TESTS AND TARGET Catch2::Catch2WithMain)
  include(CTest)
  add_subdirectory(tests)
endif()

add_library(benchmark_base INTERFACE)
target_include_directories(benchmark_base INTERFACE "${CMAKE_SOURCE_DIR}/third_party" "${CMAKE_SOURCE_DIR}/common")
target_compile_options(benchmark_base INTERFACE $<$<CONFIG:Release>:-march=native>)
if(CMAKE_COMPILER_IS_GNUCXX)
  target_link_options(benchmark_base INTERFACE $<$<CONFIG:Release>:-flto>)
endif()
target_compile_options(benchmark_base INTERFACE ${WARNING_FLAGS})
target_compile_definitions(benchmark_base INTERFACE L1_CACHE_LINESIZE=${L1_CACHE_LINESIZE})
target_compile_definitions(benchmark_base INTERFACE PAGESIZE=${PAGESIZE})
target_compile_features(benchmark_base INTERFACE cxx_std_17)

if(MULTIQUEUE_EXP_PC)
  if(PAPI_FOUND)
    target_link_libraries(benchmark_base INTERFACE PAPI::PAPI)
    target_compile_definitions(benchmark_base INTERFACE USE_PAPI)
  else()
    message(SEND_ERROR "PAPI library not found")
  endif()
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

add_library(using_multiqueue INTERFACE)
target_link_libraries(using_multiqueue INTERFACE multiqueue::multiqueue)
target_compile_definitions(using_multiqueue INTERFACE PQ_MQ)

function(add_mq_variant ibs dbs d)
  set(name mq_ibs_${ibs}_dbs_${dbs}_d_${d})
  add_library(${name} INTERFACE)
  target_compile_definitions(${name} INTERFACE MQ_STICK_POLICY_NONE MQ_INSERTION_BUFFERSIZE=${ibs} MQ_DELETION_BUFFERSIZE=${dbs} MQ_HEAP_ARITY=${d})
  target_link_libraries(${name} INTERFACE using_multiqueue)
  list(APPEND mq_targets ${name})
  list(APPEND mq_explore_targets ${name})
  set(mq_explore_targets ${mq_explore_targets} PARENT_SCOPE)
endfunction()

foreach(ibs 4 8 16 64 128)
  foreach(dbs 4 8 16 64 128)
    add_mq_variant(${ibs} ${dbs} 8)
  endforeach()
endforeach()
foreach(d 2 4 16 64)
  add_mq_variant(16 16 ${d})
endforeach()

add_library(mq_std INTERFACE)
target_compile_definitions(mq_std INTERFACE MQ_STICK_POLICY_NONE MQ_PQ_STD MQ_DISABLE_BUFFERING)
target_link_libraries(mq_std INTERFACE using_multiqueue)
list(APPEND mq_explore_targets mq_std)

add_library(mq_std_buffered INTERFACE)
target_compile_definitions(mq_std_buffered INTERFACE MQ_STICK_POLICY_NONE MQ_PQ_STD)
target_link_libraries(mq_std_buffered INTERFACE using_multiqueue)
list(APPEND mq_explore_targets mq_std_buffered)

add_library(mq_none INTERFACE)
target_compile_definitions(mq_none INTERFACE MQ_STICK_POLICY_NONE)
target_link_libraries(mq_none INTERFACE using_multiqueue)
list(APPEND mq_targets mq_none)

add_library(mq_random INTERFACE)
target_compile_definitions(mq_random INTERFACE MQ_STICK_POLICY_RANDOM)
target_link_libraries(mq_random INTERFACE using_multiqueue)
list(APPEND mq_targets mq_random)

add_library(mq_random_strict INTERFACE)
target_compile_definitions(mq_random_strict INTERFACE MQ_STICK_POLICY_RANDOM_STRICT)
target_link_libraries(mq_random_strict INTERFACE using_multiqueue)
list(APPEND mq_targets mq_random_strict)

add_library(mq_swapping INTERFACE)
target_compile_definitions(mq_swapping INTERFACE MQ_STICK_POLICY_SWAPPING)
target_link_libraries(mq_swapping INTERFACE using_multiqueue)
list(APPEND mq_targets mq_swapping)

add_library(mq_swapping_lazy INTERFACE)
target_compile_definitions(mq_swapping_lazy INTERFACE MQ_STICK_POLICY_SWAPPING_LAZY)
target_link_libraries(mq_swapping_lazy INTERFACE using_multiqueue)
list(APPEND mq_targets mq_swapping_lazy)

add_library(mq_swapping_blocking INTERFACE)
target_compile_definitions(mq_swapping_blocking INTERFACE MQ_STICK_POLICY_SWAPPING_BLOCKING)
target_link_libraries(mq_swapping_blocking INTERFACE using_multiqueue)
list(APPEND mq_targets mq_swapping_blocking)

add_library(mq_perm INTERFACE)
target_compile_definitions(mq_perm INTERFACE MQ_STICK_POLICY_PERMUTATION)
target_link_libraries(mq_perm INTERFACE using_multiqueue)
list(APPEND mq_targets mq_perm)

add_library(mf INTERFACE)
target_compile_definitions(mf INTERFACE PQ_MF)
target_link_libraries(mf INTERFACE multififo::multififo)

set(wrapper_targets capq klsm256 klsm1024 klsm4096 linden spraylist tbb_pq)

add_library(throughput INTERFACE)
target_sources(throughput INTERFACE throughput.cpp)
target_link_libraries(throughput INTERFACE thread_coordination benchmark_base)

foreach(target ${mq_explore_targets} ${mq_targets} ${wrapper_targets})
  add_executable(throughput_${target})
  target_link_libraries(throughput_${target} PRIVATE ${target} throughput)
endforeach()

add_library(quality INTERFACE)
target_sources(quality INTERFACE quality.cpp)
target_link_libraries(quality INTERFACE thread_coordination benchmark_base)

foreach(target ${mq_targets} ${wrapper_targets})
  add_executable(quality_${target})
  target_link_libraries(quality_${target} PRIVATE ${target} quality)
  if(BUILD_TESTING)
    add_test(NAME quality_${target}_small COMMAND /bin/bash -c "$<TARGET_FILE:quality_${target}> -p 10000 -n 10000 -j 1 | $<TARGET_FILE:evaluate_quality> -v")
    add_test(NAME quality_${target}_large COMMAND /bin/bash -c "$<TARGET_FILE:quality_${target}> -j 8 | $<TARGET_FILE:evaluate_quality> -v")
  endif()
endforeach()

add_library(insert_delete INTERFACE)
target_sources(insert_delete INTERFACE insert_delete.cpp)
target_link_libraries(insert_delete INTERFACE thread_coordination benchmark_base)

foreach(target ${mq_explore_targets} ${mq_targets} ${wrapper_targets})
  add_executable(insert_delete_${target})
  target_link_libraries(insert_delete_${target} PRIVATE ${target} throughput)
endforeach()

add_library(sssp INTERFACE)
target_sources(sssp INTERFACE shortest_path.cpp)
target_link_libraries(sssp INTERFACE thread_coordination benchmark_base)

foreach(target ${mq_targets} ${wrapper_targets} tbb_q mf)
  add_executable(sssp_${target})
  target_link_libraries(sssp_${target} PRIVATE ${target} sssp)
  if(BUILD_TESTING)
    add_test(NAME sssp_${target}_small COMMAND /bin/bash -c "$<TARGET_FILE:sssp_${target}> -j 1 -f data/NY_graph.gr -c data/NY_solution.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME sssp_${target}_large COMMAND /bin/bash -c "$<TARGET_FILE:sssp_${target}> -j 8 -f data/CAL_graph.gr -c data/CAL_solution.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
endforeach()

add_library(knapsack_packed INTERFACE)
target_sources(knapsack_packed INTERFACE knapsack.cpp)
target_link_libraries(knapsack_packed INTERFACE thread_coordination benchmark_base)
target_compile_definitions(knapsack_packed INTERFACE HINT_MODE PACKED_VALUE)

add_library(knapsack_heap INTERFACE)
target_sources(knapsack_heap INTERFACE knapsack.cpp)
target_link_libraries(knapsack_heap INTERFACE thread_coordination benchmark_base)
target_compile_definitions(knapsack_heap INTERFACE HINT_MODE HEAP_VALUE)

add_library(knapsack_explicit INTERFACE)
target_sources(knapsack_explicit INTERFACE knapsack.cpp)
target_link_libraries(knapsack_explicit INTERFACE thread_coordination benchmark_base)
target_compile_definitions(knapsack_explicit INTERFACE HINT_MODE EXPLICIT_VALUE)

foreach(target ${targets})
  add_executable(knapsack_packed_${target})
  target_link_libraries(knapsack_packed_${target} PRIVATE ${target} knapsack_packed)
  add_dependencies(knapsack_all knapsack_packed_${target})
  if(BUILD_TESTING)
    add_test(NAME knapsack_packed_${target}_small COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_packed_${target}> -j 1 -f  data/knapsack_small.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME knapsack_packed_${target}_large COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_packed_${target}> -j 8 -f data/knapsack_large.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
endforeach()

foreach(target ${targets})
  add_executable(knapsack_heap_${target})
  target_link_libraries(knapsack_heap_${target} PRIVATE ${target} knapsack_heap tbbmalloc)
  add_dependencies(knapsack_all knapsack_heap_${target})
  if(BUILD_TESTING)
    add_test(NAME knapsack_heap_${target}_small COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_heap_${target}> -j 1 -f  data/knapsack_small.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME knapsack_heap_${target}_large COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_heap_${target}> -j 8 -f data/knapsack_large.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
endforeach()

foreach(target ${targets})
  add_executable(knapsack_explicit_${target})
  target_link_libraries(knapsack_explicit_${target} PRIVATE ${target} knapsack_explicit)
  add_dependencies(knapsack_all knapsack_explicit_${target})
  if(BUILD_TESTING)
    add_test(NAME knapsack_explicit_${target}_small COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_explicit_${target}> -j 1 -f  data/knapsack_small.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME knapsack_explicit_${target}_large COMMAND /bin/bash -c "$<TARGET_FILE:knapsack_explicit_${target}> -j 8 -f data/knapsack_large.txt" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
endforeach()
