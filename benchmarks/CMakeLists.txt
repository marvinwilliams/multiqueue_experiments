add_library(benchmark_base INTERFACE)
target_include_directories(
  benchmark_base INTERFACE "${CMAKE_SOURCE_DIR}/third_party"
                           "${CMAKE_CURRENT_SOURCE_DIR}/util")
target_compile_options(
  benchmark_base
  INTERFACE
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:${WARNING_FLAGS}>$<$<CXX_COMPILER_ID:MSVC>:
    /W4>)
target_compile_definitions(benchmark_base
                           INTERFACE L1_CACHE_LINESIZE=${L1_CACHE_LINESIZE})
target_compile_definitions(benchmark_base INTERFACE PAGESIZE=${PAGESIZE})
target_compile_features(benchmark_base INTERFACE cxx_std_17)

if(MULTIQUEUE_EXP_PERFORMANCE_COUNTERS)
  target_link_libraries(benchmark_base INTERFACE PAPI::PAPI)
  target_compile_definitions(benchmark_base INTERFACE WITH_PAPI)
endif()
if(MULTIQUEUE_EXP_PGO)
  if(MULTIQUEUE_EXP_PGO_USE)
    target_compile_options(benchmark_base INTERFACE "-fprofile-use"
                                                    "-fprofile-correction")
  else()
    target_compile_options(benchmark_base INTERFACE "-fprofile-generate")
    target_link_options(benchmark_base INTERFACE "-fprofile-generate")
  endif()
endif()
if(MULTIQUEUE_EXP_COVERAGE)
  target_compile_options(benchmark_base INTERFACE "--coverage")
  target_link_libraries(benchmark_base INTERFACE "gcov")
endif()
if(MULTIQUEUE_EXP_ASAN)
  target_compile_options(benchmark_base INTERFACE -fsanitize=address)
  target_link_options(benchmark_base INTERFACE -fsanitize=address)
endif()
if(MULTIQUEUE_EXP_UBSAN)
  target_compile_options(benchmark_base INTERFACE -fsanitize=undefined)
  target_link_options(benchmark_base INTERFACE -fsanitize=undefined)
endif()
if(MULTIQUEUE_EXP_TSAN)
  target_compile_options(benchmark_base INTERFACE -fsanitize=thread)
  target_link_options(benchmark_base INTERFACE -fsanitize=thread)
endif()

add_library(mq_interface INTERFACE)
target_link_libraries(mq_interface INTERFACE multiqueue::multiqueue)
target_compile_definitions(mq_interface INTERFACE PQ_MQ)

function(add_tuning_variant buffer_size degree)
  set(name mq_bs_${buffer_size}_d_${degree})
  add_library(${name} INTERFACE)
  target_compile_definitions(
    ${name} INTERFACE STICK_POLICY_NONE INSERTION_BUFFERSIZE=${buffer_size}
                      DELETION_BUFFERSIZE=${buffer_size} HEAP_ARITY=${degree})
  target_link_libraries(${name} INTERFACE mq_interface)
  list(APPEND MQ_TUNING_TARGETS ${name})
  set(MQ_TUNING_TARGETS
      ${MQ_TUNING_TARGETS}
      PARENT_SCOPE)
endfunction()

foreach(
  buffer_size
  0
  4
  8
  16
  64
  128
  512)
  foreach(degree 2 4 8 16 64)
    add_tuning_variant(${buffer_size} ${degree})
  endforeach()
endforeach()

add_library(mq_bs_0_std INTERFACE)
target_compile_definitions(
  mq_bs_0_std INTERFACE STICK_POLICY_NONE USE_STD_PQ INSERTION_BUFFERSIZE=0
                        DELETION_BUFFERSIZE=0)
target_link_libraries(mq_bs_0_std INTERFACE mq_interface)
list(APPEND MQ_TUNING_TARGETS mq_bs_0_std)

add_library(mq_none INTERFACE)
target_compile_definitions(mq_none INTERFACE STICK_POLICY_NONE)
target_link_libraries(mq_none INTERFACE mq_interface)
list(APPEND MQ_TARGETS mq_none)

add_library(mq_none1 INTERFACE)
target_compile_definitions(mq_none1 INTERFACE STICK_POLICY_NONE
                                              MQ_NUM_COMPARE=1)
target_link_libraries(mq_none1 INTERFACE mq_interface)
list(APPEND MQ_TUNING_TARGETS mq_none1)

add_library(mq_none2 INTERFACE)
target_compile_definitions(mq_none2 INTERFACE STICK_POLICY_NONE
                                              MQ_NUM_COMPARE=2)
target_link_libraries(mq_none2 INTERFACE mq_interface)
list(APPEND MQ_TUNING_TARGETS mq_none2)

add_library(mq_none4 INTERFACE)
target_compile_definitions(mq_none4 INTERFACE STICK_POLICY_NONE
                                              MQ_NUM_COMPARE=4)
target_link_libraries(mq_none4 INTERFACE mq_interface)
list(APPEND MQ_TUNING_TARGETS mq_none4)

add_library(mq_none8 INTERFACE)
target_compile_definitions(mq_none8 INTERFACE STICK_POLICY_NONE
                                              MQ_NUM_COMPARE=8)
target_link_libraries(mq_none8 INTERFACE mq_interface)
list(APPEND MQ_TUNING_TARGETS mq_none8)

add_library(mq_random INTERFACE)
target_compile_definitions(mq_random INTERFACE STICK_POLICY_RANDOM)
target_link_libraries(mq_random INTERFACE mq_interface)
list(APPEND MQ_TARGETS mq_random)

add_library(mq_swapping INTERFACE)
target_compile_definitions(mq_swapping INTERFACE STICK_POLICY_SWAPPING)
target_link_libraries(mq_swapping INTERFACE mq_interface)
list(APPEND MQ_TARGETS mq_swapping)

add_library(mq_perm INTERFACE)
target_compile_definitions(mq_perm INTERFACE STICK_POLICY_PERMUTATION)
target_link_libraries(mq_perm INTERFACE mq_interface)
list(APPEND MQ_TARGETS mq_perm)

set(WRAPPER_TARGETS
    klsm4
    klsm256
    klsm1024
    klsm4096
    linden
    capq
    spraylist
    tbb_pq)

add_library(throughput INTERFACE)
target_sources(throughput INTERFACE throughput.cpp util/threading.cpp)
target_link_libraries(throughput INTERFACE benchmark_base)

foreach(target ${MQ_TUNING_TARGETS} ${MQ_TARGETS} ${WRAPPER_TARGETS})
  add_executable(throughput_${target})
  target_link_libraries(throughput_${target} PRIVATE ${target} throughput)
endforeach()

foreach(target ${MQ_TUNING_TARGETS} ${MQ_TARGETS} klsm4 klsm256 klsm1024 klsm4096 tbb_pq)
  add_executable(throughput_${target}_uint8)
  target_link_libraries(throughput_${target}_uint8 PRIVATE ${target} throughput)
  target_compile_definitions(throughput_${target}_uint8 PRIVATE -DUSE_UINT8)
endforeach()

add_custom_target(throughput_mq_tuning)
foreach(target ${MQ_TUNING_TARGETS})
  add_dependencies(throughput_mq_tuning throughput_${target})
endforeach()

add_custom_target(throughput_mq)
foreach(target ${MQ_TARGETS})
  add_dependencies(throughput_mq throughput_${target})
endforeach()

add_custom_target(throughput_wrapper)
foreach(target ${WRAPPER_TARGETS})
  add_dependencies(throughput_wrapper throughput_${target})
endforeach()

add_custom_target(throughput_all)
foreach(target ${MQ_TARGETS} ${WRAPPER_TARGETS})
  add_dependencies(throughput_all throughput_${target})
endforeach()

add_library(quality INTERFACE)
target_sources(quality INTERFACE throughput.cpp util/threading.cpp
                                 util/evaluate_quality.cpp)
target_compile_definitions(quality INTERFACE -DQUALITY)
target_link_libraries(quality INTERFACE benchmark_base Threads::Threads)

foreach(target ${MQ_TUNING_TARGETS} ${MQ_TARGETS} ${WRAPPER_TARGETS})
  add_executable(quality_${target})
  target_link_libraries(quality_${target} PRIVATE ${target} quality)
  if(BUILD_TESTING)
    add_test(
      NAME quality_${target}_small
      COMMAND
        /bin/bash -c
        "$<TARGET_FILE:quality_${target}> -p 10000 -n 10000 -j 1 | $<TARGET_FILE:evaluate_quality> -v"
    )
    add_test(
      NAME quality_${target}_large
      COMMAND
        /bin/bash -c
        "$<TARGET_FILE:quality_${target}> -j 8 | $<TARGET_FILE:evaluate_quality> -v"
    )
  endif()
endforeach()

add_custom_target(quality_mq)
foreach(target ${MQ_TARGETS})
  add_dependencies(quality_mq quality_${target})
endforeach()

add_custom_target(quality_all)
foreach(target ${MQ_TARGETS} ${WRAPPER_TARGETS})
  add_dependencies(quality_all quality_${target})
endforeach()

add_library(shortest_path INTERFACE)
target_sources(shortest_path INTERFACE shortest_path.cpp util/threading.cpp)
target_link_libraries(shortest_path INTERFACE benchmark_base Threads::Threads)

foreach(target ${MQ_TARGETS} ${WRAPPER_TARGETS})
  add_executable(shortest_path_${target})
  target_link_libraries(shortest_path_${target} PRIVATE ${target} shortest_path)
  if(BUILD_TESTING)
    add_test(
      NAME shortest_path_${target}_small
      COMMAND
        /bin/bash -c
        "$<TARGET_FILE:shortest_path_${target}> data/NY_graph.gr data/NY_solution.txt -j 1"
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(
      NAME shortest_path_${target}_large
      COMMAND
        /bin/bash -c
        "$<TARGET_FILE:shortest_path_${target}> data/CAL_graph.gr data/CAL_solution.txt -j 8"
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
endforeach()

add_custom_target(shortest_path_all)
foreach(target ${MQ_TARGETS} ${WRAPPER_TARGETS})
  add_dependencies(shortest_path_all shortest_path_${target})
endforeach()
add_dependencies(shortest_path_all shortest_path_seq)

add_executable(shortest_path_seq shortest_path_seq.cpp)
target_link_libraries(shortest_path_seq PRIVATE benchmark_base)

add_library(knapsack_packed INTERFACE)
target_sources(knapsack_packed INTERFACE knapsack.cpp util/threading.cpp)
target_link_libraries(knapsack_packed INTERFACE benchmark_base Threads::Threads)
target_compile_definitions(knapsack_packed INTERFACE HINT_MODE PACKED_VALUE)

add_library(knapsack_heap INTERFACE)
target_sources(knapsack_heap INTERFACE knapsack.cpp utils/threading.cpp)
target_link_libraries(knapsack_heap INTERFACE benchmark_base Threads::Threads)
target_compile_definitions(knapsack_heap INTERFACE HINT_MODE HEAP_VALUE)

add_library(knapsack_explicit INTERFACE)
target_sources(knapsack_explicit INTERFACE knapsack.cpp utils/threading.cpp)
target_link_libraries(knapsack_explicit INTERFACE benchmark_base
                                                  Threads::Threads)
target_compile_definitions(knapsack_explicit INTERFACE HINT_MODE EXPLICIT_VALUE)

foreach(target ${targets})
  add_executable(knapsack_packed_${target})
  target_link_libraries(knapsack_packed_${target} PRIVATE ${target}
                                                          knapsack_packed)
  add_dependencies(knapsack_all knapsack_packed_${target})
  if(BUILD_TESTING)
    add_test(
      NAME knapsack_packed_${target}_small
      COMMAND
        /bin/bash -c
        "$<TARGET_FILE:knapsack_packed_${target}> -j 1 -f  data/knapsack_small.txt"
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(
      NAME knapsack_packed_${target}_large
      COMMAND
        /bin/bash -c
        "$<TARGET_FILE:knapsack_packed_${target}> -j 8 -f data/knapsack_large.txt"
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
endforeach()

foreach(target ${targets})
  add_executable(knapsack_heap_${target})
  target_link_libraries(knapsack_heap_${target} PRIVATE ${target} knapsack_heap
                                                        tbbmalloc)
  add_dependencies(knapsack_all knapsack_heap_${target})
  if(BUILD_TESTING)
    add_test(
      NAME knapsack_heap_${target}_small
      COMMAND
        /bin/bash -c
        "$<TARGET_FILE:knapsack_heap_${target}> -j 1 -f  data/knapsack_small.txt"
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(
      NAME knapsack_heap_${target}_large
      COMMAND
        /bin/bash -c
        "$<TARGET_FILE:knapsack_heap_${target}> -j 8 -f data/knapsack_large.txt"
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
endforeach()

foreach(target ${targets})
  add_executable(knapsack_explicit_${target})
  target_link_libraries(knapsack_explicit_${target} PRIVATE ${target}
                                                            knapsack_explicit)
  add_dependencies(knapsack_all knapsack_explicit_${target})
  if(BUILD_TESTING)
    add_test(
      NAME knapsack_explicit_${target}_small
      COMMAND
        /bin/bash -c
        "$<TARGET_FILE:knapsack_explicit_${target}> -j 1 -f  data/knapsack_small.txt"
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(
      NAME knapsack_explicit_${target}_large
      COMMAND
        /bin/bash -c
        "$<TARGET_FILE:knapsack_explicit_${target}> -j 8 -f data/knapsack_large.txt"
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
endforeach()

add_executable(knapsack_seq knapsack_seq.cpp)
target_link_libraries(knapsack_seq PRIVATE benchmark_base)
