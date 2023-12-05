add_library(benchmark_base INTERFACE)
target_include_directories(
  benchmark_base INTERFACE "${CMAKE_SOURCE_DIR}/third_party"
                           "${CMAKE_SOURCE_DIR}/util")
target_compile_options(
  benchmark_base
  INTERFACE
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:${WARNING_FLAGS}>$<$<CXX_COMPILER_ID:MSVC>:
    /W4>)
target_compile_definitions(benchmark_base
                           INTERFACE L1_CACHE_LINESIZE=${L1_CACHE_LINESIZE})
target_compile_definitions(benchmark_base INTERFACE PAGESIZE=${PAGESIZE})
target_compile_definitions(benchmark_base
                           INTERFACE CORES_PER_NUMA_NODE=${CORES_PER_NUMA_NODE})
target_compile_definitions(benchmark_base
                           INTERFACE NUM_NUMA_NODES=${NUM_NUMA_NODES})
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

add_library(throughput INTERFACE)
target_sources(throughput INTERFACE synthetic.cpp
                                    "${CMAKE_SOURCE_DIR}/util/threading.cpp")
target_link_libraries(throughput INTERFACE benchmark_base Threads::Threads)

foreach(target ${MQ_TUNING_VARIANTS} ${MQ_VARIANTS} ${COMPETITORS})
  add_executable(throughput_${target})
  target_link_libraries(throughput_${target} PRIVATE ${target} throughput)
endforeach()

foreach(
  target
  ${MQ_VARIANTS}
  klsm4
  klsm256
  klsm1024
  klsm4096
  smq
  tbb_pq)
  add_executable(throughput_${target}_uint8)
  target_link_libraries(throughput_${target}_uint8 PRIVATE ${target} throughput)
  target_compile_definitions(throughput_${target}_uint8 PRIVATE -DUSE_UINT8)
endforeach()

add_custom_target(throughput_mq_tuning)
foreach(target ${MQ_TUNING_VARIANTS})
  add_dependencies(throughput_mq_tuning throughput_${target})
endforeach()

add_custom_target(throughput_mq)
foreach(target ${MQ_VARIANTS})
  add_dependencies(throughput_mq throughput_${target})
endforeach()

add_custom_target(throughput_competitors)
foreach(target ${COMPETITORS})
  add_dependencies(throughput_competitors throughput_${target})
endforeach()

add_custom_target(throughput_all)
foreach(target ${MQ_VARIANTS} ${COMPETITORS})
  add_dependencies(throughput_all throughput_${target})
endforeach()

add_library(quality INTERFACE)
target_sources(
  quality INTERFACE synthetic.cpp "${CMAKE_SOURCE_DIR}/util/threading.cpp"
                    "${CMAKE_SOURCE_DIR}/util/operation_log.cpp")
target_compile_definitions(
  quality INTERFACE -DQUALITY $<$<CONFIG:Debug>:-DREPLAY_TREE_DEBUG>)
target_link_libraries(quality INTERFACE benchmark_base Threads::Threads)

foreach(target ${MQ_TUNING_VARIANTS} ${MQ_VARIANTS} ${COMPETITORS})
  add_executable(quality_${target})
  target_link_libraries(quality_${target} PRIVATE ${target} quality)
  if(BUILD_TESTING)
    add_test(
      NAME quality_${target}_seq
      COMMAND
        /bin/bash -c
        "$<TARGET_FILE:quality_${target}> -p 10000 -n 10000 -j 1 > /dev/null")
    add_test(
      NAME quality_${target}
      COMMAND
        /bin/bash -c
        "$<TARGET_FILE:quality_${target}> -p 10000 -n 10000 -j 8 > /dev/null")
  endif()
endforeach()

add_custom_target(quality_mq_tuning)
foreach(target ${MQ_TUNING_VARIANTS})
  add_dependencies(quality_mq_tuning quality_${target})
endforeach()

add_custom_target(quality_mq)
foreach(target ${MQ_VARIANTS})
  add_dependencies(quality_mq quality_${target})
endforeach()

add_custom_target(quality_all)
foreach(target ${MQ_VARIANTS} ${COMPETITORS})
  add_dependencies(quality_all quality_${target})
endforeach()

add_library(sssp_dijkstra INTERFACE)
target_sources(sssp_dijkstra INTERFACE sssp_dijkstra.cpp
                                       "${CMAKE_SOURCE_DIR}/util/threading.cpp")
target_link_libraries(sssp_dijkstra INTERFACE benchmark_base Threads::Threads)

foreach(target ${MQ_VARIANTS} ${COMPETITORS})
  add_executable(sssp_dijkstra_${target})
  target_link_libraries(sssp_dijkstra_${target} PRIVATE ${target} sssp_dijkstra)
  if(BUILD_TESTING)
    add_test(
      NAME sssp_dijkstra_${target}_seq
      COMMAND /bin/bash -c
              "$<TARGET_FILE:sssp_dijkstra_${target}> data/NY_graph.gr
         -j 1"
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(
      NAME sssp_dijkstra_${target}
      COMMAND /bin/bash -c
              "$<TARGET_FILE:sssp_dijkstra_${target}> data/NY_graph.gr
         -j 8"
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
endforeach()

add_executable(sssp_dijkstra_seq sssp_dijkstra_seq.cpp)
target_link_libraries(sssp_dijkstra_seq PRIVATE benchmark_base)

add_executable(sssp_dijkstra_seq_fifo sssp_dijkstra_seq.cpp)
target_link_libraries(sssp_dijkstra_seq_fifo PRIVATE benchmark_base)
target_compile_definitions(sssp_dijkstra_seq_fifo PRIVATE -DUSE_FIFO)

add_executable(sssp_dijkstra_seq_rev sssp_dijkstra_seq.cpp)
target_link_libraries(sssp_dijkstra_seq_rev PRIVATE benchmark_base)
target_compile_definitions(sssp_dijkstra_seq_rev PRIVATE -DREVERSE_PRIORITY)

add_custom_target(sssp_dijkstra_all)
foreach(target ${MQ_VARIANTS} ${COMPETITORS})
  add_dependencies(sssp_dijkstra_all sssp_dijkstra_${target})
endforeach()
add_dependencies(sssp_dijkstra_all sssp_dijkstra_seq sssp_dijkstra_seq_fifo)

add_library(knapsack INTERFACE)
target_sources(knapsack INTERFACE knapsack.cpp
                                  "${CMAKE_SOURCE_DIR}/util/threading.cpp")
target_link_libraries(knapsack INTERFACE benchmark_base Threads::Threads)

add_custom_target(knapsack_all)
foreach(target ${MQ_VARIANTS} ${COMPETITORS})
  add_executable(knapsack_${target})
  target_link_libraries(knapsack_${target} PRIVATE ${target} knapsack)
  add_dependencies(knapsack_all knapsack_${target})
endforeach()

add_library(knapsack_node INTERFACE)
target_sources(knapsack_node INTERFACE knapsack_node.cpp
                                  "${CMAKE_SOURCE_DIR}/util/threading.cpp")
target_link_libraries(knapsack_node INTERFACE benchmark_base Threads::Threads)

foreach(target ${MQ_VARIANTS} ${GENERIC_COMPETITORS})
  add_executable(knapsack_node_f_${target})
  target_link_libraries(knapsack_node_f_${target} PRIVATE ${target} knapsack_node)
  add_dependencies(knapsack_all knapsack_node_f_${target})
  add_executable(knapsack_node_i_${target})
  target_link_libraries(knapsack_node_i_${target} PRIVATE ${target} knapsack_node)
  target_compile_definitions(knapsack_node_i_${target} PRIVATE INTEGER_INSTANCE)
  add_dependencies(knapsack_all knapsack_node_i_${target})
endforeach()

add_executable(knapsack_seq_f knapsack_seq.cpp)
target_link_libraries(knapsack_seq_f PRIVATE benchmark_base)
add_executable(knapsack_seq_i knapsack_seq.cpp)
target_link_libraries(knapsack_seq_i PRIVATE benchmark_base)
target_compile_definitions(knapsack_seq_i PRIVATE INTEGER_INSTANCE)
add_dependencies(knapsack_all knapsack_seq_f knapsack_seq_i)
