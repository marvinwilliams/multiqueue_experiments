add_subdirectory(wrapper)

add_library(benchmark_base INTERFACE)
target_include_directories(benchmark_base INTERFACE ".")
target_include_directories(benchmark_base SYSTEM INTERFACE "${CMAKE_SOURCE_DIR}/external")
target_compile_definitions(benchmark_base
                           INTERFACE L1_CACHE_LINESIZE=${L1_CACHE_LINESIZE})
target_compile_definitions(benchmark_base INTERFACE PAGESIZE=${PAGESIZE})
target_compile_definitions(benchmark_base
                           INTERFACE CORES_PER_NUMA_NODE=${CORES_PER_NUMA_NODE})
target_compile_definitions(benchmark_base
                           INTERFACE NUM_NUMA_NODES=${NUM_NUMA_NODES})
target_compile_features(benchmark_base INTERFACE cxx_std_17)
if(PAPI_FOUND)
  target_link_libraries(benchmark_base INTERFACE PAPI::PAPI)
  target_compile_definitions(benchmark_base INTERFACE WITH_PAPI)
endif()

add_library(stress_test_base INTERFACE)
target_sources(stress_test_base INTERFACE stress_test.cpp "util/threading.cpp")
target_link_libraries(stress_test_base INTERFACE benchmark_base
                                                 Threads::Threads)

add_library(dijkstra_base INTERFACE)
target_sources(dijkstra_base INTERFACE dijkstra.cpp "util/threading.cpp")
target_link_libraries(dijkstra_base INTERFACE benchmark_base Threads::Threads)

add_library(knapsack_base INTERFACE)
target_sources(knapsack_base INTERFACE knapsack.cpp "util/threading.cpp")
target_link_libraries(knapsack_base INTERFACE benchmark_base Threads::Threads)

add_library(knapsack_node_base INTERFACE)
target_sources(knapsack_node_base INTERFACE knapsack_node.cpp
                                            "util/threading.cpp")
target_link_libraries(knapsack_node_base INTERFACE benchmark_base
                                                   Threads::Threads)
add_library(knapsack_node_f_base INTERFACE)
target_link_libraries(knapsack_node_f_base INTERFACE knapsack_node_base)
target_compile_definitions(knapsack_node_f_base INTERFACE FLOAT_INSTANCE)

if(EXP_TUNING_TARGETS)
  add_custom_target(tuning_all)
  add_custom_target(tuning_log_all)
endif()

foreach(benchmark stress_test dijkstra knapsack knapsack_node knapsack_node_f)
  add_custom_target(${benchmark}_all)
  add_custom_target(${benchmark}_log_all)
  if(EXP_TUNING_TARGETS)
    add_custom_target(${benchmark}_tuning_all)
    add_custom_target(${benchmark}_tuning_log_all)
  endif()
endforeach()

function(add_benchmark benchmark target)
  add_executable(${benchmark}_${target})
  target_link_libraries(${benchmark}_${target} PRIVATE ${target}
                                                       ${benchmark}_base)

  add_executable(${benchmark}_${target}_log)
  target_link_libraries(${benchmark}_${target}_log PRIVATE ${target}
                                                           ${benchmark}_base)
  target_compile_definitions(${benchmark}_${target}_log PRIVATE LOG_OPERATIONS)
endfunction()

foreach(target ${COMPETITORS})
  foreach(benchmark stress_test dijkstra knapsack)
    add_benchmark(${benchmark} ${target})
    add_dependencies(${benchmark}_all ${benchmark}_${target})
    add_dependencies(${benchmark}_log_all ${benchmark}_${target}_log)
  endforeach()
endforeach()

foreach(target ${GENERIC_COMPETITORS})
  foreach(benchmark knapsack_node knapsack_node_f)
    add_benchmark(${benchmark} ${target})
    add_dependencies(${benchmark}_all ${benchmark}_${target})
    add_dependencies(${benchmark}_log_all ${benchmark}_${target}_log)
  endforeach()
endforeach()

if(EXP_TUNING_TARGETS)
  foreach(target ${TUNING_TARGETS})
    foreach(benchmark stress_test dijkstra knapsack knapsack_node
                      knapsack_node_f)
      add_benchmark(${benchmark} ${target})
      add_dependencies(${benchmark}_tuning_all ${benchmark}_${target})
      add_dependencies(${benchmark}_tuning_log_all ${benchmark}_${target}_log)
    endforeach()
  endforeach()
endif()

add_executable(dijkstra_seq dijkstra_seq.cpp)
target_link_libraries(dijkstra_seq PRIVATE benchmark_base)
add_executable(knapsack_seq knapsack_seq.cpp)
target_link_libraries(knapsack_seq PRIVATE benchmark_base)
add_executable(knapsack_seq_f knapsack_seq.cpp)
target_link_libraries(knapsack_seq_f PRIVATE benchmark_base)
target_compile_definitions(knapsack_seq_f PRIVATE FLOAT_INSTANCE)
