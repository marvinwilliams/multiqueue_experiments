add_library(benchmark_base INTERFACE)
target_include_directories(benchmark_base INTERFACE "." "..")
target_include_directories(benchmark_base SYSTEM
                           INTERFACE "${CMAKE_SOURCE_DIR}/third_party")
target_compile_definitions(benchmark_base
                           INTERFACE L1_CACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE})
target_compile_definitions(benchmark_base INTERFACE PAGE_SIZE=${PAGE_SIZE})
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

add_library(knapsack_float_base INTERFACE)
target_sources(knapsack_float_base INTERFACE knapsack.cpp "util/threading.cpp")
target_link_libraries(knapsack_float_base INTERFACE benchmark_base
                                                    Threads::Threads)
target_compile_definitions(knapsack_float_base INTERFACE FLOAT_INSTANCE)

if(EXP_TUNING_TARGETS)
  add_custom_target(tuning_all)
  add_custom_target(tuning_log_all)
endif()

foreach(benchmark stress_test dijkstra knapsack knapsack_float)
  add_custom_target(${benchmark}_all)
  add_custom_target(${benchmark}_log_all)
  if(EXP_TUNING_TARGETS)
    add_custom_target(${benchmark}_tuning_all)
    add_custom_target(${benchmark}_log_tuning_all)
  endif()
endforeach()

function(add_benchmark benchmark target)
  add_executable(${benchmark}_${target})
  target_link_libraries(${benchmark}_${target} PRIVATE ${target}
                                                       ${benchmark}_base)

  add_executable(${benchmark}_log_${target})
  target_link_libraries(${benchmark}_log_${target} PRIVATE ${target}
                                                           ${benchmark}_base)
  target_compile_definitions(${benchmark}_log_${target} PRIVATE LOG_OPERATIONS)
endfunction()

foreach(target ${COMPETITORS})
  foreach(benchmark stress_test dijkstra knapsack)
    add_benchmark(${benchmark} ${target})
    add_dependencies(${benchmark}_all ${benchmark}_${target})
    add_dependencies(${benchmark}_log_all ${benchmark}_log_${target})
  endforeach()
endforeach()

foreach(target ${GENERIC_COMPETITORS})
  add_benchmark(knapsack_float ${target})
  add_dependencies(knapsack_float_all knapsack_float_${target})
  add_dependencies(knapsack_float_log_all knapsack_float_log_${target})
endforeach()

if(EXP_TUNING_TARGETS)
  foreach(target ${TUNING_TARGETS})
    foreach(benchmark stress_test dijkstra knapsack knapsack_float)
      add_benchmark(${benchmark} ${target})
      add_dependencies(${benchmark}_tuning_all ${benchmark}_${target})
      add_dependencies(${benchmark}_log_tuning_all ${benchmark}_log_${target})
    endforeach()
  endforeach()
endif()

add_executable(dijkstra_seq dijkstra_seq.cpp)
target_link_libraries(dijkstra_seq PRIVATE benchmark_base)
add_executable(knapsack_seq knapsack_seq.cpp)
target_link_libraries(knapsack_seq PRIVATE benchmark_base)
add_executable(knapsack_float_seq knapsack_seq.cpp)
target_link_libraries(knapsack_float_seq PRIVATE benchmark_base)
target_compile_definitions(knapsack_float_seq PRIVATE FLOAT_INSTANCE)
