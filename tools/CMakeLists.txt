add_executable(evaluate_quality evaluate_quality.cpp)
target_link_libraries(evaluate_quality PRIVATE multiqueue_internal)
target_include_directories(evaluate_quality PRIVATE "${CMAKE_SOURCE_DIR}/external" "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_features(
  evaluate_quality PRIVATE
  cxx_std_17
)

add_executable(sequential_shortest_path sequential_shortest_path.cpp)
target_include_directories(sequential_shortest_path PRIVATE "${CMAKE_SOURCE_DIR}/external")
target_link_libraries(sequential_shortest_path PRIVATE multiqueue_internal)
target_compile_features(
  sequential_shortest_path PRIVATE
  cxx_std_17
)

add_executable(sequential_knapsack sequential_knapsack.cpp)
target_include_directories(sequential_knapsack PRIVATE "${CMAKE_SOURCE_DIR}/external")
target_link_libraries(sequential_knapsack PRIVATE multiqueue_internal)
target_compile_features(
  sequential_knapsack PRIVATE
  cxx_std_17
)
