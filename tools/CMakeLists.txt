add_executable(knapsack_generator knapsack_generator.cpp)
target_include_directories(knapsack_generator PRIVATE "${EXTERNAL_DIR}")
target_link_libraries(knapsack_generator PRIVATE multiqueue_internal)
target_compile_features(
  knapsack_generator PRIVATE
  cxx_std_17
)

add_executable(sequential_sssp sequential_shortest_path.cpp)
target_include_directories(sequential_sssp PRIVATE "${EXTERNAL_DIR}")
target_link_libraries(sequential_sssp PRIVATE multiqueue_internal)
target_compile_features(
  sequential_sssp PRIVATE
  cxx_std_17
)

add_executable(sequential_knapsack sequential_knapsack.cpp)
target_include_directories(sequential_knapsack PRIVATE "${EXTERNAL_DIR}")
target_link_libraries(sequential_knapsack PRIVATE multiqueue_internal)
target_compile_features(
  sequential_knapsack PRIVATE
  cxx_std_17
)
