add_executable(knapsack_generator knapsack_generator.cpp)
target_include_directories(
  knapsack_generator PRIVATE "${CMAKE_SOURCE_DIR}/third_party")

target_link_libraries(knapsack_generator PRIVATE multiqueue_internal)
target_compile_features(
  knapsack_generator PRIVATE
  cxx_std_17
)
