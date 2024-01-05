add_executable(knapsack_generator knapsack_generator.cpp)
target_include_directories(
  knapsack_generator PRIVATE "${CMAKE_SOURCE_DIR}/third_party")

target_link_libraries(knapsack_generator PRIVATE multiqueue_internal)
target_compile_features(
  knapsack_generator PRIVATE
  cxx_std_17
)

add_executable(analyze_quality analyze_quality.cpp)
target_include_directories(
  analyze_quality PRIVATE ".")
target_compile_definitions(analyze_quality
                           PRIVATE $<$<CONFIG:Debug>:-DREPLAY_TREE_DEBUG>)

