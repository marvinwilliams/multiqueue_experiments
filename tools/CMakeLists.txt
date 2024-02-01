add_executable(knapsack_generator knapsack_generator.cpp)
target_include_directories(knapsack_generator
                           SYSTEM PRIVATE "${CMAKE_SOURCE_DIR}/third_party")

target_compile_features(knapsack_generator PRIVATE cxx_std_17)

add_executable(analyze_quality analyze_quality.cpp)
target_include_directories(analyze_quality PRIVATE ".")
