# The multiqueue subproject vendors its own Catch2 and adds it when the target
# does not exist yet, so only look for an installed one if that did not happen.
if(NOT TARGET Catch2::Catch2WithMain)
  find_package(Catch2 3 REQUIRED)
endif()
# catch_discover_tests lives in Catch2's extras: next to the vendored sources
# when it came from add_subdirectory, on CMAKE_MODULE_PATH when it came from
# find_package.
if(DEFINED Catch2_SOURCE_DIR)
  list(APPEND CMAKE_MODULE_PATH "${Catch2_SOURCE_DIR}/extras")
endif()
include(Catch)

add_executable(replay_tree_test replay_tree.cpp)
target_link_libraries(replay_tree_test PRIVATE Catch2::Catch2WithMain)
target_include_directories(replay_tree_test PRIVATE "..")
target_compile_features(replay_tree_test PRIVATE cxx_std_17)

if(BUILD_TESTING)
  catch_discover_tests(replay_tree_test)
endif()
