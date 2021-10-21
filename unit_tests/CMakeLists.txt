include(Catch)

add_executable(replay_tree_test replay_tree.cpp)
target_link_libraries_system(replay_tree_test PRIVATE Catch2::Catch2 Catch2::Catch2WithMain)
target_include_directories(replay_tree_test PRIVATE "..")

if(BUILD_TESTING)
  catch_discover_tests(replay_tree_test)
endif()

add_executable(wrapper_test wrapper.cpp)
target_link_libraries_system(wrapper_test PRIVATE Catch2::Catch2 Catch2::Catch2WithMain)
target_link_libraries(wrapper_test PRIVATE capq_wrapper klsm_wrapper linden_wrapper spraylist_wrapper)

if(BUILD_TESTING)
  catch_discover_tests(wrapper_test)
endif()

