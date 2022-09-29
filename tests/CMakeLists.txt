add_executable(replay_tree_test replay_tree.cpp)
target_link_libraries(replay_tree_test PRIVATE Catch2::Catch2WithMain)
target_include_directories(replay_tree_test PRIVATE "..")

add_executable(wrapper_test wrapper.cpp)
target_link_libraries(wrapper_test PRIVATE Catch2::Catch2WithMain)
target_link_libraries(wrapper_test PRIVATE capq klsm linden spraylist)

if(BUILD_TESTING)
  catch_discover_tests(replay_tree_test)
  catch_discover_tests(wrapper_test)
endif()
