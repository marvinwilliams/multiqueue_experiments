include(Catch)

add_executable(unit_tests replay_tree.cpp)
target_link_libraries(unit_tests PRIVATE Threads::Threads)
target_link_libraries_system(unit_tests PRIVATE Catch2::Catch2 Catch2::Catch2WithMain)

if(BUILD_TESTING)
  catch_discover_tests(unit_tests)
endif()
