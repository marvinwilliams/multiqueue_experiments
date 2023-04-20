add_library(threading STATIC EXCLUDE_FROM_ALL threading.cpp)
target_compile_features(threading PUBLIC cxx_std_17)
target_include_directories(threading PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(threading PUBLIC Threads::Threads)

add_library(thread_coordination INTERFACE)
target_link_libraries(thread_coordination INTERFACE threading)
target_sources(thread_coordination INTERFACE thread_coordination.hpp)
