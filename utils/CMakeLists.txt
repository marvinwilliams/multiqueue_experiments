add_library(threading STATIC EXCLUDE_FROM_ALL)
target_include_directories(threading PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(threading PUBLIC Threads::Threads)
target_compile_features(
  threading PUBLIC
  cxx_std_17
)
target_sources(threading
  PRIVATE
    threading.cpp
  PUBLIC
    threading.hpp
)

add_library(thread_coordination INTERFACE)
target_link_libraries(thread_coordination INTERFACE threading)
target_sources(thread_coordination
  INTERFACE
    thread_coordination.hpp
)
