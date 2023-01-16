set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_library(wrapper INTERFACE)
target_include_directories(wrapper INTERFACE "..")

add_library(
  klsm STATIC
  ${EXTERNAL_DIR}/klsm/src/util/thread_local_ptr.cpp
)
target_link_libraries(klsm INTERFACE wrapper)
target_include_directories(
  klsm SYSTEM INTERFACE
  ${EXTERNAL_DIR}/klsm/src
  ${EXTERNAL_DIR}/klsm/src/util
)
add_library(klsm256 INTERFACE)
target_link_libraries(klsm256 INTERFACE klsm)
set_property(TARGET klsm256 PROPERTY pq_name klsm256)
target_compile_definitions(klsm256 INTERFACE PQ_KLSM256)
add_library(klsm1024 INTERFACE)
target_link_libraries(klsm1024 INTERFACE klsm)
set_property(TARGET klsm1024 PROPERTY pq_name klsm1024)
target_compile_definitions(klsm1024 INTERFACE PQ_KLSM1024)
add_library(klsm4096 INTERFACE)
target_link_libraries(klsm4096 INTERFACE klsm)
set_property(TARGET klsm4096 PROPERTY pq_name klsm4096)
target_compile_definitions(klsm4096 INTERFACE PQ_KLSM4096)

add_library(
  linden STATIC
  ${EXTERNAL_DIR}/klsm/lib/spraylist_linden/linden_common.c
  ${EXTERNAL_DIR}/klsm/lib/spraylist_linden/linden.c
  ${EXTERNAL_DIR}/klsm/lib/spraylist_linden/gc/gc.c
  ${EXTERNAL_DIR}/klsm/lib/spraylist_linden/gc/ptst.c
  linden.cpp
)
target_link_libraries(linden PUBLIC wrapper)
target_compile_options(linden PRIVATE -w -std=gnu++11)
target_compile_definitions(linden PRIVATE INTEL CACHE_LINE_SIZE=${L1_CACHE_LINESIZE})
target_include_directories(
  linden PRIVATE
  ${EXTERNAL_DIR}/klsm/lib
  ${EXTERNAL_DIR}/klsm/lib/spraylist_linden/atomic_ops
  ${EXTERNAL_DIR}/klsm/lib/spraylist_linden/include
)
set_property(TARGET linden PROPERTY pq_name linden)
target_compile_definitions(linden INTERFACE PQ_LINDEN)

add_library(
  capq STATIC
  ${EXTERNAL_DIR}/klsm/lib/capq/capq.c
  ${EXTERNAL_DIR}/klsm/lib/capq/fat_skiplist.c
  ${EXTERNAL_DIR}/klsm/lib/capq/qdlocksrc/locks/qd_lock.c
  ${EXTERNAL_DIR}/klsm/lib/capq/qdlocksrc/locks/tatas_lock.c
  ${EXTERNAL_DIR}/klsm/lib/capq/gc/gc.c
  ${EXTERNAL_DIR}/klsm/lib/capq/gc/ptst.c
  capq.cpp
)
target_link_libraries(capq PUBLIC wrapper linden)
target_include_directories(
  capq PRIVATE
  ${EXTERNAL_DIR}/klsm/lib
  ${EXTERNAL_DIR}/klsm/lib/capq
  ${EXTERNAL_DIR}/klsm/lib/capq/qdlocksrc
)
target_compile_options(capq PRIVATE -w -std=gnu++11)
target_compile_definitions(capq PRIVATE INTEL CACHE_LINE_SIZE=${L1_CACHE_LINESIZE})
target_compile_definitions(capq INTERFACE PQ_CAPQ)

add_library(
  spraylist STATIC
  ${EXTERNAL_DIR}/klsm/lib/spraylist_linden/fraser.c
  ${EXTERNAL_DIR}/klsm/lib/spraylist_linden/intset.c
  ${EXTERNAL_DIR}/klsm/lib/spraylist_linden/linden_common.c
  ${EXTERNAL_DIR}/klsm/lib/spraylist_linden/linden.c
  ${EXTERNAL_DIR}/klsm/lib/spraylist_linden/measurements.c
  ${EXTERNAL_DIR}/klsm/lib/spraylist_linden/pqueue.c
  ${EXTERNAL_DIR}/klsm/lib/spraylist_linden/skiplist.c
  ${EXTERNAL_DIR}/klsm/lib/spraylist_linden/ssalloc.c
  spraylist.cpp
)
target_link_libraries(spraylist PUBLIC wrapper)
target_compile_options(spraylist PRIVATE -w -std=gnu11)
target_compile_definitions(spraylist PRIVATE INTEL PRIVATE LOCKFREE PRIVATE CACHE_LINE_SIZE=${L1_CACHE_LINESIZE} PRIVATE SSALLOC_USE_MALLOC)
target_include_directories(
  spraylist PRIVATE
  ${EXTERNAL_DIR}/klsm/lib
  ${EXTERNAL_DIR}/klsm/lib/spraylist_linden/atomic_ops
  ${EXTERNAL_DIR}/klsm/lib/spraylist_linden/include
)
set_property(TARGET spraylist PROPERTY pq_name spraylist)
target_compile_definitions(spraylist INTERFACE PQ_SPRAYLIST)

add_library(
  tbb_q INTERFACE
  )
target_link_libraries(tbb_q INTERFACE wrapper tbb)
target_compile_definitions(tbb_q INTERFACE PQ_TBB_Q)

add_library(
  tbb_pq INTERFACE
  )
target_link_libraries(tbb_pq INTERFACE wrapper tbb)
target_compile_definitions(tbb_pq INTERFACE PQ_TBB_PQ)
