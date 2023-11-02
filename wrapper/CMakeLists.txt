set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

set(KLSM_DIR "${CMAKE_SOURCE_DIR}/third_party/klsm")

add_library(wrapper INTERFACE)
target_include_directories(wrapper SYSTEM INTERFACE ${CMAKE_SOURCE_DIR}
                                             ${CMAKE_SOURCE_DIR}/third_party)

add_library(mq INTERFACE)
target_link_libraries(mq INTERFACE multiqueue::multiqueue wrapper)
target_compile_definitions(mq INTERFACE PQ_MQ)
if(MULTIQUEUE_EXP_COMPARE_STRICT)
  target_compile_definitions(mq INTERFACE MQ_COMPARE_STRICT)
endif()
if(MULTIQUEUE_EXP_COUNT_STATS)
  target_compile_definitions(mq INTERFACE MQ_COUNT_STATS)
endif()
if(MULTIQUEUE_EXP_NO_SCAN_ON_FAILED_POP)
  target_compile_definitions(mq INTERFACE MQ_NO_SCAN_ON_FAILED_POP)
endif()

add_library(klsm ${KLSM_DIR}/src/util/thread_local_ptr.cpp)
target_link_libraries(klsm INTERFACE wrapper)
target_include_directories(klsm SYSTEM INTERFACE ${KLSM_DIR}/src
                                                 ${KLSM_DIR}/src/util)
add_library(klsm4 INTERFACE)
target_link_libraries(klsm4 INTERFACE klsm)
set_property(TARGET klsm4 PROPERTY pq_name klsm4)
target_compile_definitions(
  klsm4
  INTERFACE PQ_KLSM
  INTERFACE KLSM_K=4)
add_library(klsm256 INTERFACE)
target_link_libraries(klsm256 INTERFACE klsm)
set_property(TARGET klsm256 PROPERTY pq_name klsm256)
target_compile_definitions(
  klsm256
  INTERFACE PQ_KLSM
  INTERFACE KLSM_K=256)
add_library(klsm1024 INTERFACE)
target_link_libraries(klsm1024 INTERFACE klsm)
set_property(TARGET klsm1024 PROPERTY pq_name klsm1024)
target_compile_definitions(
  klsm1024
  INTERFACE PQ_KLSM
  INTERFACE KLSM_K=1024)
add_library(klsm4096 INTERFACE)
target_link_libraries(klsm4096 INTERFACE klsm)
set_property(TARGET klsm4096 PROPERTY pq_name klsm4096)
target_compile_definitions(
  klsm4096
  INTERFACE PQ_KLSM
  INTERFACE KLSM_K=4096)

add_library(smq INTERFACE)
target_link_libraries(smq INTERFACE wrapper)
set_property(TARGET smq PROPERTY pq_name smq)
target_compile_definitions(smq INTERFACE PQ_SMQ)

add_library(
  linden
  ${KLSM_DIR}/lib/spraylist_linden/linden_common.c
  ${KLSM_DIR}/lib/spraylist_linden/linden.c
  ${KLSM_DIR}/lib/spraylist_linden/gc/gc.c
  ${KLSM_DIR}/lib/spraylist_linden/gc/ptst.c
  )
target_link_libraries(linden INTERFACE wrapper)
target_compile_options(linden PRIVATE -w -std=gnu++17)
target_compile_definitions(linden PUBLIC INTEL
                                          CACHE_LINE_SIZE=${L1_CACHE_LINESIZE})
target_include_directories(
  linden PRIVATE ${KLSM_DIR}/lib ${KLSM_DIR}/lib/spraylist_linden/atomic_ops
                 ${KLSM_DIR}/lib/spraylist_linden/include)
set_property(TARGET linden PROPERTY pq_name linden)
target_compile_definitions(linden INTERFACE PQ_LINDEN)

add_library(
  capq STATIC
  ${KLSM_DIR}/lib/capq/capq.c
  ${KLSM_DIR}/lib/capq/fat_skiplist.c
  ${KLSM_DIR}/lib/capq/qdlocksrc/locks/qd_lock.c
  ${KLSM_DIR}/lib/capq/qdlocksrc/locks/tatas_lock.c
  ${KLSM_DIR}/lib/capq/gc/gc.c
  ${KLSM_DIR}/lib/capq/gc/ptst.c
  )
target_link_libraries(capq INTERFACE wrapper)
target_include_directories(
  capq PUBLIC ${KLSM_DIR}/lib ${EXTERNAL_DIR}/klsm/lib/capq
               ${KLSM_DIR}/lib/capq/qdlocksrc)
target_compile_options(capq PRIVATE -w -std=gnu++17)
target_compile_definitions(capq PUBLIC INTEL
                                        CACHE_LINE_SIZE=${L1_CACHE_LINESIZE})
target_compile_definitions(capq INTERFACE PQ_CAPQ)

add_library(
  spraylist
  ${KLSM_DIR}/lib/spraylist_linden/fraser.c
  ${KLSM_DIR}/lib/spraylist_linden/intset.c
  ${KLSM_DIR}/lib/spraylist_linden/linden_common.c
  ${KLSM_DIR}/lib/spraylist_linden/linden.c
  ${KLSM_DIR}/lib/spraylist_linden/measurements.c
  ${KLSM_DIR}/lib/spraylist_linden/pqueue.c
  ${KLSM_DIR}/lib/spraylist_linden/skiplist.c
  ${KLSM_DIR}/lib/spraylist_linden/ssalloc.c
  )
target_link_libraries(spraylist INTERFACE wrapper)
target_compile_options(spraylist PRIVATE -w -std=gnu17)
target_compile_definitions(spraylist PUBLIC INTEL LOCKFREE CACHE_LINE_SIZE=${L1_CACHE_LINESIZE})
target_include_directories(
  spraylist PUBLIC ${KLSM_DIR}/lib ${KLSM_DIR}/lib/spraylist_linden/atomic_ops
                    ${KLSM_DIR}/lib/spraylist_linden/include)
set_property(TARGET spraylist PROPERTY pq_name spraylist)
target_compile_definitions(spraylist INTERFACE PQ_SPRAYLIST)

add_library(tbb_pq INTERFACE)
target_link_libraries(tbb_pq INTERFACE wrapper tbb)
target_compile_definitions(tbb_pq INTERFACE PQ_TBB_PQ)
