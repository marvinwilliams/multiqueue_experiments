set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

set(KLSM_DIR "${CMAKE_SOURCE_DIR}/third_party/klsm")

add_library(wrapper INTERFACE)
target_include_directories(
  wrapper SYSTEM INTERFACE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/third_party)

add_library(mq INTERFACE)
target_link_libraries(mq INTERFACE multiqueue::multiqueue wrapper)
target_compile_definitions(mq INTERFACE PQ_MQ)
if(MULTIQUEUE_EXP_NO_POP_STALE)
  target_compile_definitions(mq INTERFACE MQ_NO_POP_STALE)
endif()
# if(MULTIQUEUE_EXP_COUNT_STATS) target_compile_definitions(mq INTERFACE
# MQ_COUNT_STATS) endif()
if(MULTIQUEUE_EXP_NO_SCAN_IF_EMPTY)
  target_compile_definitions(mq INTERFACE MQ_NO_SCAN_IF_EMPTY)
endif()

add_library(mq_r INTERFACE)
target_compile_definitions(mq_r INTERFACE MQ_OP_POLICY=0)
target_link_libraries(mq_r INTERFACE mq)

add_library(mq_sri INTERFACE)
target_compile_definitions(mq_sri INTERFACE MQ_OP_POLICY=1)
target_link_libraries(mq_sri INTERFACE mq)

add_library(mq_src INTERFACE)
target_compile_definitions(mq_src INTERFACE MQ_OP_POLICY=2)
target_link_libraries(mq_src INTERFACE mq)

add_library(mq_spi INTERFACE)
target_compile_definitions(mq_spi INTERFACE MQ_OP_POLICY=3)
target_link_libraries(mq_spi INTERFACE mq)

add_library(mq_spg INTERFACE)
target_compile_definitions(mq_spg INTERFACE MQ_OP_POLICY=4)
target_link_libraries(mq_spg INTERFACE mq)

list(
  APPEND
  MQ_VARIANTS
  mq_r
  mq_sri
  mq_src
  mq_spi
  mq_spg)

set(MQ_VARIANTS
    ${MQ_VARIANTS}
    PARENT_SCOPE)

function(add_pop_pq_variant num_pop_pqs)
  set(name mq_pop_pqs_${num_pop_pqs})
  add_library(${name} INTERFACE)
  target_compile_definitions(
    ${name}
    INTERFACE MQ_BUFFERED_PQ_INSERTION_BUFFER_SIZE=0
              MQ_BUFFERED_PQ_DELETION_BUFFER_SIZE=0 MQ_HEAP_ARITY=2
              MQ_NUM_POP_PQS=${num_pop_pqs})
  target_link_libraries(${name} INTERFACE mq_r)
  list(APPEND MQ_TUNING_VARIANTS ${name})
  set(MQ_TUNING_VARIANTS
      ${MQ_TUNING_VARIANTS}
      PARENT_SCOPE)
endfunction()

foreach(num_pop_pqs 1 2 4 8)
  add_pop_pq_variant(${num_pop_pqs})
endforeach()

function(add_tuning_variant buffer_size degree)
  set(name mq_tuning_bs_${buffer_size}_d_${degree})
  add_library(${name} INTERFACE)
  target_compile_definitions(
    ${name}
    INTERFACE MQ_BUFFERED_PQ_INSERTION_BUFFER_SIZE=${buffer_size}
              MQ_BUFFERED_PQ_DELETION_BUFFER_SIZE=${buffer_size}
              MQ_HEAP_ARITY=${degree})
  target_link_libraries(${name} INTERFACE mq_r)
  list(APPEND MQ_TUNING_VARIANTS ${name})
  set(MQ_TUNING_VARIANTS
      ${MQ_TUNING_VARIANTS}
      PARENT_SCOPE)
endfunction()

function(add_tuning_std_pq_variant buffer_size)
  set(name mq_tuning_std_bs_${buffer_size})
  add_library(${name} INTERFACE)
  target_compile_definitions(
    ${name}
    INTERFACE MQ_BUFFERED_PQ_INSERTION_BUFFER_SIZE=${buffer_size}
              MQ_BUFFERED_PQ_DELETION_BUFFER_SIZE=${buffer_size} MQ_USE_STD_PQ)
  target_link_libraries(${name} INTERFACE mq_r)
  list(APPEND MQ_TUNING_VARIANTS ${name})
  set(MQ_TUNING_VARIANTS
      ${MQ_TUNING_VARIANTS}
      PARENT_SCOPE)
endfunction()

foreach(
  buffer_size
  0
  2
  4
  8
  16
  32
  64
  128
  256
  512
  1024)
  foreach(degree 2 4 8 16 32 64)
    add_tuning_variant(${buffer_size} ${degree})
  endforeach()
  add_tuning_std_pq_variant(${buffer_size})
endforeach()

add_library(mq_tuning_btree INTERFACE)
target_compile_definitions(mq_tuning_btree INTERFACE MQ_USE_BTREE)
target_link_libraries(mq_tuning_btree INTERFACE mq_r)
list(APPEND MQ_TUNING_VARIANTS mq_tuning_btree)

set(MQ_TUNING_VARIANTS
    ${MQ_TUNING_VARIANTS}
    PARENT_SCOPE)

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
  ${KLSM_DIR}/lib/spraylist_linden/gc/ptst.c)
target_link_libraries(linden INTERFACE wrapper)
target_compile_options(linden PRIVATE -w -std=gnu++17)
target_compile_definitions(linden PUBLIC INTEL
                                         CACHE_LINE_SIZE=${L1_CACHE_LINESIZE})
target_include_directories(
  linden PUBLIC ${KLSM_DIR}/lib ${KLSM_DIR}/lib/spraylist_linden/atomic_ops
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
  ${KLSM_DIR}/lib/capq/gc/ptst.c)
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
  ${KLSM_DIR}/lib/spraylist_linden/ssalloc.c)
target_link_libraries(spraylist INTERFACE wrapper)
target_compile_options(spraylist PRIVATE -w -std=gnu17)
target_compile_definitions(
  spraylist PUBLIC INTEL LOCKFREE CACHE_LINE_SIZE=${L1_CACHE_LINESIZE})
target_include_directories(
  spraylist PUBLIC ${KLSM_DIR}/lib ${KLSM_DIR}/lib/spraylist_linden/atomic_ops
                   ${KLSM_DIR}/lib/spraylist_linden/include)
set_property(TARGET spraylist PROPERTY pq_name spraylist)
target_compile_definitions(spraylist INTERFACE PQ_SPRAYLIST)

add_library(tbb_pq INTERFACE)
target_link_libraries(tbb_pq INTERFACE wrapper tbb)
target_compile_definitions(tbb_pq INTERFACE PQ_TBB_PQ)

add_library(locked_pq INTERFACE)
target_link_libraries(locked_pq INTERFACE wrapper locked_pq)
target_compile_definitions(locked_pq INTERFACE PQ_LOCKED_PQ)

set(COMPETITORS
    klsm4
    klsm256
    klsm1024
    klsm4096
    smq
    linden
    capq
    spraylist
    tbb_pq
    locked_pq
    PARENT_SCOPE)

set(GENERIC_COMPETITORS
    smq tbb_pq locked_pq
    PARENT_SCOPE)
