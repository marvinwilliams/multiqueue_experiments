set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_library(mq INTERFACE)
# target_include_directories(mq INTERFACE ".")
target_link_libraries(mq INTERFACE multiqueue::multiqueue)
target_compile_definitions(mq INTERFACE PQ_MQ)

add_library(mq_classic INTERFACE)
target_compile_definitions(mq_classic INTERFACE MQ_OP_POLICY=0)
target_link_libraries(mq_classic INTERFACE mq)

add_library(mq_classic_strict INTERFACE)
target_compile_definitions(mq_classic_strict INTERFACE MQ_OP_POLICY=0)
target_link_libraries(mq_classic_strict INTERFACE mq)
target_compile_definitions(mq_classic_strict INTERFACE MQ_NO_POP_STALE)

add_library(mq_stick_random INTERFACE)
target_compile_definitions(mq_stick_random INTERFACE MQ_OP_POLICY=1)
target_link_libraries(mq_stick_random INTERFACE mq)

add_library(mq_stick_random_shared INTERFACE)
target_compile_definitions(mq_stick_random_shared INTERFACE MQ_OP_POLICY=2)
target_link_libraries(mq_stick_random_shared INTERFACE mq)

add_library(mq_stick_swap INTERFACE)
target_compile_definitions(mq_stick_swap INTERFACE MQ_OP_POLICY=3)
target_link_libraries(mq_stick_swap INTERFACE mq)

add_library(mq_stick_parametric INTERFACE)
target_compile_definitions(mq_stick_parametric INTERFACE MQ_OP_POLICY=4)
target_link_libraries(mq_stick_parametric INTERFACE mq)

list(
  APPEND
  MQ_VARIANTS
  mq_classic
  mq_classic_strict
  mq_stick_random
  mq_stick_random_shared
  mq_stick_swap
  mq_stick_parametric)

if(EXP_TUNING_TARGETS)
  add_library(tuning_baseline INTERFACE)
  target_compile_definitions(
    tuning_baseline
    INTERFACE MQ_OP_POLICY=0 MQ_INSERTION_BUFFER_SIZE=0
              MQ_DELETION_BUFFER_SIZE=0 MQ_HEAP_ARITY=2 MQ_NUM_POP_PQS=2)
  target_link_libraries(tuning_baseline INTERFACE mq)

  foreach(i ${EXP_NUM_POP_PQS})
    add_library(mq_tuning_pop_${i} INTERFACE)
    target_link_libraries(mq_tuning_pop_${i} INTERFACE tuning_baseline)
    target_compile_definitions(mq_tuning_pop_${i} INTERFACE MQ_NUM_POP_PQS=${i})
    list(APPEND TUNING_TARGETS mq_tuning_pop_${i})
  endforeach()
  foreach(bs ${EXP_BUFFER_SIZES})
    foreach(d ${EXP_ARITIES})
      add_library(mq_tuning_bs_${bs}_d_${d} INTERFACE)
      target_link_libraries(mq_tuning_bs_${bs}_d_${d} INTERFACE tuning_baseline)
      target_compile_definitions(
        mq_tuning_bs_${bs}_d_${d}
        INTERFACE MQ_INSERTION_BUFFER_SIZE=${bs} MQ_DELETION_BUFFER_SIZE=${bs}
                  MQ_HEAP_ARITY=${d})
      list(APPEND TUNING_TARGETS mq_tuning_bs_${bs}_d_${d})
    endforeach()
    add_library(mq_tuning_bs_${bs}_std INTERFACE)
    target_link_libraries(mq_tuning_bs_${bs}_std INTERFACE tuning_baseline)
    target_compile_definitions(
      mq_tuning_bs_${bs}_std
      INTERFACE MQ_INSERTION_BUFFER_SIZE=${bs} MQ_DELETION_BUFFER_SIZE=${bs}
                MQ_USE_STD_PQ)
    list(APPEND TUNING_TARGETS mq_tuning_bs_${bs}_std)
  endforeach()
  add_library(mq_tuning_btree INTERFACE)
  target_link_libraries(mq_tuning_btree INTERFACE tuning_baseline)
  target_compile_definitions(mq_tuning_btree INTERFACE MQ_USE_BTREE)
  list(APPEND TUNING_TARGETS mq_tuning_btree)

  set(TUNING_TARGETS
      ${TUNING_TARGETS}
      PARENT_SCOPE)
endif()

set(KLSM_DIR "${CMAKE_SOURCE_DIR}/external/klsm")

add_library(klsm ${KLSM_DIR}/src/util/thread_local_ptr.cpp)
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
set_property(TARGET smq PROPERTY pq_name smq)
target_compile_definitions(smq INTERFACE PQ_SMQ)

add_library(
  linden
  ${KLSM_DIR}/lib/spraylist_linden/linden_common.c
  ${KLSM_DIR}/lib/spraylist_linden/linden.c
  ${KLSM_DIR}/lib/spraylist_linden/gc/gc.c
  ${KLSM_DIR}/lib/spraylist_linden/gc/ptst.c)
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
target_include_directories(capq PUBLIC ${KLSM_DIR}/lib ${KLSM_DIR}/lib/capq
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
target_compile_options(spraylist PRIVATE -w -std=gnu17)
target_compile_definitions(
  spraylist PUBLIC INTEL LOCKFREE CACHE_LINE_SIZE=${L1_CACHE_LINESIZE})
target_include_directories(
  spraylist PUBLIC ${KLSM_DIR}/lib ${KLSM_DIR}/lib/spraylist_linden/atomic_ops
                   ${KLSM_DIR}/lib/spraylist_linden/include)
set_property(TARGET spraylist PROPERTY pq_name spraylist)
target_compile_definitions(spraylist INTERFACE PQ_SPRAYLIST)

find_package(TBB)

if (TBB_FOUND)
  add_library(tbb_pq INTERFACE)
  target_link_libraries(tbb_pq INTERFACE tbb)
  target_compile_definitions(tbb_pq INTERFACE PQ_TBB_PQ)
else()
  message(WARNING "TBB not found, not building TBB priority queue")
endif()

add_library(locked_pq INTERFACE)
target_compile_definitions(locked_pq INTERFACE PQ_LOCKED_PQ)

set(COMPETITORS
    ${MQ_VARIANTS}
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
    ${MQ_VARIANTS} smq tbb_pq locked_pq
    PARENT_SCOPE)
