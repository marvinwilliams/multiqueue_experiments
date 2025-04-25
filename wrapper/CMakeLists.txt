set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

find_package(GSL)

add_library(mq INTERFACE)
target_link_libraries(mq INTERFACE multiqueue::multiqueue)
target_compile_definitions(mq INTERFACE PQ_MQ)

add_library(mq_random INTERFACE)
target_compile_definitions(mq_random INTERFACE MQ_MODE_RANDOM)
target_link_libraries(mq_random INTERFACE mq)

add_library(mq_random_strict INTERFACE)
target_compile_definitions(mq_random_strict INTERFACE MQ_MODE_RANDOM_STRICT)
target_link_libraries(mq_random_strict INTERFACE mq)

add_library(mq_stick_random INTERFACE)
target_compile_definitions(mq_stick_random INTERFACE MQ_MODE_STICK_RANDOM)
target_link_libraries(mq_stick_random INTERFACE mq)

add_library(mq_stick_swap INTERFACE)
target_compile_definitions(mq_stick_swap INTERFACE MQ_MODE_STICK_SWAP)
target_link_libraries(mq_stick_swap INTERFACE mq)

add_library(mq_stick_mark INTERFACE)
target_compile_definitions(mq_stick_mark INTERFACE MQ_MODE_STICK_MARK)
target_link_libraries(mq_stick_mark INTERFACE mq)

add_library(mq_stick_parametric INTERFACE)
target_compile_definitions(mq_stick_parametric
                           INTERFACE MQ_MODE_STICK_PARAMETRIC)
target_link_libraries(mq_stick_parametric INTERFACE mq)

list(APPEND MQ_VARIANTS mq_random mq_random_strict mq_stick_random
     mq_stick_swap mq_stick_mark)

if(EXP_TUNING_TARGETS)
  foreach(i ${EXP_TUNING_POP_PQS})
    add_library(mq_tuning_pop_${i} INTERFACE)
    target_link_libraries(mq_tuning_pop_${i} INTERFACE mq)
    target_compile_definitions(
      mq_tuning_pop_${i}
      INTERFACE MQ_MODE_RANDOM MQ_INSERTION_BUFFER_SIZE=0
                MQ_DELETION_BUFFER_SIZE=0 MQ_HEAP_ARITY=2 MQ_NUM_POP_PQS=${i})
    list(APPEND TUNING_TARGETS mq_tuning_pop_${i})
  endforeach()
  foreach(ibs ${EXP_TUNING_BUFFER_SIZES})
    foreach(dbs ${EXP_TUNING_BUFFER_SIZES})
      foreach(d ${EXP_TUNING_HEAP_ARITIES})
        add_library(mq_tuning_ibs_${ibs}_dbs_${dbs}_d_${d} INTERFACE)
        target_link_libraries(mq_tuning_ibs_${ibs}_dbs_${dbs}_d_${d}
                              INTERFACE mq)
        target_compile_definitions(
          mq_tuning_ibs_${ibs}_dbs_${dbs}_d_${d}
          INTERFACE MQ_MODE_RANDOM MQ_INSERTION_BUFFER_SIZE=${ibs}
                    MQ_DELETION_BUFFER_SIZE=${dbs} MQ_HEAP_ARITY=${d}
                    MQ_NUM_POP_PQS=2)
        list(APPEND TUNING_TARGETS mq_tuning_ibs_${ibs}_dbs_${dbs}_d_${d})
      endforeach()
      add_library(mq_tuning_ibs_${ibs}_dbs_${dbs}_std INTERFACE)
      target_link_libraries(mq_tuning_ibs_${ibs}_dbs_${dbs}_std INTERFACE mq)
      target_compile_definitions(
        mq_tuning_ibs_${ibs}_dbs_${dbs}_std
        INTERFACE MQ_MODE_RANDOM MQ_INSERTION_BUFFER_SIZE=${ibs}
                  MQ_DELETION_BUFFER_SIZE=${dbs} MQ_HEAP_ARITY=2
                  MQ_NUM_POP_PQS=2 MQ_USE_STD_PQ)
      list(APPEND TUNING_TARGETS mq_tuning_ibs_${ibs}_dbs_${dbs}_std)
    endforeach()
  endforeach()
  add_library(mq_tuning_btree INTERFACE)
  target_link_libraries(mq_tuning_btree INTERFACE mq)
  target_compile_definitions(
    mq_tuning_btree
    INTERFACE MQ_USE_BTREE MQ_MODE_RANDOM MQ_NUM_POP_PQS=2)
  list(APPEND TUNING_TARGETS mq_tuning_btree)

  set(TUNING_TARGETS
      ${TUNING_TARGETS}
      PARENT_SCOPE)
endif()

set(KLSM_DIR "${CMAKE_SOURCE_DIR}/third_party/klsm")

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
                                         CACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} _GNU_SOURCE)
target_include_directories(
  linden SYSTEM
  PUBLIC ${KLSM_DIR}/lib ${KLSM_DIR}/lib/spraylist_linden/atomic_ops
         ${KLSM_DIR}/lib/spraylist_linden/include)
set_property(TARGET linden PROPERTY pq_name linden)
target_compile_definitions(linden INTERFACE PQ_LINDEN)
target_link_libraries(linden PUBLIC GSL::gsl)

add_library(
  capq STATIC
  ${KLSM_DIR}/lib/capq/capq.c
  ${KLSM_DIR}/lib/capq/fat_skiplist.c
  ${KLSM_DIR}/lib/capq/qdlocksrc/locks/qd_lock.c
  ${KLSM_DIR}/lib/capq/qdlocksrc/locks/tatas_lock.c
  ${KLSM_DIR}/lib/capq/gc/gc.c
  ${KLSM_DIR}/lib/capq/gc/ptst.c)
target_include_directories(
  capq SYSTEM PUBLIC ${KLSM_DIR}/lib ${KLSM_DIR}/lib/capq
                     ${KLSM_DIR}/lib/capq/qdlocksrc)
target_compile_options(capq PRIVATE -w -std=gnu11)
target_compile_definitions(capq PUBLIC INTEL
                                        CACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE})
target_compile_definitions(capq INTERFACE PQ_CAPQ)
target_link_libraries(capq PUBLIC GSL::gsl)

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
  spraylist PUBLIC INTEL LOCKFREE _GNU_SOURCE CACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE})
target_include_directories(
  spraylist SYSTEM
  PUBLIC ${KLSM_DIR}/lib ${KLSM_DIR}/lib/spraylist_linden/atomic_ops
         ${KLSM_DIR}/lib/spraylist_linden/include)
set_property(TARGET spraylist PROPERTY pq_name spraylist)
target_compile_definitions(spraylist INTERFACE PQ_SPRAYLIST)
target_link_libraries(spraylist PUBLIC GSL::gsl)

find_package(TBB)

if(TBB_FOUND)
  add_library(tbb_pq INTERFACE)
  target_link_libraries(tbb_pq INTERFACE tbb)
  target_compile_definitions(tbb_pq INTERFACE PQ_TBB_PQ)
else()
  message(WARNING "TBB not found, not building TBB priority queue")
endif()

add_library(locked_pq INTERFACE)
target_compile_definitions(locked_pq INTERFACE PQ_LOCKED_PQ)

add_library(mq_pq INTERFACE)
target_link_libraries(mq_pq INTERFACE multiqueue::multiqueue)
target_compile_definitions(mq_pq INTERFACE PQ_MQ_PQ)

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
    mq_pq
    PARENT_SCOPE)

set(GENERIC_COMPETITORS
    ${MQ_VARIANTS}
    klsm4
    klsm256
    klsm1024
    klsm4096
    smq
    tbb_pq
    locked_pq
    mq_pq
    PARENT_SCOPE)
