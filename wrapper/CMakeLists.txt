set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

set(MQ_CONFIG MQ_NUM_POP_PQS=${EXP_MQ_POP_CANDIDATES})
if(EXP_MQ_PQ STREQUAL "heap")
  list(APPEND MQ_CONFIG MQ_HEAP_ARITY=${EXP_MQ_HEAP_ARITY} MQ_INSERTION_BUFFER_SIZE=${EXP_MQ_INSERTION_BUFFER_SIZE}
       MQ_DELETION_BUFFER_SIZE=${EXP_MQ_DELETION_BUFFER_SIZE})
elseif(EXP_MQ_PQ STREQUAL "std")
  list(APPEND MQ_CONFIG MQ_USE_STD_PQ MQ_INSERTION_BUFFER_SIZE=${EXP_MQ_INSERTION_BUFFER_SIZE}
       MQ_DELETION_BUFFER_SIZE=${EXP_MQ_DELETION_BUFFER_SIZE})
elseif(EXP_MQ_PQ STREQUAL "btree")
  list(APPEND MQ_CONFIG MQ_USE_BTREE)
elseif(EXP_MQ_PQ STREQUAL "merge_heap")
  list(APPEND MQ_CONFIG MQ_USE_MERGE_HEAP MQ_MERGE_HEAP_NODE_SIZE=${EXP_MQ_MERGE_HEAP_NODE_SIZE})
else()
  message(FATAL_ERROR "Invalid EXP_MQ_PQ '${EXP_MQ_PQ}' (heap, std, btree, merge_heap)")
endif()
message(STATUS "MultiQueue configuration: ${MQ_CONFIG}")

add_library(mq_base INTERFACE)
target_link_libraries(mq_base INTERFACE multiqueue::multiqueue)
target_compile_definitions(mq_base INTERFACE PQ_MQ)

add_library(mq INTERFACE)
target_link_libraries(mq INTERFACE mq_base)
target_compile_definitions(mq INTERFACE ${MQ_CONFIG})

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

add_library(mq_stick_replace INTERFACE)
target_compile_definitions(mq_stick_replace INTERFACE MQ_MODE_STICK_REPLACE)
target_link_libraries(mq_stick_replace INTERFACE mq)

add_library(mq_random_merge_heap INTERFACE)
target_compile_definitions(
  mq_random_merge_heap INTERFACE MQ_MODE_RANDOM MQ_USE_MERGE_HEAP MQ_MERGE_HEAP_NODE_SIZE=${EXP_MQ_MERGE_HEAP_NODE_SIZE}
                                 MQ_NUM_POP_PQS=${EXP_MQ_POP_CANDIDATES})
target_link_libraries(mq_random_merge_heap INTERFACE mq_base)

list(APPEND MQ_VARIANTS mq_random mq_random_strict mq_stick_random mq_stick_swap
     mq_stick_mark mq_stick_replace mq_random_merge_heap)

if(EXP_TUNING_TARGETS)
  function(add_tuning_target name)
    add_library(${name} INTERFACE)
    target_link_libraries(${name} INTERFACE mq_base)
    target_compile_definitions(${name} INTERFACE MQ_MODE_RANDOM ${ARGN})
    set(TUNING_TARGETS
        ${TUNING_TARGETS} ${name}
        PARENT_SCOPE)
  endfunction()

  foreach(i ${EXP_TUNING_POP_PQS})
    add_tuning_target(mq_tuning_pop_${i} MQ_INSERTION_BUFFER_SIZE=0 MQ_DELETION_BUFFER_SIZE=0 MQ_HEAP_ARITY=2
                      MQ_NUM_POP_PQS=${i})
  endforeach()
  foreach(ibs ${EXP_TUNING_BUFFER_SIZES})
    foreach(dbs ${EXP_TUNING_BUFFER_SIZES})
      if((ibs EQUAL 0) AND NOT (dbs EQUAL 0) OR NOT (ibs EQUAL 0) AND (dbs EQUAL 0))
        continue()
      endif()
      foreach(d ${EXP_TUNING_HEAP_ARITIES})
        add_tuning_target(mq_tuning_ibs_${ibs}_dbs_${dbs}_d_${d} MQ_INSERTION_BUFFER_SIZE=${ibs}
                          MQ_DELETION_BUFFER_SIZE=${dbs} MQ_HEAP_ARITY=${d} MQ_NUM_POP_PQS=2)
      endforeach()
      add_tuning_target(mq_tuning_ibs_${ibs}_dbs_${dbs}_std MQ_INSERTION_BUFFER_SIZE=${ibs}
                        MQ_DELETION_BUFFER_SIZE=${dbs} MQ_NUM_POP_PQS=2 MQ_USE_STD_PQ)
    endforeach()
  endforeach()
  add_tuning_target(mq_tuning_btree MQ_USE_BTREE MQ_NUM_POP_PQS=2)
  foreach(n ${EXP_TUNING_MERGE_HEAP_NODE_SIZES})
    add_tuning_target(mq_tuning_mh_${n} MQ_USE_MERGE_HEAP MQ_MERGE_HEAP_NODE_SIZE=${n} MQ_NUM_POP_PQS=2)
  endforeach()

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

add_library(locked_pq INTERFACE)
target_compile_definitions(locked_pq INTERFACE PQ_LOCKED_PQ)

add_library(mq_pq INTERFACE)
target_link_libraries(mq_pq INTERFACE multiqueue::multiqueue)
target_compile_definitions(mq_pq INTERFACE PQ_MQ_PQ ${MQ_CONFIG})

# GENERIC_COMPETITORS also support non-integer keys (knapsack_float)
set(GENERIC_COMPETITORS
    ${MQ_VARIANTS}
    klsm4
    klsm256
    klsm1024
    klsm4096
    smq
    locked_pq
    mq_pq)
set(COMPETITORS ${GENERIC_COMPETITORS})

find_package(GSL)

if (GSL_FOUND)
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

  list(APPEND COMPETITORS linden capq spraylist)
else()
  message(WARNING "GSL not found, not building linden, capq, spraylist")
endif()

find_package(TBB)

if(TBB_FOUND)
  add_library(tbb_pq INTERFACE)
  target_link_libraries(tbb_pq INTERFACE TBB::tbb)
  target_compile_definitions(tbb_pq INTERFACE PQ_TBB_PQ)

  list(APPEND COMPETITORS tbb_pq)
  list(APPEND GENERIC_COMPETITORS tbb_pq)
else()
  message(WARNING "TBB not found, not building TBB priority queue")
endif()

set(COMPETITORS
    ${COMPETITORS}
    PARENT_SCOPE)
set(GENERIC_COMPETITORS
    ${GENERIC_COMPETITORS}
    PARENT_SCOPE)
