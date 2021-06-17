#ifndef WRAPPER_SPRAYLIST_HPP_INCLUDED
#define WRAPPER_SPRAYLIST_HPP_INCLUDED

// Adapted from klsm

extern "C" {
#include "spraylist_linden/include/random.h"
#include "spraylist_linden/intset.h"
#include "spraylist_linden/linden.h"
#include "spraylist_linden/pqueue.h"
#include "ssalloc.h"
}

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>

struct sl_intset;

/** See documentation of --elasticity in spraylist/test.c. */
#define READ_ADD_REM_ELASTIC_TX (4)

namespace wrapper {

constexpr unsigned int INITIAL_SIZE = 1000000;
extern thread_local thread_data_t* d;

template <typename Key, typename T>
class spraylist {
    using pq_t = sl_intset;

    pq_t* pq_;

   public:
    struct Handle {};
    explicit spraylist(size_t const num_threads) {
        init_thread(num_threads);
        *levelmax = floor_log_2(INITIAL_SIZE);
        pq_ = sl_set_new();
    }
    virtual ~spraylist() {
        sl_set_delete(pq_);
        delete d;
    }

    constexpr Handle get_handle(unsigned int) {
        return Handle{};
    }

    void init_thread(size_t const num_threads) {
        static thread_local bool initialized = false;
        static thread_local thread_data_t* d;
        if (!initialized) {
            ssalloc_init(num_threads);
            seeds = seed_rand();

            d = new thread_data_t;
            d->seed = rand();
            d->seed2 = rand();

            initialized = true;
        }

        d->nb_threads = num_threads;
    }

    void push(Handle, std::pair<uint32_t, uint32_t> const& value) {
        sl_add_val(pq_, value.first, value.second, TRANSACTIONAL);
    }
    bool extract_top(Handle, std::pair<uint32_t, uint32_t>& retval) {
        slkey_t k_ret;
        val_t v_ret;
        int ret;
        do {
            ret = spray_delete_min_key(pq_, &k_ret, &v_ret, d);
        } while (ret == 0 && k_ret != -1);
        retval.first = static_cast<uint32_t>(k_ret);
        retval.second = static_cast<uint32_t>(v_ret);
        return k_ret != -1;
    }

    static std::string description() {
        return "spraylist";
    }
};

}  // namespace wrapper

#endif
