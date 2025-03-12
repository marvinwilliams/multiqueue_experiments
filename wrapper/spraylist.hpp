#pragma once

// Adapted from klsm

extern "C" {
#include "spraylist_linden/include/random.h"
#include "spraylist_linden/intset.h"
#include "ssalloc.h"
}

#undef min
#undef max

#include "util.hpp"

#include <cxxopts.hpp>

#include <cstdint>
#include <limits>
#include <mutex>
#include <optional>
#include <ostream>
#include <utility>

inline __thread unsigned long* seeds;  // NOLINT

namespace wrapper::spraylist {

template <bool Min, typename Key = unsigned long, typename T = Key>
class Spraylist {
    static_assert(std::is_same_v<Key, unsigned long> && std::is_same_v<T, unsigned long>,
                  "Spraylist only supports unsigned long as key and value type");

   public:
    using key_type = unsigned long;
    using mapped_type = unsigned long;
    using value_type = std::pair<key_type, mapped_type>;

   private:
    // No unique_ptr because of alignment
    ::sl_intset_t* pq_;
    int num_threads_;

    static constexpr key_type sentinel = std::numeric_limits<key_type>::max();

    class Handle {
        friend Spraylist;
        ::sl_intset_t* pq_;
        // No unique_ptr because of alignment
        ::thread_data_t* data_;

        explicit Handle(::sl_intset_t& pq, ::thread_data_t* data) : pq_{&pq}, data_{std::move(data)} {
        }

       public:
        ~Handle() {
            delete data_;
        }

        Handle(Handle const&) = delete;
        Handle& operator=(Handle const&) = delete;
        Handle(Handle&& other) noexcept {
            pq_ = std::exchange(other.pq_, nullptr);
            data_ = std::exchange(other.data_, nullptr);
        }
        Handle& operator=(Handle&& other) noexcept {
            pq_ = std::exchange(other.pq_, nullptr);
            data_ = std::exchange(other.data_, nullptr);
            return *this;
        }

        void push(value_type const& value) {
            ::sl_add_val(pq_, Min ? value.first : sentinel - value.first - 1, value.second, TRANSACTIONAL);
        }

        /* std::optional<value_type> try_pop() { */
        /*     unsigned long key{0}; */
        /*     unsigned long value; */
        /*     int ret; */
        /*     do { */
        /*         ret = ::spray_delete_min_key(pq_, &key, &value, data_); */
        /*     } while (ret == 0 && key != sentinel); */
        /*     if (key == sentinel) { */
        /*         return std::nullopt; */
        /*     } */
        /*     return value_type{Min ? key : sentinel - key - 1, value}; */
        /* }; */
        std::optional<value_type> try_pop() {
            while (true) {
                unsigned long key{0};
                unsigned long value;
                int ret = ::spray_delete_min_key(pq_, &key, &value, data_);
                if (key == sentinel) {
                    return std::nullopt;
                }
                if (ret != 0) {
                    return std::make_optional<value_type>(Min ? key : sentinel - key - 1, value);
                }
            }
        }
    };

    [[nodiscard]] ::thread_data_t* new_thread_data() const {
        static std::mutex m;
        auto l = std::scoped_lock(m);
        ::ssalloc_init(num_threads_);
        seeds = ::seed_rand();
        auto thread_data = new thread_data_t;
        thread_data->seed = static_cast<unsigned int>(::rand());
        thread_data->seed2 = static_cast<unsigned int>(::rand());
        thread_data->nb_threads = static_cast<unsigned int>(num_threads_);
        return thread_data;
    }

    static unsigned int floor_log_2(unsigned long long n) {
        unsigned int ret = 0;
        while (n >>= 1 > 0) {
            ++ret;
        }
        return ret;
    }

   public:
    using handle_type = Handle;
    using settings_type = util::EmptySettings;

    Spraylist(int num_threads, std::size_t initial_capacity, settings_type const& /*unused*/)
        : num_threads_(num_threads) {
        ::ssalloc_init(num_threads_);
        *::levelmax = static_cast<std::uint8_t>(floor_log_2(initial_capacity));
        pq_ = sl_set_new();
    }
    Spraylist(Spraylist const&) = delete;
    Spraylist& operator=(Spraylist const&) = delete;
    Spraylist(Spraylist&& other) noexcept {
        pq_ = std::exchange(other.pq_, nullptr);
        num_threads_ = other.num_threads_;
    }
    Spraylist& operator=(Spraylist&& other) noexcept {
        pq_ = std::exchange(other.pq_, nullptr);
        return *this;
    }
    ~Spraylist() {
        ::sl_set_delete(pq_);
    }

    static void write_human_readable(std::ostream& out) {
        out << "Spraylist\n";
    }

    Handle get_handle() {
        return Handle{*pq_, new_thread_data()};
    }
};

}  // namespace wrapper::spraylist
