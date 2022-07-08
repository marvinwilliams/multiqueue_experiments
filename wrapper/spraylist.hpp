#pragma once
#ifndef WRAPPER_SPRAYLIST_HPP_INCLUDED
#define WRAPPER_SPRAYLIST_HPP_INCLUDED

// Adapted from klsm

#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <utility>

typedef struct sl_intset sl_intset_t;
struct sl_intset_deleter {
    void operator()(sl_intset_t*);
};

typedef struct thread_data thread_data_t;
struct thread_data_deleter {
    void operator()(thread_data_t*);
};

namespace wrapper {

class Spraylist {
   public:
    using key_type = unsigned long;
    using mapped_type = unsigned long;
    using value_type = std::pair<key_type, mapped_type>;

    struct Handle {
        std::unique_ptr<thread_data_t, thread_data_deleter> data_;
        sl_intset_t* pq_;

        void push(value_type const& value);
        bool try_pop(value_type& retval);
    };

    static constexpr key_type min_valid_key = std::numeric_limits<key_type>::min();
    // Use INT_MAX_32
    static constexpr key_type max_valid_key =
        std::numeric_limits<std::uint32_t>::max() - 1ul;
    static_assert(std::numeric_limits<unsigned long>::max() == -1ul, "");

   private:
    static constexpr key_type sentinel_ = std::numeric_limits<key_type>::max();

    alignas(64) std::unique_ptr<sl_intset_t, sl_intset_deleter> pq_;
    unsigned int num_threads_;
    std::unique_ptr<thread_data_t, thread_data_deleter> thread_data_;

   public:
    explicit Spraylist(unsigned int num_threads);

    Handle get_handle();

    void push(value_type const& value);
    bool try_pop(value_type& retval);

    static std::string description() { return "spraylist"; }
};

}  // namespace wrapper

#endif
