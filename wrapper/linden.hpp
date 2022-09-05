#pragma once
#ifndef WRAPPER_LINDEN_HPP_INCLUDED
#define WRAPPER_LINDEN_HPP_INCLUDED

// Adapted from klsm

#include <limits>
#include <memory>
#include <string>
#include <utility>

struct linden_pq_t;
struct linden_pq_deleter {
    void operator()(linden_pq_t*);
};

namespace wrapper {

class Linden {
   public:
    using key_type = unsigned long;
    using mapped_type = unsigned long;

    using value_type = std::pair<key_type, mapped_type>;

    struct Handle {
        linden_pq_t* pq_;

        void push(value_type const& value);
        bool try_pop(value_type& retval);
    };

    using handle_type = Handle;

    // The queue itself only supports
    // keys >= 1, so one is added on each insert
    static constexpr key_type min_valid_key =
        std::numeric_limits<key_type>::min();
    static constexpr key_type max_valid_key =
        std::numeric_limits<key_type>::max() - 3;

   private:
    static constexpr key_type sentinel_ = std::numeric_limits<key_type>::max();

    struct wrapper_type;

    alignas(64) std::unique_ptr<linden_pq_t, linden_pq_deleter> pq_;

   public:
    Linden(unsigned int num_threads);

    Handle get_handle();

    void push(value_type const& value);
    bool try_pop(value_type& retval);

    static std::string description() { return "linden"; }
};

}  // namespace wrapper

#endif
