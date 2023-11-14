#pragma once

#include "multiqueue/multiqueue.hpp"
#include "util.hpp"

namespace wrapper::multiqueue::detail::stick_policy {
#ifdef MQ_NUM_POP_PQS
static constexpr unsigned int num_pop_pqs = MQ_NUM_POP_PQS;
#else
static constexpr unsigned int num_pop_pqs = 2;
#endif
}  // namespace wrapper::multiqueue::detail::stick_policy

#ifndef MQ_STICK_POLICY
#error "MQ_STICK_POLICY not defined"
#endif

#if MQ_STICK_POLICY == 0

#include "multiqueue/stick_policy/noop.hpp"

namespace wrapper::multiqueue::detail::stick_policy {

using type = ::multiqueue::stick_policy::Noop<num_pop_pqs>;
static constexpr auto name = "noop";
static constexpr bool has_stickiness = false;

}  // namespace wrapper::multiqueue::detail::stick_policy

#elif MQ_STICK_POLICY == 1

#include "multiqueue/stick_policy/random.hpp"

namespace wrapper::multiqueue::detail::stick_policy {

using type = ::multiqueue::stick_policy::Random<num_pop_pqs>;
static constexpr auto name = "random";
static constexpr bool has_stickiness = true;
}  // namespace wrapper::multiqueue::detail::stick_policy

#elif MQ_STICK_POLICY == 2

#include "multiqueue/stick_policy/swap.hpp"

namespace wrapper::multiqueue::detail::stick_policy {

using type = ::multiqueue::stick_policy::Swap<num_pop_pqs>;
static constexpr auto name = "swap";
static constexpr bool has_stickiness = true;

}  // namespace wrapper::multiqueue::detail::stick_policy

#elif MQ_STICK_POLICY == 3

#include "multiqueue/stick_policy/parametric.hpp"

namespace wrapper::multiqueue::detail::stick_policy {

using type = ::multiqueue::stick_policy::Parametric<num_pop_pqs>;
static constexpr auto name = "parametric";
static constexpr bool has_stickiness = true;

}  // namespace wrapper::multiqueue::detail::stick_policy

#else
#error "Invalid MQ_STICK_POLICY"
#endif

#include "cxxopts.hpp"
#ifdef MQ_USE_BTREE
#include "tlx_btree.hpp"
#endif

#include <ostream>
#include <queue>
#include <utility>

namespace wrapper::multiqueue {

namespace detail {

#ifdef MQ_BUFFERED_PQ_INSERTION_BUFFER_SIZE
static constexpr std::size_t insertion_buffersize = MQ_BUFFERED_PQ_INSERTION_BUFFER_SIZE;
#else
static constexpr std::size_t insertion_buffersize = 16;
#endif
#ifdef MQ_BUFFERED_PQ_DELETION_BUFFER_SIZE
static constexpr std::size_t deletion_buffersize = MQ_BUFFERED_PQ_DELETION_BUFFER_SIZE;
#else
static constexpr std::size_t deletion_buffersize = 16;
#endif
#ifdef MQ_HEAP_ARITY
static constexpr unsigned int heap_arity = MQ_HEAP_ARITY;
#else
static constexpr unsigned int heap_arity = 8;
#endif

struct Traits {
    using stick_policy_type = stick_policy::type;
#ifdef MQ_COMPARE_STRICT
    static constexpr bool strict_comparison = true;
#else
    static constexpr bool strict_comparison = false;
#endif
#ifdef MQ_COUNT_STATS
    static constexpr bool count_stats = true;
#else
    static constexpr bool count_stats = false;
#endif
#ifdef MQ_NUM_POP_TRIES
    static constexpr unsigned int num_pop_tries = MQ_NUM_POP_TRIES;
#else
    static constexpr unsigned int num_pop_tries = 1;
#endif
#ifdef MQ_NOSCAN_ON_FAILED_POP
    static constexpr bool scan_on_failed_pop = false;
#else
    static constexpr bool scan_on_failed_pop = true;
#endif
};

template <bool Min, typename Key, typename Value, typename KeyOfValue>
struct MultiQueueBuilder {
    using key_type = Key;
    using value_type = Value;
    using key_compare = std::conditional_t<Min, std::greater<key_type>, std::less<key_type>>;
    using key_of_value = KeyOfValue;

    struct value_compare {
        key_compare comp;

        bool operator()(value_type const &lhs, value_type const &rhs) const noexcept {
            return comp(key_of_value::get(lhs), key_of_value::get(rhs));
        }
    };

#ifdef MQ_USE_STD_PQ
    using pq_type = std::multiqueue::BufferedPQ<std::priority_queue<value_type, std::vector<value_type>, value_compare>,
                                                insertion_buffersize, deletion_buffersize>;
#elif defined MQ_USE_BTREE
    class BTreePQWrapper {
        using btree_type = tlx::BTree<key_type, value_type, key_of_value, key_compare,
                                      tlx::btree_default_traits<key_type, value_type>, true>;

       public:
        using key_type = btree_type::key_type;
        using value_type = btree_type::value_type;
        using size_type = btree_type::size_type;
        using key_compare = btree_type::key_compare;
        using value_compare = btree_type::value_compare;

       private:
        tlx::BTree<key_type, value_type, KeyOfValue, compare_type, tlx::btree_default_traits<key_type, value_type>,
                   true>
            btree_;

       public:
        BTreePQWrapper() = default;
        explicit BTreePQWrapper(compare_type const &comp) : btree_(comp) {
        }

        void push(value_type const &value) {
            btree_.insert(value);
        }

        void pop() {
            // end() because comparator is reversed
            btree_.erase(--btree_.end());
        }

        value_type const &top() const {
            // end() because comparator is reversed
            return *(--btree_.end());
        }

        [[nodiscard]] bool empty() const {
            return btree_.empty();
        }

        [[nodiscard]] size_type size() const {
            return btree_.size();
        }

        void clear() {
            btree_.clear();
        }

        void reserve(size_type /*capacity*/) {
            // no-op
        }
    };

    using pq_type = BTreePQWrapper;
#else
    using pq_type = ::multiqueue::BufferedPQ<::multiqueue::Heap<value_type, value_compare, heap_arity>,
                                             insertion_buffersize, deletion_buffersize>;
#endif
    using multiqueue_type = ::multiqueue::MultiQueue<key_type, key_compare, value_type, key_of_value, Traits, pq_type>;
};

}  // namespace detail

template <bool Min = true, typename Key = unsigned long, typename Value = std::pair<unsigned long, unsigned long>,
          typename KeyOfValue = util::KeyOfValue<Key, Value>>
using PQWrapper = typename detail::MultiQueueBuilder<Min, Key, Value, KeyOfValue>::multiqueue_type;

inline void add_options(cxxopts::Options &options) {
    options.add_options()("c,factor", "The number of queues per thread", cxxopts::value<int>(), "NUMBER");
    if (detail::stick_policy::has_stickiness) {
        options.add_options()("k,stickiness", "The stickiness period", cxxopts::value<int>(), "NUMBER");
    }
}

template <bool Min = true, typename Key = unsigned long, typename Value = std::pair<unsigned long, unsigned long>,
          typename KeyOfValue = util::KeyOfValue<Key, Value>>
PQWrapper<Min, Key, Value, KeyOfValue> create(int num_threads, std::size_t initial_capacity,
                                            cxxopts::ParseResult const &result) {
    typename PQWrapper<Min, Key, Value, KeyOfValue>::stick_policy_config_type config{};
    int factor = result.count("factor") > 0 ? result["factor"].as<int>() : 2;
    if constexpr (detail::stick_policy::has_stickiness) {
        if (result.count("stickiness") > 0) {
            config.stickiness = result["stickiness"].as<int>();
        }
    }
    return PQWrapper<Min, Key, Value, KeyOfValue>{num_threads * factor, initial_capacity, config};
}

template <typename MultiQueue>
std::ostream &describe(MultiQueue const &mq, std::ostream &out) {
    out << "MultiQueue (comparisons: " << (detail::Traits::strict_comparison ? "strict" : "non-strict")
        << ", pop tries: " << detail::Traits::num_pop_tries
        << ", scan on failed pop: " << (detail::Traits::scan_on_failed_pop ? "true" : "false")
        << ", stick policy: " << detail::stick_policy::name
        << ", pop pqs: " << detail::stick_policy::num_pop_pqs;
    out << ", pq type: "
#ifdef MQ_USE_STD_PQ
        << "buffered std::priority_queue, buffer size (i/d): " << insertion_buffersize << '/' << deletion_buffersize;
#elif defined MQ_USE_BTREE
        << "tlx::btree";
#else
        << "buffered d-ary heap (d: " << detail::heap_arity << ')'
        << ", buffer size (i/d): " << detail::insertion_buffersize << '/' << detail::deletion_buffersize << ")";
#endif
    out << ", pqs: " << mq.num_pqs();
    if constexpr (detail::stick_policy::has_stickiness) {
        out << ", stickiness: " << mq.get_stick_policy_config().stickiness;
    }
    out << ')';
    return out;
}

}  // namespace wrapper::multiqueue
