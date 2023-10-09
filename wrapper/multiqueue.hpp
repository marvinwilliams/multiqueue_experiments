#pragma once

#include "multiqueue/multiqueue.hpp"

namespace wrapper::detail::queue_selection {
#ifdef MQ_NUM_POP_PQS
static constexpr unsigned int num_pop_pqs = MQ_NUM_POP_PQS;
#else
static constexpr unsigned int num_pop_pqs = 2;
#endif
}  // namespace wrapper::detail::queue_selection

#ifndef MQ_QUEUE_SELECTION_POLICY
#define MQ_QUEUE_SELECTION_POLICY 1
#endif

#if MQ_QUEUE_SELECTION_POLICY == 0

#include "multiqueue/queue_selection/random.hpp"

namespace wrapper::detail::queue_selection {

using type = multiqueue::queue_selection::Random<num_pop_pqs>;
static constexpr auto name = "random";
static constexpr bool has_stickiness = false;

}  // namespace wrapper::detail::queue_selection

#elif MQ_QUEUE_SELECTION_POLICY == 1

#include "multiqueue/queue_selection/stick_random.hpp"

namespace wrapper::detail::queue_selection {

using type = multiqueue::queue_selection::StickRandom<num_pop_pqs>;
static constexpr auto name = "stick random";
static constexpr bool has_stickiness = true;
}  // namespace wrapper::detail::queue_selection

#elif MQ_QUEUE_SELECTION_POLICY == 2

#include "multiqueue/queue_selection/swap_assignment.hpp"

namespace wrapper::detail::queue_selection {

using type = multiqueue::queue_selection::SwapAssignment<num_pop_pqs>;
static constexpr auto name = "swap assignment";
static constexpr bool has_stickiness = true;

}  // namespace wrapper::detail::queue_selection

#elif MQ_QUEUE_SELECTION_POLICY == 3

#include "multiqueue/queue_selection/global_permutation.hpp"

namespace wrapper::detail::queue_selection {

using type = multiqueue::queue_selection::GlobalPermutation<num_pop_pqs>;
static constexpr auto name = "global permutation";
static constexpr bool has_stickiness = true;

}  // namespace wrapper::detail::queue_selection

#else
#error "Invalid MQ_QUEUE_SELECTION_POLICY"
#endif

#include "wrapper/priority.hpp"

#include "cxxopts.hpp"
#ifdef MQ_USE_BTREE
#include "tlx_btree.hpp"
#endif

#include <ostream>
#include <queue>
#include <utility>

namespace wrapper::detail {

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
    using queue_selection_policy_type = queue_selection::type;
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

template <typename Key, typename T, Priority P>
struct MultiQueueBuilder {
    using key_type = Key;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;
    using compare_type = std::conditional_t<P == Priority::Min, std::greater<key_type>, std::less<key_type>>;

    struct ExtractKey {
        static constexpr auto const &get(value_type const &v) noexcept {
            return v.first;
        }
    };

    struct ValueCompare {
        compare_type comp;

        bool operator()(value_type const &lhs, value_type const &rhs) const noexcept {
            return comp(lhs.first, rhs.first);
        }
    };

#ifdef MQ_USE_STD_PQ
    using pq_type = std::multiqueue::BufferedPQ<std::priority_queue<value_type, std::vector<value_type>, ValueCompare>,
                                                insertion_buffersize, deletion_buffersize>;

    static std::ostream &describe_pq(std::ostream &os) {
        return os << "buffered std::priority_queue, buffer size (i/d): " << insertion_buffersize << '/'
                  << deletion_buffersize;
    }

#elif defined MQ_USE_BTREE
    class BTreePQWrapper {
       public:
        using key_type = Key;
        using mapped_type = T;
        using value_type = std::pair<key_type, mapped_type>;
        using size_type = std::size_t;
        using reference = value_type &;
        using const_reference = value_type const &;
        using key_compare = compare_type;
        using value_compare = ValueCompare;

       private:
        tlx::BTree<key_type, value_type, ExtractKey, compare_type, tlx::btree_default_traits<key_type, value_type>,
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

    static std::ostream &describe_pq(std::ostream &os) {
        return os << "tlx::btree";
    }

#else
    using pq_type = multiqueue::BufferedPQ<multiqueue::Heap<value_type, ValueCompare, heap_arity>, insertion_buffersize,
                                           deletion_buffersize>;

    static std::ostream &describe_pq(std::ostream &os) {
        return os << "buffered d-ary heap (d: " << heap_arity << ')' << ", buffer size (i/d): " << insertion_buffersize
                  << '/' << deletion_buffersize << ")";
    }
#endif
    using multiqueue_type = multiqueue::MultiQueue<key_type, value_type, compare_type, Traits, ExtractKey, pq_type>;
};

}  // namespace wrapper::detail

namespace wrapper {

template <typename Key, typename Value, Priority P>
class MultiQueue : public detail::MultiQueueBuilder<Key, Value, P>::multiqueue_type {
    using base_type = typename detail::MultiQueueBuilder<Key, Value, P>::multiqueue_type;

   public:
    struct config_type : base_type::queue_selection_config_type {
        int factor = 2;
    };

    static void add_options(cxxopts::Options &options, config_type &config) {
        options.add_options()("c,factor", "The number of queues per thread",
                              cxxopts::value<int>(config.factor)->default_value("2"), "NUMBER");
        if constexpr (detail::queue_selection::has_stickiness) {
            options.add_options()("k,stickiness", "The stickiness period", cxxopts::value<int>(config.stickiness),
                                  "NUMBER");
        }
    }

    MultiQueue(int num_threads, std::size_t initial_capacity, config_type const &config)
        : base_type(static_cast<std::size_t>(config.factor * num_threads), initial_capacity, config) {
    }

    std::ostream &describe(std::ostream &out) {
        out << "MultiQueue (comparisons: " << (detail::Traits::strict_comparison ? "strict" : "non-strict")
            << ", pop tries: " << detail::Traits::num_pop_tries
            << ", scan on failed pop: " << (detail::Traits::scan_on_failed_pop ? "true" : "false")
            << ", queue selection: " << detail::queue_selection::name
            << ", pop pqs: " << detail::queue_selection::num_pop_pqs;
        if constexpr (detail::queue_selection::has_stickiness) {
            out << ", stickiness: " << this->get_queue_selection_config().stickiness;
        }
        out << ", pq type: ";
        detail::MultiQueueBuilder<Key, Value, P>::describe_pq(out);
        out << ", pqs: " << this->num_pqs();
        out << ')';
        return out;
    }
};

}  // namespace wrapper
