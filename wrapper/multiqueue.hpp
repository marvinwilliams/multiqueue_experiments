#pragma once

#include "multiqueue/multiqueue.hpp"

#include "multiqueue/queue_selection/global_permutation.hpp"
#include "multiqueue/queue_selection/random.hpp"
#include "multiqueue/queue_selection/stick_random.hpp"
#include "multiqueue/queue_selection/swap_assignment.hpp"

#include "cxxopts.hpp"
#include "tlx_btree.hpp"

#include <ostream>
#include <queue>
#include <utility>

#ifndef MQ_QUEUE_SELECTION_POLICY
#define MQ_QUEUE_SELECTION_POLICY 1
#endif

namespace detail {
#ifdef MQ_NUM_POP_PQS
static constexpr unsigned int num_pop_pqs = MQ_NUM_POP_PQS;
#else
static constexpr unsigned int num_pop_pqs = 2;
#endif
#if MQ_QUEUE_SELECTION_POLICY == 0
using queue_selection_policy_type = multiqueue::queue_selection::Random<num_pop_pqs>;
static constexpr auto queue_selection_policy_name = "random";
static constexpr bool has_stickiness = false;
#elif MQ_QUEUE_SELECTION_POLICY == 1
using queue_selection_policy_type = multiqueue::queue_selection::StickRandom<num_pop_pqs>;
static constexpr auto queue_selection_policy_name = "stick random";
static constexpr bool has_stickiness = true;
#elif MQ_QUEUE_SELECTION_POLICY == 2
using queue_selection_policy_type = multiqueue::queue_selection::SwapAssignment<num_pop_pqs>;
static constexpr auto queue_selection_policy_name = "swap assignment";
static constexpr bool has_stickiness = true;
#elif MQ_QUEUE_SELECTION_POLICY == 3
using queue_selection_policy_type = multiqueue::queue_selection::GlobalPermutation<num_pop_pqs>;
static constexpr auto queue_selection_policy_name = "global permutation";
static constexpr bool has_stickiness = true;
#else
#error "Invalid MQ_QUEUE_SELECTION_POLICY"
#endif

struct Traits {
    using queue_selection_policy_type = ::detail::queue_selection_policy_type;
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

#ifdef MQ_BUFFERED_PQ_INSERTION_BUFFER_SIZE
static constexpr std::size_t insertion_buffersize = MQ_BUFFERED_PQ_INSERTION_BUFFER_SIZE;
#else
static constexpr std::size_t insertion_buffersize = 64;
#endif
#ifdef MQ_BUFFERED_PQ_DELETION_BUFFER_SIZE
static constexpr std::size_t deletion_buffersize = MQ_BUFFERED_PQ_DELETION_BUFFER_SIZE;
#else
static constexpr std::size_t deletion_buffersize = 64;
#endif

template <typename Pair>
struct PairFirst {
    static constexpr auto const &get(Pair const &p) noexcept {
        return p.first;
    }
};

template <typename Pair, typename Compare>
struct PairCompare {
    Compare comp;

    bool operator()(Pair const &lhs, Pair const &rhs) const noexcept {
        return comp(lhs.first, rhs.first);
    }
};

template <typename Key, typename Value, typename KeyOfValue, typename Compare>
class BTreePQWrapper {
    tlx::BTree<Key, Value, KeyOfValue, Compare, tlx::btree_default_traits<Key, Value>, true> btree_;

   public:
    using value_type = Value;
    using size_type = std::size_t;
    using reference = value_type &;
    using const_reference = value_type const &;
    using key_compare = Compare;
    using value_compare = Compare;

    BTreePQWrapper() = default;
    explicit BTreePQWrapper(Compare const &comp) : btree_(comp) {
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

    void reserve(size_type capacity) {
        // no-op
    }
};

#ifdef MQ_USE_STD_PQ
template <typename Key, typename Value, typename Compare, typename Traits>
using multiqueue_type = multiqueue::MultiQueue<
    Key, Value, Compare, Traits, PairFirst<Value>,
    multiqueue::BufferedPQ<std::priority_queue<Value, std::vector<Value>, PairCompare<Value, Compare>>,
                           insertion_buffersize, deletion_buffersize>>;

#elif defined MQ_USE_BTREE

template <typename Key, typename Value, typename Compare, typename Traits>
using multiqueue_type = multiqueue::MultiQueue<Key, Value, Compare, Traits, PairFirst<Value>,
                                               BTreePQWrapper<Key, Value, PairFirst<Value>, Compare>>;

#else
#ifdef MQ_HEAP_ARITY
static constexpr unsigned int heap_arity = MQ_HEAP_ARITY;
#else
static constexpr unsigned int heap_arity = 8;
#endif

template <typename Key, typename Value, typename Compare, typename Traits>
using multiqueue_type =
    multiqueue::MultiQueue<Key, Value, Compare, Traits, PairFirst<Value>,
                           multiqueue::BufferedPQ<multiqueue::Heap<Value, PairCompare<Value, Compare>, heap_arity>,
                                                  insertion_buffersize, deletion_buffersize>>;
#endif

template <typename Policy>
struct PolicyTraits {
    static constexpr bool has_stickiness = false;
    static constexpr auto name = "unknown";
};

template <unsigned int N>
struct PolicyTraits<multiqueue::queue_selection::Random<N>> {
    static constexpr bool has_stickiness = false;
    static constexpr auto name = "random";
    static constexpr unsigned int num_queues = N;
};

template <unsigned int N>
struct PolicyTraits<multiqueue::queue_selection::StickRandom<N>> {
    static constexpr bool has_stickiness = true;
    static constexpr auto name = "stick random";
    static constexpr unsigned int num_queues = N;
};

template <unsigned int N>
struct PolicyTraits<multiqueue::queue_selection::SwapAssignment<N>> {
    static constexpr bool has_stickiness = true;
    static constexpr auto name = "swap";
    static constexpr unsigned int num_queues = N;
};

template <unsigned int N>
struct PolicyTraits<multiqueue::queue_selection::GlobalPermutation<N>> {
    static constexpr bool has_stickiness = true;
    static constexpr auto name = "permute";
    static constexpr unsigned int num_queues = N;
};

template <typename PQ>
struct PQTraits {
    static std::ostream &describe(std::ostream &os) {
        return os << "unknown";
    }
};

template <typename Value, typename Compare, std::size_t InsertionBuffersize, std::size_t DeletionBuffersize,
          unsigned int Arity>
struct PQTraits<
    multiqueue::BufferedPQ<multiqueue::Heap<Value, Compare, Arity>, InsertionBuffersize, DeletionBuffersize>> {
    static std::ostream &describe(std::ostream &os) {
        os << "d-ary heap (d: " << Arity << "), buffersizes (i/d): " << InsertionBuffersize << '/'
           << DeletionBuffersize;
        return os;
    }
};

template <typename Value, typename Compare, std::size_t InsertionBuffersize, std::size_t DeletionBuffersize>
struct PQTraits<multiqueue::BufferedPQ<std::priority_queue<Value, std::vector<Value>, Compare>, InsertionBuffersize,
                                       DeletionBuffersize>> {
    static std::ostream &describe(std::ostream &os) {
        os << "std::priority_queue, buffersizes (i/d): " << InsertionBuffersize << '/' << DeletionBuffersize;
        return os;
    }
};

template <typename Key, typename Value, typename KeyOfValue, typename Compare>
struct PQTraits<detail::BTreePQWrapper<Key, Value, KeyOfValue, Compare>> {
    static std::ostream &describe(std::ostream &os) {
        os << "tlx::btree";
        return os;
    }
};

}  // namespace detail

namespace wrapper {
template <typename Key, typename Value, typename Compare, typename Traits = detail::Traits>
class MultiQueue : public detail::multiqueue_type<Key, Value, Compare, Traits> {
    using base_type = detail::multiqueue_type<Key, Value, Compare, Traits>;
    using policy_traits = detail::PolicyTraits<typename base_type::traits_type::queue_selection_policy_type>;
    using pq_traits = detail::PQTraits<typename base_type::priority_queue_type>;

   public:
    struct config_type : base_type::queue_selection_config_type {
        int factor = 2;
    };

    static void add_options(cxxopts::Options &options, config_type &config) {
        options.add_options()("c,factor", "The number of queues per thread",
                              cxxopts::value<int>(config.factor)->default_value("2"), "NUMBER");
        if constexpr (policy_traits::has_stickiness) {
            options.add_options()("k,stickiness", "The stickiness period", cxxopts::value<int>(config.stickiness),
                                  "NUMBER");
        }
    }

    MultiQueue(int num_threads, std::size_t initial_capacity, config_type const &config)
        : base_type(static_cast<std::size_t>(config.factor * num_threads), initial_capacity, config) {
    }

    std::ostream &describe(std::ostream &out) {
        out << "MultiQueue (comparisons: " << (Traits::strict_comparison ? "strict" : "non-strict")
            << ", pop tries: " << Traits::num_pop_tries
            << ", scan on failed pop: " << (Traits::scan_on_failed_pop ? "true" : "false")
            << ", queue selection: " << policy_traits::name << ", pop pqs: " << policy_traits::num_queues
            << ", inner pq: ";
        pq_traits::describe(out);
        out << ", pqs: " << this->num_pqs();
        if constexpr (policy_traits::has_stickiness) {
            out << ", stickiness: " << this->get_queue_selection_config().stickiness;
        }
        out << ')';
        return out;
    }
};

}  // namespace wrapper
