#pragma once

#include "multiqueue/multiqueue.hpp"

#include "multiqueue/queue_selection/global_permutation.hpp"
#include "multiqueue/queue_selection/random.hpp"
#include "multiqueue/queue_selection/stick_random.hpp"
#include "multiqueue/queue_selection/swap_assignment.hpp"

#include "cxxopts.hpp"

#include <ostream>
#include <utility>

#ifndef MQ_QUEUE_SELECTION_POLICY
#define MQ_QUEUE_SELECTION_POLICY 1
#endif

namespace detail {
struct Traits : multiqueue::defaults::Traits {
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
#endif
#ifdef MQ_NOSCAN_ON_FAILED_POP
    static constexpr bool scan_on_failed_pop = false;
#else
    static constexpr bool scan_on_failed_pop = true;
#endif
#ifdef MQ_HEAP_ARITY
    static constexpr unsigned int heap_arity = MQ_HEAP_ARITY;
#else
    static constexpr unsigned int heap_arity = 8;
#endif
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
#ifdef MQ_NUM_POP_PQS
    static constexpr unsigned int num_pop_pqs = MQ_NUM_POP_PQS;
#else
    static constexpr unsigned int num_pop_pqs = 2;
#endif
#if MQ_QUEUE_SELECTION_POLICY == 0
    using queue_selection_policy = multiqueue::queue_selection::Random<num_pop_pqs>;
    static constexpr auto queue_selection_policy_name = "random";
#elif MQ_QUEUE_SELECTION_POLICY == 1
    using queue_selection_policy = multiqueue::queue_selection::StickRandom<num_pop_pqs>;
    static constexpr auto queue_selection_policy_name = "stick random";
#elif MQ_QUEUE_SELECTION_POLICY == 2
    using queue_selection_policy = multiqueue::queue_selection::SwapAssignment<num_pop_pqs>;
    static constexpr auto queue_selection_policy_name = "swap assignment";
#elif MQ_QUEUE_SELECTION_POLICY == 3
    using queue_selection_policy = multiqueue::queue_selection::GlobalPermutation<num_pop_pqs>;
    static constexpr auto queue_selection_policy_name = "global permutation";
#else
#error "Invalid MQ_QUEUE_SELECTION_POLICY"
#endif
};
}  // namespace detail

namespace wrapper {

template <typename Key, typename Value, typename Compare, typename Traits = detail::Traits,
          typename Policy = typename Traits::queue_selection_policy_type>
class MultiQueue : public multiqueue::MultiQueue<Key, Value, Compare, Traits> {
    using base_type = multiqueue::MultiQueue<Key, Value, Compare, Traits>;

   public:
    struct config_type : base_type::queue_selection_config_type {
        int factor = 2;
    };

    static void add_options(cxxopts::Options &options, config_type &config) {
        options.add_options()("c,factor", "The number of queues per thread",
                              cxxopts::value<int>(config.factor)->default_value("2"), "NUMBER");
        options.add_options()("k,stickiness", "The stickiness period", cxxopts::value<int>(config.stickiness),
                              "NUMBER");
    }

    MultiQueue(int num_threads, std::size_t initial_capacity, config_type const &config)
        : base_type(static_cast<std::size_t>(config.factor * num_threads), initial_capacity, config) {
    }

    std::ostream &describe(std::ostream &out) {
        out << "MultiQueue [comparisons: " << (Traits::strict_comparison ? "strict" : "non-strict")
            << ", pop tries: " << Traits::num_pop_tries
            << ", scan on failed pop: " << (Traits::scan_on_failed_pop ? "true" : "false")
            << ", heap arity: " << Traits::heap_arity << ", buffersizes (i/d): " << Traits::insertion_buffersize << '/'
            << Traits::deletion_buffersize << ", queue selection: " << Traits::queue_selection_policy_name
            << ", pop pqs: " << Traits::num_pop_pqs << "] (pqs: " << this->num_pqs()
            << ", stickiness: " << this->get_queue_selection_config().stickiness << ')';
        return out;
    }
};

template <typename Key, typename Value, typename Compare, typename Traits, unsigned N>
class MultiQueue<Key, Value, Compare, Traits, multiqueue::queue_selection::Random<N>>
    : public multiqueue::MultiQueue<Key, Value, Compare, Traits> {
    using base_type = multiqueue::MultiQueue<Key, Value, Compare, Traits>;

   public:
    struct config_type : base_type::queue_selection_config_type {
        int factor = 2;
    };

    static void add_options(cxxopts::Options &options, config_type &config) {
        options.add_options()("c,factor", "The number of queues per thread", cxxopts::value<int>(config.factor),
                              "NUMBER");
    }

    MultiQueue(int num_threads, std::size_t initial_capacity, config_type const &config)
        : base_type(static_cast<std::size_t>(config.factor * num_threads), initial_capacity, config) {
    }

    std::ostream &describe(std::ostream &out) {
        out << "MultiQueue [comparisons: " << (Traits::strict_comparison ? "strict" : "non-strict")
            << ", pop tries: " << Traits::num_pop_tries
            << ", scan on failed pop: " << (Traits::scan_on_failed_pop ? "true" : "false")
            << ", heap arity: " << Traits::heap_arity << ", buffersizes (i/d): " << Traits::insertion_buffersize << '/'
            << Traits::deletion_buffersize << ", queue selection: " << Traits::queue_selection_policy_name
            << ", pop pqs: " << Traits::num_pop_pqs << "] (pqs: " << this->num_pqs() << ')';
        return out;
    }
};

}  // namespace wrapper
