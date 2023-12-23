#pragma once

#include "multiqueue/buffered_pq.hpp"
#include "multiqueue/multiqueue.hpp"
#include "multiqueue/utils.hpp"
#include "util.hpp"

#ifdef MQ_USE_STD_PQ
#include <queue>
#include <vector>
#elif defined MQ_USE_BTREE
#include "tlx_btree.hpp"
#endif

#ifndef MQ_OP_POLICY
#error "MQ_OP_POLICY not defined"
#endif

#if MQ_OP_POLICY == 0
#include "multiqueue/operation_policy/noop.hpp"
#elif MQ_OP_POLICY == 1
#include "multiqueue/stick_policy/random.hpp"
#elif MQ_OP_POLICY == 2
#include "multiqueue/stick_policy/swap.hpp"
#elif MQ_OP_POLICY == 3
#include "multiqueue/stick_policy/parametric.hpp"
#elif MQ_OP_POLICY == 3
#include "multiqueue/stick_policy/parametric.hpp"
#else
#error "Invalid MQ_OP_POLICY"
#endif

#include "cxxopts.hpp"

#include <ostream>
#include <utility>

namespace wrapper {

template <bool Min, typename Key = unsigned long, typename Value = std::pair<Key, Key>,
          typename KeyOfValue = util::KeyOfValue<Key, Value>>
class MultiQueue {
    using key_type = Key;
    using value_type = Value;
    using key_compare = std::conditional_t<Min, std::greater<>, std::less<>>;
    using value_compare = ::multiqueue::utils::ValueCompare<Value, KeyOfValue, key_compare>;

#ifdef MQ_NUM_POP_CANDIDATES
    static constexpr unsigned int num_pop_candidates = MQ_NUM_POP_CANDIDATES;
#else
    static constexpr unsigned int num_pop_candidates = 2;
#endif

#ifdef MQ_NO_POP_STALE
    static constexpr bool pop_stale = false;
#else
    static constexpr bool pop_stale = false;
#endif

#ifdef MQ_INSERTION_BUFFER_SIZE
    static constexpr std::size_t insertion_buffersize = MQ_INSERTION_BUFFER_SIZE;
#else
    static constexpr std::size_t insertion_buffersize = 16;
#endif

#ifdef MQ_DELETION_BUFFER_SIZE
    static constexpr std::size_t deletion_buffersize = MQ_DELETION_BUFFER_SIZE;
#else
    static constexpr std::size_t deletion_buffersize = 16;
#endif

#ifdef MQ_HEAP_ARITY
    static constexpr unsigned int heap_arity = MQ_HEAP_ARITY;
#else
    static constexpr unsigned int heap_arity = 8;
#endif

    struct OperationTraits {
#if MQ_OP_POLICY == 0
        using policy_type = ::multiqueue::operation_policy::Random<num_pop_candidates, pop_stale>;
        static constexpr auto name = "random";
        static constexpr bool has_stickiness = false;
#elif MQ_OP_POLICY == 1
        using policy_type = ::multiqueue::operation_policy::StickRandomIndependent<num_pop_candidates>;
        static constexpr auto name = "stick random independent";
        static constexpr bool has_stickiness = true;
#elif MQ_OP_POLICY == 2
        using policy_type = ::multiqueue::operation_policy::StickRandomCommon<num_pop_candidates>;
        static constexpr auto name = "stick random common";
        static constexpr bool has_stickiness = true;
#elif MQ_OP_POLICY == 3
        using policy_type = ::multiqueue::operation_policy::StickPermutationIndividual<num_pop_candidates>;
        static constexpr auto name = "stick permutation individual";
        static constexpr bool has_stickiness = true;
#elif MQ_OP_POLICY == 4
        using policy_type = ::multiqueue::operation_policy::StickPermutationGlobal<num_pop_candidates>;
        static constexpr auto name = "stick permutation global";
        static constexpr bool has_stickiness = true;
#endif
#ifdef MQ_NOSCAN_IF_EMPTY
        static constexpr bool scan_if_empty = false;
#else
        static constexpr bool scan_if_empty = true;
#endif
    };

#ifdef MQ_USE_STD_PQ
    using pq_type = ::multiqueue::BufferedPQ<std::priority_queue<value_type, std::vector<value_type>, value_compare>,
                                             insertion_buffersize, deletion_buffersize>;
#elif defined MQ_USE_BTREE
    class pq_type {
        using btree_type = tlx::BTree<key_type, value_type, KeyOfValue, key_compare,
                                      tlx::btree_default_traits<key_type, value_type>, true>;

       public:
        using key_type = btree_type::key_type;
        using value_type = btree_type::value_type;
        using size_type = btree_type::size_type;
        using key_compare = btree_type::key_compare;
        using value_compare = btree_type::value_compare;

       private:
        btree_type btree_;

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

#else
    using pq_type = ::multiqueue::BufferedPQ<
        ::multiqueue::Heap<Value, ::multiqueue::utils::ValueCompare<value_type, KeyOfValue, key_compare>, heap_arity>,
        insertion_buffersize, deletion_buffersize>;
#endif

    using multiqueue_type =
        ::multiqueue::MultiQueue<key_type, value_type, KeyOfValue, key_compare, OperationTraits, pq_type>;

    auto parse_config(cxxopts::ParseResult const &result) {
        typename multiqueue_type::operation_config_type config{};
        if constexpr (OperationTraits::has_stickiness) {
            if (result.count("stickiness") > 0) {
                config.stickiness = result["stickiness"].as<int>();
            }
        }
        return config;
    }

    multiqueue_type mq_;

   public:
    using handle_type = typename multiqueue_type::handle_type;

    explicit MultiQueue(int num_threads, std::size_t initial_capacity, cxxopts::ParseResult const &result)
        : mq_{num_threads * (result.count("factor") > 0 ? result["factor"].as<int>() : 2), initial_capacity,
              parse_config(result)} {
    }

    static void add_options(cxxopts::Options &options) {
        options.add_options()("c,factor", "The number of queues per thread", cxxopts::value<int>(), "NUMBER");
        if (OperationTraits::has_stickiness) {
            options.add_options()("k,stickiness", "The stickiness period", cxxopts::value<int>(), "NUMBER");
        }
    }

    std::ostream &describe(std::ostream &out) {
        out << "MultiQueue ("
            << "operations: " << OperationTraits::name << ", pop stale: " << (pop_stale ? "true" : "false")
            << ", pop candidates: " << num_pop_candidates
            << ", scan if empty: " << (OperationTraits::scan_if_empty ? "true" : "false") << ", pq type: "
#ifdef MQ_USE_STD_PQ
            << "buffered std::priority_queue, buffer size (i/d): " << insertion_buffersize << '/'
            << deletion_buffersize;
#elif defined MQ_USE_BTREE
            << "tlx::btree";
#else
            << "buffered d-ary heap"
            << ", arity: " << heap_arity << ", buffer size [ins/del]: " << insertion_buffersize << '/'
            << deletion_buffersize;
#endif
        out << ", pqs: " << mq_.num_pqs();
        if constexpr (OperationTraits::has_stickiness) {
            out << ", stickiness: " << mq_.get_stick_policy_config().stickiness;
        }
        out << ')';
        return out;
    }

    auto get_handle() {
        return mq_.get_handle();
    }
};

}  // namespace wrapper
