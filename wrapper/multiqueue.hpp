#pragma once

#include "multiqueue/buffered_pq.hpp"
#include "multiqueue/multiqueue.hpp"
#include "multiqueue/utils.hpp"

#include "util/base.hpp"

#ifdef MQ_USE_STD_PQ
#include <queue>
#include <vector>
#elif defined MQ_USE_BTREE
#include "tlx_btree.hpp"
#elif defined MQ_USE_MERGE_HEAP
#include "util/merge_heap.hpp"
#endif

#if defined MQ_MODE_RANDOM || defined MQ_MODE_RANDOM_STRICT
#include "multiqueue/modes/random.hpp"
#elif defined MQ_MODE_STICK_RANDOM
#include "multiqueue/modes/stick_random.hpp"
#elif defined MQ_MODE_STICK_SWAP
#include "multiqueue/modes/stick_swap.hpp"
#elif defined MQ_MODE_STICK_MARK
#include "multiqueue/modes/stick_mark.hpp"
#elif defined MQ_MODE_STICK_REPLACE
#include "multiqueue/modes/stick_replace.hpp"
#else
#error "No valid mode specified"
#endif

#include <cxxopts.hpp>

#include <type_traits>

#include <ostream>
#include <utility>

namespace wrapper::multiqueue {

#ifdef MQ_NUM_POP_PQS
static constexpr unsigned int num_pop_candidates = MQ_NUM_POP_PQS;
#else
static constexpr unsigned int num_pop_candidates = 2;
#endif

#ifdef MQ_INSERTION_BUFFER_SIZE
static constexpr std::size_t insertion_buffer_size = MQ_INSERTION_BUFFER_SIZE;
#else
static constexpr std::size_t insertion_buffer_size = 16;
#endif

#ifdef MQ_DELETION_BUFFER_SIZE
static constexpr std::size_t deletion_buffer_size = MQ_DELETION_BUFFER_SIZE;
#else
static constexpr std::size_t deletion_buffer_size = 16;
#endif

#ifdef MQ_HEAP_ARITY
static constexpr unsigned int heap_arity = MQ_HEAP_ARITY;
#else
static constexpr unsigned int heap_arity = 8;
#endif

#if defined MQ_MODE_RANDOM
using mode_type = ::multiqueue::mode::Random<num_pop_candidates, true>;
static constexpr auto mode_name = "random";
static constexpr bool has_stickiness = false;
#elif defined MQ_MODE_RANDOM_STRICT
using mode_type = ::multiqueue::mode::Random<num_pop_candidates, false>;
static constexpr auto mode_name = "random_strict";
static constexpr bool has_stickiness = false;
#elif defined MQ_MODE_STICK_RANDOM
using mode_type = ::multiqueue::mode::StickRandom<num_pop_candidates>;
static constexpr auto mode_name = "stick_random";
static constexpr bool has_stickiness = true;
#elif defined MQ_MODE_STICK_SWAP
using mode_type = ::multiqueue::mode::StickSwap<num_pop_candidates>;
static constexpr auto mode_name = "stick_swap";
static constexpr bool has_stickiness = true;
#elif defined MQ_MODE_STICK_MARK
using mode_type = ::multiqueue::mode::StickMark<num_pop_candidates>;
static constexpr auto mode_name = "stick_mark";
static constexpr bool has_stickiness = true;
#elif defined MQ_MODE_STICK_REPLACE
using mode_type = ::multiqueue::mode::StickReplace<num_pop_candidates>;
static constexpr auto mode_name = "stick_replace";
static constexpr bool has_stickiness = true;
#endif

struct Policy {
    using mode_type = ::wrapper::multiqueue::mode_type;
    static constexpr int pop_tries = 1;
    static constexpr bool scan = true;
};

#ifdef MQ_USE_BTREE
template <typename Key, typename Value, typename KeyOfValue, typename Compare>
class BTreeWrapper {
    using btree_type = tlx::BTree<Key, Value, KeyOfValue, Compare, tlx::btree_default_traits<Key, Value>, true>;

   public:
    using key_type = typename btree_type::key_type;
    using value_type = typename btree_type::value_type;
    using size_type = typename btree_type::size_type;
    using key_compare = typename btree_type::key_compare;
    using value_compare = typename btree_type::value_compare;

   private:
    btree_type btree_;

   public:
    BTreeWrapper() = default;
    explicit BTreeWrapper(key_compare const &comp) : btree_(comp) {
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
#endif

#ifdef MQ_USE_MERGE_HEAP
#ifdef MQ_MERGE_HEAP_NODE_SIZE
static constexpr std::size_t merge_heap_node_size = MQ_MERGE_HEAP_NODE_SIZE;
#else
static constexpr std::size_t merge_heap_node_size = 64;
#endif

// merge_heap pops the key ordered first by its comparator, the library pops the largest key by `Compare`
template <typename Key, typename Value, typename KeyOfValue, typename Compare>
class MergeHeapWrapper {
    struct ExtractKey {
        Key const &operator()(Value const &value) const noexcept {
            return KeyOfValue::get(value);
        }
    };

    struct OrderedFirst {
        [[no_unique_address]] Compare comp;

        bool operator()(Key const &lhs, Key const &rhs) const noexcept {
            return comp(rhs, lhs);
        }
    };

    using heap_type = ::multiqueue::merge_heap<Value, Key, ExtractKey, OrderedFirst, merge_heap_node_size>;

   public:
    using key_type = Key;
    using value_type = Value;
    using size_type = typename heap_type::size_type;
    using key_compare = Compare;
    using value_compare = base::ValueCompare<Value, KeyOfValue, Compare>;
    using reference = value_type &;
    using const_reference = value_type const &;

   private:
    heap_type heap_;

   public:
    MergeHeapWrapper() = default;
    explicit MergeHeapWrapper(key_compare const &comp) : heap_(OrderedFirst{comp}) {
    }

    void push(value_type const &value) {
        heap_.push(value);
    }

    void pop() {
        heap_.pop();
    }

    const_reference top() const {
        return heap_.top();
    }

    [[nodiscard]] bool empty() const noexcept {
        return heap_.empty();
    }

    [[nodiscard]] size_type size() const noexcept {
        return heap_.size();
    }

    void clear() noexcept {
        heap_.clear();
    }

    void reserve(size_type capacity) {
        heap_.reserve(capacity);
    }
};
#endif

static_assert((insertion_buffer_size == 0) == (deletion_buffer_size == 0),
              "Either both or none of the buffers must be disabled");

// The library reserves capacity in every PQ, which the plain heaps only provide through BufferedPQ
template <typename PQ>
struct UnbufferedPQ : PQ {
    using PQ::PQ;

    void reserve(typename PQ::size_type capacity) {
        this->c.reserve(capacity);
    }
};

inline void write_pq_description(std::ostream &out) {
#if defined MQ_USE_BTREE
    out << "  PQ: tlx::btree" << '\n';
#elif defined MQ_USE_MERGE_HEAP
    out << "  PQ: merge heap" << '\n';
    out << "  Node size: " << merge_heap_node_size << '\n';
#else
#ifdef MQ_USE_STD_PQ
    out << "  PQ: std::priority_queue" << '\n';
#else
    out << "  PQ: d-ary heap" << '\n';
    out << "  Heap arity: " << heap_arity << '\n';
#endif
    if constexpr (insertion_buffer_size == 0) {
        out << "  Buffers: none" << '\n';
    } else {
        out << "  Insertion buffer size: " << insertion_buffer_size << '\n';
        out << "  Deletion buffer size: " << deletion_buffer_size << '\n';
    }
#endif
}

template <typename JsonObject>
void write_pq_json(JsonObject &obj) {
    obj.entry("pop_candidates", num_pop_candidates);
#if defined MQ_USE_BTREE
    obj.entry("base_pq", "btree");
#elif defined MQ_USE_MERGE_HEAP
    obj.entry("base_pq", "merge_heap");
    obj.entry("node_size", merge_heap_node_size);
#else
#ifdef MQ_USE_STD_PQ
    obj.entry("base_pq", "std");
#else
    obj.entry("base_pq", "heap");
    obj.entry("heap_arity", heap_arity);
#endif
    obj.entry("insertion_buffer_size", insertion_buffer_size);
    obj.entry("deletion_buffer_size", deletion_buffer_size);
#endif
}

template <bool Min, typename Key = unsigned long, typename T = Key>
class MultiQueue {
   public:
    using key_type = Key;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;
    using key_compare = std::conditional_t<Min, std::greater<>, std::less<>>;
    using value_compare = base::ValueCompare<value_type, base::PairFirst, key_compare>;

#if defined MQ_USE_BTREE
    using pq_type = BTreeWrapper<key_type, value_type, base::PairFirst, key_compare>;
#elif defined MQ_USE_MERGE_HEAP
    using pq_type = MergeHeapWrapper<key_type, value_type, base::PairFirst, key_compare>;
#else
    using base_pq_type =
#ifdef MQ_USE_STD_PQ
        std::priority_queue<value_type, std::vector<value_type>, value_compare>;
#else
        ::multiqueue::Heap<value_type, value_compare, heap_arity>;
#endif
    using pq_type = std::conditional_t<insertion_buffer_size == 0, UnbufferedPQ<base_pq_type>,
                                       ::multiqueue::BufferedPQ<base_pq_type, insertion_buffer_size, deletion_buffer_size>>;
#endif

    using multiqueue_type = ::multiqueue::KeyValueMultiQueue<key_type, mapped_type, key_compare, Policy, pq_type>;

   private:
    multiqueue_type mq_;

    struct Settings {
        int factor = 2;
        typename multiqueue_type::config_type config{};

        void register_cmd_options(cxxopts::Options &cmd) {
            cmd.add_options()("c,queue-factor", "The number of queues per thread", cxxopts::value<int>(factor),
                              "NUMBER");
            cmd.add_options()("mq-seed", "Seed for the multiqueue", cxxopts::value<int>(config.seed), "NUMBER");
            if constexpr (has_stickiness) {
                cmd.add_options()("k,stickiness",
                                  "Number of operations the selected queues are used for (1 reselects every "
                                  "operation)",
                                  cxxopts::value<int>(config.stickiness), "NUMBER");
            }
        }

        [[nodiscard]] bool validate() const {
            if (factor <= 1) {
                std::cerr << "Error: Queue factor must be at least 2\n";
                return false;
            }
            if constexpr (has_stickiness) {
                if (config.stickiness <= 0) {
                    std::cerr << "Error: Stickiness must be positive\n";
                    return false;
                }
            }
            return true;
        }

        void write_human_readable(std::ostream &out) const {
            out << "Queue factor: " << factor << '\n';
            out << "MQ seed: " << config.seed << '\n';
            if constexpr (has_stickiness) {
                out << "Stickiness: " << config.stickiness << '\n';
            }
        }

        template <typename JsonObject>
        void write_json(JsonObject &obj) const {
            obj.entry("mode", mode_name);
            obj.entry("queue_factor", factor);
            obj.entry("seed", config.seed);
            if constexpr (has_stickiness) {
                obj.entry("stickiness", config.stickiness);
            }
            write_pq_json(obj);
        }
    };

    class Handle : public multiqueue_type::handle_type {
        friend MultiQueue;
        explicit Handle(multiqueue_type &pq) : multiqueue_type::handle_type{pq.get_handle()} {
        }

       public:
        using value_type = typename multiqueue_type::value_type;

        bool push(typename multiqueue_type::value_type const &value) {
            multiqueue_type::handle_type::push(value);
            return true;
        }

        std::optional<typename multiqueue_type::value_type> try_pop() {
            return multiqueue_type::handle_type::try_pop();
        }
    };

   public:
    using handle_type = Handle;
    using settings_type = Settings;

    explicit MultiQueue(int num_threads, std::size_t initial_capacity, Settings const &settings)
        : mq_{static_cast<std::size_t>(num_threads * settings.factor), initial_capacity, settings.config} {
    }

    static void write_human_readable(std::ostream &out) {
        out << "MultiQueue\n";
        out << "  Mode: " << mode_name << '\n';
        out << "  Pop candidates: " << num_pop_candidates << '\n';
        write_pq_description(out);
    }

    auto get_handle() {
        return Handle{mq_};
    }
};

}  // namespace wrapper::multiqueue
