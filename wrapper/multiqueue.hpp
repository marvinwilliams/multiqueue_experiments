#pragma once

#include "multiqueue/buffered_pq.hpp"
#include "multiqueue/multiqueue.hpp"
#include "multiqueue/utils.hpp"

#ifdef MQ_USE_STD_PQ
#include <queue>
#include <vector>
#elif defined MQ_USE_BTREE
#include "tlx_btree.hpp"
#endif

#if defined MQ_MODE_RANDOM || defined MQ_MODE_RANDOM_STRICT
#include "multiqueue/modes/random.hpp"
#elif defined MQ_MODE_STICK_RANDOM
#include "multiqueue/modes/stick_random.hpp"
#elif defined MQ_MODE_STICK_SWAP
#include "multiqueue/modes/stick_swap.hpp"
#elif defined MQ_MODE_STICK_MARK
#include "multiqueue/modes/stick_mark.hpp"
#elif defined MQ_MODE_STICK_PARAMETRIC
#include "multiqueue/modes/stick_parametric.hpp"
#else
#error "No valid mode specified"
#endif

#include <cxxopts.hpp>

#include <iomanip>
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
#elif defined MQ_MODE_STICK_PARAMETRIC
using mode_type = ::multiqueue::mode::StickParametric<num_pop_candidates>;
static constexpr auto mode_name = "stick_parametric";
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
    using key_type = btree_type::key_type;
    using value_type = btree_type::value_type;
    using size_type = btree_type::size_type;
    using key_compare = btree_type::key_compare;
    using value_compare = btree_type::value_compare;

   private:
    btree_type btree_;

   public:
    BTreePQWrapper() = default;
    explicit BTreePQWrapper(key_compare const &comp) : btree_(comp) {
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

template <bool Min, typename Key = unsigned long, typename T = Key>
class MultiQueue {
   public:
    using key_type = Key;
    using mapped_type = T;
    using value_type = std::pair<key_type, mapped_type>;
    using key_compare = std::conditional_t<Min, std::greater<>, std::less<>>;
    using value_compare = ::multiqueue::utils::ValueCompare<value_type, ::multiqueue::utils::PairFirst, key_compare>;

#ifdef MQ_USE_BTREE
    using pq_type = BTreeWrapper<key_type, value_type, KeyOfValue, key_compare>;
#else
    using pq_type = ::multiqueue::BufferedPQ<
#ifdef MQ_USE_STD_PQ
        std::priority_queue<value_type, std::vector<value_type>, value_compare>
#else
        ::multiqueue::Heap<value_type, value_compare, heap_arity>
#endif
        ,
        insertion_buffer_size, deletion_buffer_size>;
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
                cmd.add_options()("k,stickiness", "The stickiness period", cxxopts::value<int>(config.stickiness),
                                  "NUMBER");
            }
        }

        [[nodiscard]] bool validate() const {
            if (factor <= 1) {
                std::cerr << "Error: Queue factor must be at least 2\n";
                return false;
            }
            if constexpr (has_stickiness) {
                if (config.stickiness <= 0) {
                    std::cerr << "Error: Stickiness must be at least 1\n";
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

        void write_json(std::ostream &out) const {
            out << '{';
            out << std::quoted("queue_factor") << ':' << factor << ',';
            out << std::quoted("seed") << ':' << config.seed;
            if constexpr (has_stickiness) {
                out << ',';
                out << std::quoted("stickiness") << ':' << config.stickiness;
            }
            out << '}';
        }
    };

    class Handle : public multiqueue_type::handle_type {
        friend MultiQueue;
        explicit Handle(multiqueue_type &pq) : multiqueue_type::handle_type{pq.get_handle()} {
        }

       public:
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
#ifdef MQ_USE_BTREE
        out << "  PQ: tlx::btree" << '\n';
#else
#ifdef MQ_USE_STD_PQ
        out << "  PQ: std::priority_queue" << '\n';
#else
        out << "  PQ: d-ary heap" << '\n';
        out << "  Heap arity: " << heap_arity << '\n';
#endif
        out << "  Insertion buffer size: " << insertion_buffer_size << '\n';
        out << "  Deletion buffer size: " << deletion_buffer_size << '\n';
#endif
    }

    auto get_handle() {
        return Handle{mq_};
    }
};

}  // namespace wrapper::multiqueue
