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

#if defined MQ_MODE_RANDOM || defined MQ_MODE_RANDOM_STRICT
#include "multiqueue/modes/random.hpp"
#elif defined MQ_MODE_STICK_RANDOM
#include "multiqueue/modes/stick_random.hpp"
#elif defined MQ_MODE_STICK_SWAP
#include "multiqueue/modes/stick_swap.hpp"
#elif defined MQ_MODE_STICK_PARAMETRIC
#include "multiqueue/modes/stick_parametric.hpp"
#else
#error "No valid mode specified"
#endif

#include "cxxopts.hpp"

#include <iomanip>
#include <ostream>
#include <sstream>
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

    struct Policy {
#if defined MQ_MODE_RANDOM
        using mode_type = ::multiqueue::mode::Random<num_pop_candidates, true>;
        static constexpr auto name = "random";
        static constexpr bool has_stickiness = false;
#elif defined MQ_MODE_RANDOM_STRICT
        using mode_type = ::multiqueue::mode::Random<num_pop_candidates, false>;
        static constexpr auto name = "random_strict";
        static constexpr bool has_stickiness = false;
#elif defined MQ_MODE_STICK_RANDOM
        using mode_type = ::multiqueue::mode::StickRandom<num_pop_candidates>;
        static constexpr auto name = "stick_random";
        static constexpr bool has_stickiness = true;
#elif defined MQ_MODE_STICK_SWAP
        using mode_type = ::multiqueue::mode::StickSwap<num_pop_candidates>;
        static constexpr auto name = "stick_swap";
        static constexpr bool has_stickiness = true;
#elif defined MQ_MODE_STICK_PARAMETRIC
        using mode_type = ::multiqueue::mode::StickParametric<num_pop_candidates>;
        static constexpr auto name = "stick_parametric";
        static constexpr bool has_stickiness = true;
#endif
        static constexpr int pop_tries = 1;
        static constexpr bool scan = true;
    };

#ifdef MQ_USE_STD_PQ
    using pq_type = ::multiqueue::BufferedPQ<std::priority_queue<value_type, std::vector<value_type>, value_compare>,
                                             insertion_buffer_size, deletion_buffer_size>;
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
        insertion_buffer_size, deletion_buffer_size>;
#endif

    using multiqueue_type = ::multiqueue::MultiQueue<key_type, value_type, KeyOfValue, key_compare, Policy, pq_type>;

    multiqueue_type mq_;

   public:
    struct Settings {
        int factor = 2;
        typename multiqueue_type::config_type config{};

        static void add_cmd_options(cxxopts::Options &cmd) {
            cmd.add_options()("c,factor", "The number of queues per thread", cxxopts::value<int>(), "NUMBER");
            if (Policy::has_stickiness) {
                cmd.add_options()("k,stickiness", "The stickiness period", cxxopts::value<int>(), "NUMBER");
            }
        }

        static std::optional<Settings> from_cmd(cxxopts::ParseResult const &args) {
            Settings settings{};
            if (args.count("factor") > 0) {
                settings.factor = args["factor"].as<int>();
                if (settings.factor <= 1) {
                    std::cerr << "Error: Queue factor must be at least 2\n";
                    return {};
                }
            }
            if constexpr (Policy::has_stickiness) {
                if (args.count("stickiness") > 0) {
                    settings.config.stickiness = args["stickiness"].as<int>();
                    if (settings.config.stickiness <= 0) {
                        std::cerr << "Error: Stickiness must be at least 1\n";
                        return {};
                    }
                }
            }
            return settings;
        }

        void describe(std::ostream &out) const {
            out << "Queue factor: " << factor << '\n';
            if constexpr (Policy::has_stickiness) {
                out << "Stickiness: " << config.stickiness << '\n';
            }
        }

        void write_json(std::ostream &out) const {
            out << '{';
            out << std::quoted("factor") << ':' << factor;
            if constexpr (Policy::has_stickiness) {
                out << ',' << std::quoted("stickiness") << ':' << config.stickiness;
            }
            out << '}';
        }
    };

    using handle_type = typename multiqueue_type::handle_type;

    explicit MultiQueue(int num_threads, std::size_t initial_capacity, Settings const &settings)
        : mq_{static_cast<std::size_t>(num_threads * settings.factor), initial_capacity, settings.config} {
    }

    static void describe(std::ostream &out) {
        out << "MultiQueue\n";
        out << "  Mode: " << Policy::name << '\n';
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
        return mq_.get_handle();
    }
};

}  // namespace wrapper
