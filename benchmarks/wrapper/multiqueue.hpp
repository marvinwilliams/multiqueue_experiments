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

    using multiqueue_type = ::multiqueue::MultiQueue<key_type, value_type, KeyOfValue, key_compare, Policy, pq_type>;

    auto parse_config(cxxopts::ParseResult const &result) {
        typename multiqueue_type::config_type config{};
        if constexpr (Policy::has_stickiness) {
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
        : mq_{static_cast<std::size_t>(num_threads * (result.count("factor") > 0 ? result["factor"].as<int>() : 2)),
              initial_capacity, parse_config(result)} {
    }

    static void add_cmd_options(cxxopts::Options &options) {
        options.add_options()("c,factor", "The number of queues per thread", cxxopts::value<int>(), "NUMBER");
        if (Policy::has_stickiness) {
            options.add_options()("k,stickiness", "The stickiness period", cxxopts::value<int>(), "NUMBER");
        }
    }

    void write_json(std::ostream &out) const {
        out << '{';
        out << std::quoted("name") << ':' << std::quoted("MultiQueue");
        out << ',' << std::quoted("configuration") << ':';
        out << '{';
        out << std::quoted("mode") << ':' << std::quoted(Policy::name);
        out << ',' << std::quoted("pop_candidates") << ':' << num_pop_candidates;
#ifdef MQ_USE_STD_PQ
        out << ',' << std::quoted("pq") << ':';
        out << '{';
        out << std::quoted("name") << ':' << std::quoted("buffered");
        out << ',' << std::quoted("type") << ':' << std::quoted("std::priority_queue");
        out << ',' << std::quoted("insertion_buffersize") << ':' << insertion_buffersize;
        out << ',' << std::quoted("deletion_buffersize") << ':' << deletion_buffersize;
        out << '}';
#elif defined MQ_USE_BTREE
        out << ',' << std::quoted("pq") << ':' << std::quoted("tlx::btree");
#else
        out << ',' << std::quoted("pq") << ':';
        out << '{';
        out << std::quoted("name") << ':' << std::quoted("buffered");
        out << ',' << std::quoted("type") << ':';
        out << '{';
        out << std::quoted("name") << ':' << std::quoted("heap");
        out << ',' << std::quoted("arity") << ':' << heap_arity;
        out << '}';
        out << ',' << std::quoted("insertion_buffersize") << ':' << insertion_buffersize;
        out << ',' << std::quoted("deletion_buffersize") << ':' << deletion_buffersize;
        out << '}';
#endif
        out << '}';
        out << ',' << std::quoted("options") << ':';
        out << '{' << std::quoted("num_pqs") << ':' << mq_.num_pqs();
        if constexpr (Policy::has_stickiness) {
            out << ',' << std::quoted("stickiness") << ':' << mq_.config().stickiness;
        }
        out << '}';
        out << '}';
    }

    auto get_handle() {
        return mq_.get_handle();
    }
};

}  // namespace wrapper
