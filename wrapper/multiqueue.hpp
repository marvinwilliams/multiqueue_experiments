#pragma once

#include "multiqueue/multiqueue.hpp"

#include "cxxopts.hpp"

#include <ostream>
#include <utility>

namespace multiqueue::queue_selection {
template <std::size_t>
class Random;
template <std::size_t>
class StickRandom;
template <std::size_t>
class SwapAssignment;
template <std::size_t>
class GlobalPermutation;
}  // namespace multiqueue::queue_selection

namespace detail {
template <typename QueueSelectionPolicy>
struct QueueSelectionPolicyTraits {
    static constexpr auto name = "unknown";
};
template <std::size_t N>
struct QueueSelectionPolicyTraits<multiqueue::queue_selection::Random<N>> {
    static constexpr auto name = "random";
};
template <std::size_t N>
struct QueueSelectionPolicyTraits<multiqueue::queue_selection::StickRandom<N>> {
    static constexpr auto name = "stick random";
};
template <std::size_t N>
struct QueueSelectionPolicyTraits<multiqueue::queue_selection::SwapAssignment<N>> {
    static constexpr auto name = "swap assignment";
};
template <std::size_t N>
struct QueueSelectionPolicyTraits<multiqueue::queue_selection::GlobalPermutation<N>> {
    static constexpr auto name = "global permutation";
};
}  // namespace detail

namespace wrapper {

template <typename Key, typename T, bool Min = true,
          typename QueueSelectionPolicy = typename multiqueue::MultiQueue<
              Key, std::pair<Key, T>,
              std::conditional_t<Min, std::greater<Key>, std::less<Key>>>::traits_type::queue_selection_policy_type>
class MultiQueue : public multiqueue::MultiQueue<Key, std::pair<Key, T>,
                                                 std::conditional_t<Min, std::greater<Key>, std::less<Key>>> {
    using base_type =
        multiqueue::MultiQueue<Key, std::pair<Key, T>, std::conditional_t<Min, std::greater<Key>, std::less<Key>>>;

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
        : base_type(config.factor * num_threads, initial_capacity, config) {
    }

    std::ostream &describe(std::ostream &out) {
        out << "MultiQueue (pqs: " << this->num_pqs()
            << ", queue selection: " << detail::QueueSelectionPolicyTraits<QueueSelectionPolicy>::name
            << ", stickiness: " << this->get_queue_selection_config().stickiness << ')';
        return out;
    }
};

template <typename Key, typename T, bool Min, std::size_t N>
class MultiQueue<Key, T, Min, multiqueue::queue_selection::Random<N>>
    : public multiqueue::MultiQueue<Key, std::pair<Key, T>,
                                    std::conditional_t<Min, std::greater<Key>, std::less<Key>>> {
    using base_type =
        multiqueue::MultiQueue<Key, std::pair<Key, T>, std::conditional_t<Min, std::greater<Key>, std::less<Key>>>;

   public:
    struct config_type : base_type::queue_selection_config_type {
        int factor = 2;
    };

    static void add_options(cxxopts::Options &options, config_type &config) {
        options.add_options()("c,factor", "The number of queues per thread",
                              cxxopts::value<int>(config.factor), "NUMBER");
    }

    MultiQueue(int num_threads, std::size_t initial_capacity, config_type const &config)
        : base_type(config.factor * num_threads, initial_capacity, config) {
    }

    std::ostream &describe(std::ostream &out) {
        out << "MultiQueue (pqs: " << this->num_pqs() << ", queue selection: random)";
        return out;
    }
};

}  // namespace wrapper
