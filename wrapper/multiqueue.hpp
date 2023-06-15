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

namespace wrapper {

template <typename Key, typename T, bool Min = true,
          typename QueueSelectionPolicy = typename multiqueue::MultiQueue<
              Key, std::pair<Key, T>,
              std::conditional_t<Min, std::greater<Key>, std::less<Key>>>::traits_type::queue_selection_policy_type>
class MultiQueue
    : public multiqueue::MultiQueue<Key, std::pair<Key, T>, std::conditional_t<Min, std::greater<Key>, std::less<Key>>> {
    using base_type =
        multiqueue::MultiQueue<Key, std::pair<Key, T>, std::conditional_t<Min, std::greater<Key>, std::less<Key>>>;
    using config_type = typename base_type::queue_selection_config_type;

   public:
    static void add_options(cxxopts::Options &options) {
        options.add_options()("c,factor", "The multiple of threads used for the queues", cxxopts::value<int>(),
                              "NUMBER");
        options.add_options()("k,stickiness", "The stickiness period", cxxopts::value<int>(), "NUMBER");
    }

    MultiQueue(int num_threads, std::size_t initial_capacity, cxxopts::ParseResult const &options)
        : base_type((options.count("factor") > 0 ? options["factor"].as<int>() : 2) * num_threads, initial_capacity,
                    options.count("stickiness") > 0 ? config_type{1, options["factor"].as<int>()} : config_type{}) {
    }

    std::ostream &describe(std::ostream &out) {
        out << "MultiQueue (pqs: " << this->num_pqs()
            << ", stickiness: " << this->get_queue_selection_config().stickiness << ')';
        return out;
    }
};

template <typename Key, typename T, bool Min, std::size_t N>
class MultiQueue<Key, T, Min, multiqueue::queue_selection::Random<N>>
    : public multiqueue::MultiQueue<Key, std::pair<Key, T>, std::conditional_t<Min, std::greater<Key>, std::less<Key>>> {
    using base_type =
        multiqueue::MultiQueue<Key, std::pair<Key, T>, std::conditional_t<Min, std::greater<Key>, std::less<Key>>>;
    using config_type = typename base_type::queue_selection_config_type;

   public:
    static void add_options(cxxopts::Options &options) {
        options.add_options()("c,factor", "The multiple of threads used for the queues", cxxopts::value<int>(),
                              "NUMBER");
    }

    MultiQueue(int num_threads, std::size_t initial_capacity, cxxopts::ParseResult const &options)
        : base_type((options.count("factor") > 0 ? options["factor"].as<int>() : 2) * num_threads, initial_capacity) {
    }

    std::ostream &describe(std::ostream &out) {
        out << "MultiQueue (pqs: " << this->num_pqs() << ')';
        return out;
    }
};

}  // namespace wrapper
