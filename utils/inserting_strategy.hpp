/**
******************************************************************************
* @file:   inserting_strategy.hpp
*
* @author: Marvin Williams
* @date:   2021/06/17 10:31
* @brief:
*******************************************************************************
**/
#pragma once
#ifndef INSERTING_STRATEGY_HPP_INCLUDED
#define INSERTING_STRATEGY_HPP_INCLUDED

#include <array>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <random>

enum class InsertPolicy : std::size_t { Uniform, Split, Producer, Alternating };
constexpr char const* get_insert_policy_name(InsertPolicy policy) {
    switch (policy) {
        case InsertPolicy::Uniform:
            return "uniform";
        case InsertPolicy::Split:
            return "split";
        case InsertPolicy::Producer:
            return "producer";
        case InsertPolicy::Alternating:
            return "alternating";
        default:
            return "unknown";
    }
}
enum class KeyDistribution : std::size_t { Uniform, Dijkstra, Ascending, Descending, ThreadId };
constexpr char const* get_key_distribution_name(KeyDistribution distribution) {
    switch (distribution) {
        case KeyDistribution::Uniform:
            return "uniform";
        case KeyDistribution::ThreadId:
            return "threadid";
        case KeyDistribution::Dijkstra:
            return "dijkstra";
        case KeyDistribution::Ascending:
            return "ascending";
        case KeyDistribution::Descending:
            return "descending";
        default:
            return "unknown";
    }
}

template <typename Key>
struct InsertConfig {
    InsertPolicy insert_policy;
    KeyDistribution key_distribution;
    Key min_key;
    Key max_key;
    Key dijkstra_min_increase;
    Key dijkstra_max_increase;
};

template <typename Key>
struct InsertingStrategy {
    static constexpr unsigned int log_bits = 6;
    static constexpr unsigned int num_bits = 1u << log_bits;
    static_assert(std::numeric_limits<unsigned long long>::digits == num_bits, "num_bits wrong, adjust log_bits");

    InsertPolicy const policy;
    KeyDistribution const distribution;

    std::bitset<num_bits> random_bits;
    std::uint8_t bit_pos : 6;

    std::mt19937 gen;
    std::uniform_int_distribution<unsigned long long> insert_dist;
    std::uniform_int_distribution<Key> key_dist;

    Key current;
    Key limit;

    explicit InsertingStrategy(unsigned int id, InsertConfig<Key> const& config, std::uint32_t seed)
        : policy{config.insert_policy},
          distribution{config.key_distribution},
          random_bits{0},
          bit_pos{0},
          current{0},
          limit{config.max_key} {
        std::seed_seq seq{seed};
        gen.seed(seq);
        switch (policy) {
            case InsertPolicy::Split:
            case InsertPolicy::Alternating:
                random_bits[0] = (id % 2 == 0);
                break;
            case InsertPolicy::Producer:
                random_bits[0] = (id == 0);
                break;
            default:;
        }
        switch (distribution) {
            case KeyDistribution::Uniform:
                key_dist = std::uniform_int_distribution<Key>(config.min_key, config.max_key);
                break;
            case KeyDistribution::Ascending:
                current = config.min_key;
                break;
            case KeyDistribution::Descending:
                current = config.max_key;
                limit = config.min_key;
                break;
            case KeyDistribution::Dijkstra:
                current = config.min_key;
                key_dist =
                    std::uniform_int_distribution<Key>(config.dijkstra_min_increase, config.dijkstra_max_increase);
                break;
            case KeyDistribution::ThreadId:
                key_dist = std::uniform_int_distribution<Key>(id * 1000, (id + 2) * 1000);
                break;
        }
    }

    bool insert() {
        switch (policy) {
            case InsertPolicy::Uniform: {
                if (bit_pos == 0) {
                    random_bits = insert_dist(gen);
                }
                return random_bits[bit_pos++];
            }
            case InsertPolicy::Split:
            case InsertPolicy::Producer:
                return random_bits[0];
            case InsertPolicy::Alternating:
                return random_bits[0].flip();
            default:
                assert(false);
                return false;
        }
    }

    Key get_key() {
        switch (distribution) {
            case KeyDistribution::Uniform:
            case KeyDistribution::ThreadId:
                return key_dist(gen);
            case KeyDistribution::Dijkstra:
                return current + key_dist.b() <= limit ? current++ + key_dist(gen) : limit;
            case KeyDistribution::Ascending:
                return current < limit ? current++ : limit;
            case KeyDistribution::Descending:
                return current > limit ? current-- : limit;
            default:
                assert(false);
                return Key{};
        }
    }
};

#endif  //! INSERTING_STRATEGY_HPP_INCLUDED
