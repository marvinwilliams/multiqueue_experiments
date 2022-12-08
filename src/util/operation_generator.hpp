/**
******************************************************************************
* @file:   operation_generator.hpp
*
* @author: Marvin Williams
* @date:   2021/06/17 10:31
* @brief:
*******************************************************************************
**/
#pragma once
#ifndef OPERATION_GENERATOR_HPP_INCLUDED
#define OPERATION_GENERATOR_HPP_INCLUDED

#include <array>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <limits>
#include <random>

enum class KeyDistribution : std::size_t { Uniform, Ascending, Descending };

constexpr char const* get_key_distribution_name(KeyDistribution distribution) {
    switch (distribution) {
        case KeyDistribution::Uniform:
            return "uniform";
        case KeyDistribution::Ascending:
            return "ascending";
        case KeyDistribution::Descending:
            return "descending";
        default:
            return "unknown";
    }
}

template <typename Key>
struct OperationGenerator {
    KeyDistribution key_distribution;
    Key min_key;
    Key max_key;
};

template <typename Key>
struct KeyGenerator {
    static constexpr unsigned int log_bits = 6;
    static constexpr unsigned int num_bits = 1u << log_bits;
    static_assert(std::numeric_limits<unsigned long long>::digits == num_bits,
                  "num_bits wrong, adjust log_bits");

    KeyDistribution distribution;

    std::bitset<num_bits> random_bits;
    std::uint8_t bit_pos : log_bits;

    std::uniform_int_distribution<unsigned long long> insert_dist;
    std::uniform_int_distribution<Key> key_dist;

    Key current;
    Key limit;

    KeyGenerator() {}

    explicit KeyGenerator(OperationGenerator<Key> const& config)
        : distribution{config.key_distribution},
          random_bits{0},
          bit_pos{0},
          current{0},
          limit{config.max_key} {
        switch (distribution) {
            case KeyDistribution::Uniform:
                key_dist = std::uniform_int_distribution<Key>(config.min_key,
                                                              config.max_key);
                break;
            case KeyDistribution::Ascending:
                current = config.min_key;
                break;
            case KeyDistribution::Descending:
                current = config.max_key;
                limit = config.min_key;
                break;
        }
    }

    template <typename Generator>
    bool insert(Generator& rng) {
        if (bit_pos == 0) {
            random_bits = insert_dist(rng);
        }
        return random_bits[bit_pos++];
    }

    template <typename Generator>
    Key get_key(Generator& rng) {
        switch (distribution) {
            case KeyDistribution::Uniform:
                return key_dist(rng);
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

#endif  //! OPERATION_GENERATOR_HPP_INCLUDED
