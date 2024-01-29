#pragma once

#include <iostream>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <random>
#include <vector>

template <typename WeightType = long long, typename ValueType = WeightType>
class KnapsackInstance {
    static_assert(std::is_arithmetic_v<WeightType> && std::is_arithmetic_v<ValueType>,
                  "WeightType and ValueType must be arithmetic types");

   public:
    struct Item {
        WeightType weight;
        ValueType value;
    };

   private:
    std::vector<Item> prefix_sum_;
    WeightType capacity_{};

    void to_prefix_sum() {
        std::sort(prefix_sum_.begin() + 1, prefix_sum_.end(), [](auto const& lhs, auto const& rhs) {
            return (static_cast<double>(lhs.value) / static_cast<double>(lhs.weight)) >
                (static_cast<double>(rhs.value) / static_cast<double>(rhs.weight));
        });
        std::inclusive_scan(prefix_sum_.begin() + 1, prefix_sum_.end(), prefix_sum_.begin() + 1,
                            [](auto const& lhs, auto const& rhs) {
                                return Item{lhs.weight + rhs.weight, lhs.value + rhs.value};
                            });
    }

   public:
    KnapsackInstance() = default;
    KnapsackInstance(std::filesystem::path const& file) {
        std::ifstream in{file};
        if (!in) {
            throw std::runtime_error{"Could not open file"};
        }
        std::size_t n{};
        in >> n;
        if (!in || in.eof()) {
            throw std::runtime_error{"Could not get number of items"};
        }
        in >> capacity_;
        if (!in || (n > 0 && in.eof())) {
            throw std::runtime_error{"Could not get capacity"};
        }
        prefix_sum_.reserve(n + 1);
        prefix_sum_.push_back(Item{0, 0});

        for (std::size_t i = 0; i < n; ++i) {
            if (!in || in.eof()) {
                throw std::runtime_error{"Unexpected end of file"};
            }
            auto item = Item{};
            in >> item.value;
            if (!in || in.eof()) {
                throw std::runtime_error{"Could not read item value"};
            }
            in >> item.weight;
            if (!in) {
                throw std::runtime_error{"Could not read item weight"};
            }
            prefix_sum_.push_back(item);
        }
        to_prefix_sum();
    }

    KnapsackInstance(long long n, WeightType min_weight, WeightType max_weight, ValueType min_add, ValueType max_add,
                     double capacity_factor, unsigned int seed)
        : prefix_sum_(static_cast<std::size_t>(n + 1)),
          capacity_(static_cast<WeightType>(static_cast<double>(n) * static_cast<double>(max_weight - min_weight) *
                                            capacity_factor)) {
        std::default_random_engine rng(seed);
        std::generate(prefix_sum_.begin() + 1, prefix_sum_.end(), [&rng, min_weight, max_weight, min_add, max_add]() {
            if constexpr (std::is_floating_point_v<WeightType>) {
                std::uniform_real_distribution<WeightType> random_weight(min_weight, max_weight);
                std::uniform_real_distribution<ValueType> random_add(min_add, max_add);
                auto weight = random_weight(rng);
                return Item{weight, static_cast<ValueType>(weight) + random_add(rng)};
            } else {
                std::uniform_int_distribution<WeightType> random_weight(min_weight, max_weight);
                std::uniform_int_distribution<ValueType> random_add(min_add, max_add);
                auto weight = random_weight(rng);
                return Item{weight, static_cast<ValueType>(weight) + random_add(rng)};
            }
        });
        to_prefix_sum();
    }

    [[nodiscard]] std::size_t size() const noexcept {
        return prefix_sum_.size() - 1;
    }

    [[nodiscard]] WeightType capacity() const noexcept {
        return capacity_;
    }

    [[nodiscard]] WeightType weight(std::size_t index) const noexcept {
        return prefix_sum_[index + 1].weight - prefix_sum_[index].weight;
    }

    [[nodiscard]] ValueType value(std::size_t index) const noexcept {
        return prefix_sum_[index + 1].value - prefix_sum_[index].value;
    }

    [[nodiscard]] std::pair<ValueType, ValueType> compute_bounds_linear(WeightType capacity,
                                                                        std::size_t index) const noexcept {
        assert(index < prefix_sum_.size());
        auto value_offset = prefix_sum_[index].value;
        capacity += prefix_sum_[index].weight;
        while (index < size() && capacity >= prefix_sum_[index + 1].weight) {
            ++index;
        }
        auto lower_bound = prefix_sum_[index].value - value_offset;
        if (index == size() || capacity == prefix_sum_[index].weight) {
            return {lower_bound, lower_bound};
        }
        auto residual_capacity = capacity - prefix_sum_[index].weight;
        auto fractional_value = (value(index) * residual_capacity) / weight(index);
        return {lower_bound, lower_bound + fractional_value};
    }

    [[nodiscard]] std::pair<ValueType, ValueType> compute_bounds_binary(WeightType capacity,
                                                                        std::size_t index) const noexcept {
        assert(index < prefix_sum_.size());
        auto it = std::upper_bound(std::next(prefix_sum_.begin(), static_cast<std::ptrdiff_t>(index) + 1),
                                   prefix_sum_.end(), prefix_sum_[index].weight + capacity,
                                   [](auto lhs, auto const& rhs) { return lhs < rhs.weight; });
        ValueType lower_bound = std::prev(it)->value - prefix_sum_[index].value;
        if (it == prefix_sum_.end()) {
            return {lower_bound, lower_bound};
        }
        auto residual_capacity = capacity - (std::prev(it)->weight - prefix_sum_[index].weight);
        auto fractional_value =
            ((it->value - std::prev(it)->value) * residual_capacity) / (it->weight - std::prev(it)->weight);
        return {lower_bound, lower_bound + fractional_value};
    }
};
