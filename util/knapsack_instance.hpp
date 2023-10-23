#pragma once

#include <iostream>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <vector>

class KnapsackInstance {
   public:
    struct Item {
        long long weight;
        long long value;
    };

   private:
    std::vector<Item> prefix_sum_;
    long long capacity_{};

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
        prefix_sum_.clear();
        prefix_sum_.reserve(n + 1);
        in >> capacity_;
        if (!in || (n > 0 && in.eof())) {
            throw std::runtime_error{"Could not get capacity"};
        }

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
        std::sort(prefix_sum_.begin(), prefix_sum_.end(), [](auto const& lhs, auto const& rhs) {
            return (static_cast<double>(lhs.value) / static_cast<double>(lhs.weight)) >
                (static_cast<double>(rhs.value) / static_cast<double>(rhs.weight));
        });
        std::inclusive_scan(prefix_sum_.begin(), prefix_sum_.end(), prefix_sum_.begin(),
                            [](auto const& lhs, auto const& rhs) {
                                return Item{lhs.weight + rhs.weight, lhs.value + rhs.value};
                            });
        prefix_sum_.insert(prefix_sum_.begin(), Item{0, 0});
    }

    [[nodiscard]] std::size_t size() const noexcept {
        return prefix_sum_.size() - 1;
    }

    [[nodiscard]] long long capacity() const noexcept {
        return capacity_;
    }

    [[nodiscard]] long long weight(std::size_t index) const noexcept {
        return prefix_sum_[index + 1].weight - prefix_sum_[index].weight;
    }

    [[nodiscard]] long long value(std::size_t index) const noexcept {
        return prefix_sum_[index + 1].value - prefix_sum_[index].value;
    }

    [[nodiscard]] std::pair<long long, long long> compute_bounds_linear(long long capacity,
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

    [[nodiscard]] std::pair<long long, long long> compute_bounds_binary(long long capacity,
                                                                        std::size_t index) const noexcept {
        assert(index < prefix_sum_.size());
        auto it = std::upper_bound(std::next(prefix_sum_.begin(), static_cast<std::ptrdiff_t>(index) + 1),
                                   prefix_sum_.end(), prefix_sum_[index].weight + capacity,
                                   [](auto lhs, auto const& rhs) { return lhs < rhs.weight; });
        long long lower_bound = std::prev(it)->value - prefix_sum_[index].value;
        if (it == prefix_sum_.end()) {
            return {lower_bound, lower_bound};
        }
        auto residual_capacity = capacity - (std::prev(it)->weight - prefix_sum_[index].weight);
        auto fractional_value =
            ((it->value - std::prev(it)->value) * residual_capacity) / (it->weight - std::prev(it)->weight);
        return {lower_bound, lower_bound + fractional_value};
    }
};
