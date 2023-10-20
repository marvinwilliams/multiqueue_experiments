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
    std::vector<Item> items_;
    long long capacity_{};
    std::vector<Item> prefix_sum_;

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
        items_.clear();
        items_.reserve(n);
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
            items_.push_back(item);
        }
        std::sort(items_.begin(), items_.end(), [](auto const& lhs, auto const& rhs) {
            return (static_cast<double>(lhs.value) / static_cast<double>(lhs.weight)) >
                (static_cast<double>(rhs.value) / static_cast<double>(rhs.weight));
        });
        prefix_sum_.reserve(items_.size() + 1);
        prefix_sum_.push_back(Item{0, 0});
        std::inclusive_scan(items_.begin(), items_.end(), std::back_inserter(prefix_sum_),
                            [](auto const& lhs, auto const& rhs) {
                                return Item{lhs.weight + rhs.weight, lhs.value + rhs.value};
                            });
    }

    [[nodiscard]] std::size_t size() const noexcept {
        return items_.size();
    }
    [[nodiscard]] long long capacity() const noexcept {
        return capacity_;
    }
    [[nodiscard]] std::vector<Item> const& items() const noexcept {
        return items_;
    }
    [[nodiscard]] std::vector<Item> const& prefix_sum() const noexcept {
        return prefix_sum_;
    }

    [[nodiscard]] std::pair<long long, long long> compute_bounds_linear(long long capacity,
                                                                        std::size_t index) const noexcept {
        assert(index <= items_.size());
        long long lower_bound = 0;
        while (index != items_.size() && capacity >= items_[index].weight) {
            capacity -= items_[index].weight;
            lower_bound += items_[index].value;
            ++index;
        }
        if (index == items_.size() || capacity == 0) {
            return {lower_bound, lower_bound};
        }
        auto fractional_value = (items_[index].value * capacity) / items_[index].weight;
        return {lower_bound, lower_bound + fractional_value};
    }

    [[nodiscard]] std::pair<long long, long long> compute_bounds_hint(long long capacity, std::size_t index,
                                                                      std::size_t& hint) const noexcept {
        assert(index <= items_.size());
        capacity += prefix_sum_[index].weight;
        while (hint != prefix_sum_.size() && capacity >= prefix_sum_[hint].weight) {
            ++hint;
        }
        auto lower_bound = prefix_sum_[hint - 1].value - prefix_sum_[index].value;
        if (hint == prefix_sum_.size() || prefix_sum_[hint - 1].weight == capacity) {
            return {lower_bound, lower_bound};
        }
        auto residual_capacity = capacity - prefix_sum_[hint - 1].weight;
        auto fractional_value = ((prefix_sum_[hint].value - prefix_sum_[hint - 1].value) * residual_capacity) /
            (prefix_sum_[hint].weight - prefix_sum_[hint - 1].weight);
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
