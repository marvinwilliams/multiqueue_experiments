#pragma once

#include <algorithm>
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
        prefix_sum_ = {Item{0, 0}};
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
};
