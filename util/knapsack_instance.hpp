#pragma once

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <vector>

struct KnapsackInstance {
    struct Item {
        unsigned long weight;
        unsigned long value;
    };
    std::vector<Item> items;
    unsigned long capacity;
    std::vector<Item> prefix_sum;

    void from_file(std::filesystem::path const& file) {
        std::ifstream in{file};
        if (!in) {
            throw std::runtime_error{"Could not open file"};
        }
        std::size_t n{};
        in >> n;
        if (!in || in.eof()) {
            throw std::runtime_error{"Could not get number of items"};
        }
        items.clear();
        items.reserve(n);
        in >> capacity;
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
            items.push_back(item);
        }
        std::sort(items.begin(), items.end(), [](auto const& lhs, auto const& rhs) {
            return (static_cast<double>(lhs.value) / static_cast<double>(lhs.weight)) >
                (static_cast<double>(rhs.value) / static_cast<double>(rhs.weight));
        });
        prefix_sum.reserve(items.size() + 1);
        prefix_sum = {Item{0, 0}};
        std::inclusive_scan(items.begin(), items.end(), std::back_inserter(prefix_sum),
                            [](auto const& lhs, auto const& rhs) {
                                return Item{lhs.weight + rhs.weight, lhs.value + rhs.value};
                            });
    }
};
