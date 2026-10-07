#pragma once

/**
******************************************************************************
* @file:   merge_heap.hpp
*
* @author: Marvin Williams
* @date:   2021/03/02 16:21
* @brief:
*******************************************************************************
**/

#include <algorithm>
#include <array>
#include <cassert>
#include <iterator>
#include <memory>       // allocator
#include <type_traits>  // is_constructible, enable_if
#include <utility>      // move, forward
#include <vector>

namespace multiqueue {

template <typename T, typename Key, typename KeyExtractor, typename Comparator>
struct heap_base : private KeyExtractor, private Comparator {
    using value_type = T;
    using key_type = Key;
    using key_extractor = KeyExtractor;
    using comp_type = Comparator;
    using reference = value_type &;
    using const_reference = value_type const &;

    static_assert(std::is_invocable_r_v<key_type const &, key_extractor const &, value_type const &>,
                  "Keys must be extractable from values using the signature "
                  "`Key const& KeyExtractor(Value const&) const &`");
    static_assert(std::is_invocable_r_v<bool, comp_type const &, key_type const &, key_type const &>,
                  "Keys must be comparable using the signature `bool "
                  "Comparator(Key const&, Key const&) const &`");
    static_assert(std::is_default_constructible_v<key_extractor>, "`KeyExtractor` must be default-constructible");

   protected:
    static constexpr bool is_key_extract_noexcept = noexcept(std::declval<key_extractor>()(std::declval<value_type>()));
    static constexpr bool is_compare_noexcept =
        noexcept(std::declval<comp_type>()(std::declval<key_type>(), std::declval<key_type>()));

    heap_base() = default;

    explicit heap_base(comp_type const &comp) noexcept(std::is_nothrow_default_constructible_v<key_extractor>)
        : key_extractor(), comp_type(comp) {
    }

    constexpr key_extractor const &to_key_extractor() const noexcept {
        return *this;
    }

    constexpr comp_type const &to_comparator() const noexcept {
        return *this;
    }

    constexpr key_type const &extract_key(value_type const &value) const noexcept(is_key_extract_noexcept) {
        return to_key_extractor()(value);
    }

    constexpr bool compare(key_type const &lhs, key_type const &rhs) const noexcept(is_compare_noexcept) {
        return to_comparator()(lhs, rhs);
    }

    constexpr bool value_compare(value_type const &lhs, value_type const &rhs) const noexcept(is_compare_noexcept) {
        return compare(extract_key(lhs), extract_key(rhs));
    }
};

template <typename InputIt, typename InOutIt, typename OutputIt, typename Comparator>
void inplace_merge(InputIt input, InOutIt in_out, OutputIt output, std::size_t n, Comparator comp) {
    auto const input_end = input + static_cast<typename std::iterator_traits<InputIt>::difference_type>(n);
    auto const in_out_end = in_out + static_cast<typename std::iterator_traits<InOutIt>::difference_type>(n);
    auto in_out_copy = in_out;
    while (n > 0) {
        if (comp(*input, *in_out)) {
            *output++ = std::move(*(input++));
        } else {
            *output++ = std::move(*(in_out++));
        }
        --n;
    }
    if (input != input_end && in_out != in_out_end) {
        while (true) {
            if (comp(*input, *in_out)) {
                *in_out_copy++ = std::move(*(input++));
                if (input == input_end) {
                    assert(in_out == in_out_copy);
                    return;
                }
            } else {
                *in_out_copy++ = std::move(*(in_out++));
                if (in_out == in_out_end) {
                    break;
                }
            }
        }
    }
    in_out_copy = std::move(input, input_end, in_out_copy);
    assert(in_out == in_out_copy);
}

template <typename T, typename Key, typename KeyExtractor, typename Comparator, std::size_t NodeSize,
          typename Allocator = std::allocator<T>>
class merge_heap : private heap_base<T, Key, KeyExtractor, Comparator> {
    using base_type = heap_base<T, Key, KeyExtractor, Comparator>;
    using base_type::compare;
    using base_type::extract_key;
    using base_type::value_compare;

   public:
    using value_type = typename base_type::value_type;
    using key_type = typename base_type::key_type;
    using key_extractor = typename base_type::key_extractor;
    using comp_type = typename base_type::comp_type;
    using reference = typename base_type::reference;
    using const_reference = typename base_type::const_reference;

    using node_type = std::array<value_type, NodeSize>;
    using allocator_type = Allocator;
    using container_type = std::vector<node_type>;
    using iterator = typename container_type::const_iterator;
    using const_iterator = typename container_type::const_iterator;
    using difference_type = typename container_type::difference_type;
    using size_type = std::size_t;

    static_assert(NodeSize > 0 && (NodeSize & (NodeSize - 1)) == 0,
                  "NodeSize must be greater than 0 and a power of two");

   private:
    static constexpr auto node_size_ = NodeSize;
    node_type push_buffer_;
    std::array<value_type, 2 * NodeSize> pop_buffer_;
    std::size_t push_buffer_size_ = 0;
    std::size_t pop_buffer_size_ = 0;
    container_type data_;

   private:
    static constexpr size_type parent_index(size_type const index) noexcept {
        assert(index > 0);
        return (index - 1) >> 1;
    }

    static constexpr size_type first_child_index(size_type const index) noexcept {
        return (index << 1) + 1;
    }

    constexpr bool compare_last(size_type const lhs_index, size_type const rhs_index) const
        noexcept(base_type::is_key_extract_noexcept && base_type::is_compare_noexcept) {
        return value_compare(data_[lhs_index].back(), data_[rhs_index].back());
    }

    size_type min_child_index(size_type index) const noexcept(noexcept(compare_last(0, 0))) {
        assert(index < data_.size());
        index = first_child_index(index);
        assert(index + 1 < data_.size());
        return compare_last(index, index + 1) ? index : index + 1;
    }

    void insert_push_buffer(value_type const &value) {
        if (push_buffer_size_ == NodeSize) {
            insert(push_buffer_.begin(), push_buffer_.end());
            push_buffer_size_ = 0;
        }
        if (push_buffer_size_ == 0) {
            push_buffer_.back() = value;
            ++push_buffer_size_;
            return;
        }
        if (value_compare(push_buffer_.back(), value)) {
            std::move(push_buffer_.end() - push_buffer_size_, push_buffer_.end(),
                      push_buffer_.end() - push_buffer_size_ - 1);
            push_buffer_.back() = value;
            ++push_buffer_size_;
        } else {
            auto pos = push_buffer_.end() - push_buffer_size_;
            while (value_compare(*pos, value)) {
                *(pos - 1) = std::move(*pos);
                ++pos;
            }
            *(pos - 1) = value;
            ++push_buffer_size_;
        }
    }

#ifndef NDEBUG
    bool is_heap() const {
        if (empty()) {
            return true;
        };
        auto value_comparator = [&](const_reference lhs, const_reference rhs) { return value_compare(lhs, rhs); };
        if (!std::is_sorted(data_.front().begin(), data_.front().end(), value_comparator)) {
            return false;
        }
        for (size_type i = 0; i < data_.size(); ++i) {
            for (auto j = first_child_index(i); j < first_child_index(i) + 2u; ++j) {
                if (j >= data_.size()) {
                    return true;
                }
                if (!std::is_sorted(data_[j].begin(), data_[j].end(), value_comparator)) {
                    return false;
                }
                if (value_compare(data_[j].front(), data_[i].back())) {
                    return false;
                }
            }
        }
        return true;
    }
#endif

   public:
    merge_heap() = default;

    explicit merge_heap(allocator_type const &alloc) noexcept(std::is_nothrow_constructible_v<base_type>)
        : base_type(), data_(alloc) {
    }

    explicit merge_heap(comp_type const &comp, allocator_type const &alloc = allocator_type()) noexcept(
        std::is_nothrow_constructible_v<base_type, comp_type>)
        : base_type(comp), data_(alloc) {
    }

    constexpr comp_type const &get_comparator() const noexcept {
        return base_type::to_comparator();
    }

    inline iterator begin() const noexcept {
        return data_.cbegin();
    }

    inline const_iterator cbegin() const noexcept {
        return data_.cbegin();
    }

    inline iterator end() const noexcept {
        return data_.cend();
    }

    inline const_iterator cend() const noexcept {
        return data_.cend();
    }

    [[nodiscard]] inline bool empty() const noexcept {
        return pop_buffer_size_ == 0;
    }

    inline size_type size() const noexcept {
        return data_.size() * NodeSize + pop_buffer_size_ + push_buffer_size_;
    }

    inline const_reference top() const {
        assert(!empty());
        return *(pop_buffer_.end() - pop_buffer_size_);
    }

    inline node_type const &top_node() const {
        assert(!empty());
        return data_.front();
    }

    inline node_type &top_node() {
        assert(!empty());
        return data_.front();
    }

    inline void reserve(std::size_t const cap) {
        data_.reserve(cap / NodeSize + (cap % NodeSize == 0 ? 0 : 1));
    }

    inline void reserve_and_touch(std::size_t const cap) {
        auto const num_nodes = cap / NodeSize + (cap % NodeSize == 0 ? 0 : 1);
        if (data_.size() < num_nodes) {
            auto const old_size = data_.size();
            data_.resize(num_nodes);
            data_.resize(old_size);
        }
    }

    void pop_node() {
        assert(!empty());
        auto value_comparator = [this](const_reference lhs, const_reference rhs) { return value_compare(lhs, rhs); };
        size_type index = 0;
        size_type const first_incomplete_parent = parent_index(data_.size());
        while (index < first_incomplete_parent) {
            auto min_child = first_child_index(index);
            auto max_child = min_child + 1;
            assert(max_child < data_.size());
            if (compare_last(max_child, min_child)) {
                std::swap(min_child, max_child);
            }
            inplace_merge(data_[min_child].begin(), data_[max_child].begin(), data_[index].begin(), NodeSize,
                          value_comparator);
            index = min_child;
        }
        if (first_child_index(index) + 1 == data_.size()) {
            std::move(data_.back().begin(), data_.back().end(), data_[index].begin());
        } else if (index + 1 < data_.size()) {
            auto reverse_value_comparator = [this](const_reference lhs, const_reference rhs) {
                return value_compare(rhs, lhs);
            };

            size_type parent;
            while (index > 0 &&
                   (parent = parent_index(index), value_compare(data_.back().front(), data_[parent].back()))) {
                inplace_merge(data_[parent].rbegin(), data_.back().rbegin(), data_[index].rbegin(), NodeSize,
                              reverse_value_comparator);
                index = parent;
            }
            std::move(data_.back().begin(), data_.back().end(), data_[index].begin());
        }
        data_.pop_back();
        assert(is_heap());
    }

    template <typename Iter>
    void extract_top_node(Iter output) {
        assert(!empty());
        std::move(data_.front().begin(), data_.front().end(), output);
        pop_node();
    }

    void push(value_type const &value) {
        if ((data_.empty() || value_compare(value, data_.front().front())) &&
            (push_buffer_size_ == 0 || value_compare(value, *(push_buffer_.end() - push_buffer_size_)))) {
            if (pop_buffer_size_ == 0) {
                pop_buffer_.back() = value;
                ++pop_buffer_size_;
                return;
            }
            if (value_compare(value, pop_buffer_.back())) {
                auto pos = pop_buffer_.end() - pop_buffer_size_;
                while (value_compare(*pos, value)) {
                    ++pos;
                }
                if (pop_buffer_size_ == 2 * NodeSize) {
                    insert_push_buffer(pop_buffer_.back());
                    std::move_backward(pos, pop_buffer_.end() - 1, pop_buffer_.end());
                    *pos = value;
                } else {
                    std::move(pop_buffer_.end() - pop_buffer_size_, pos, pop_buffer_.end() - pop_buffer_size_ - 1);
                    *(pos - 1) = value;
                    ++pop_buffer_size_;
                }
                return;
            } else if (pop_buffer_size_ < 2 * NodeSize) {
                std::move(pop_buffer_.end() - pop_buffer_size_, pop_buffer_.end(),
                          pop_buffer_.end() - pop_buffer_size_ - 1);
                pop_buffer_.back() = value;
                ++pop_buffer_size_;
                return;
            }
        }
        insert_push_buffer(value);
    }

    void pop() {
        assert(pop_buffer_size_ > 0);
        --pop_buffer_size_;
        if (pop_buffer_size_ == 0) {
            if (data_.empty()) {
                if (push_buffer_size_ > 0) {
                    std::move(push_buffer_.end() - push_buffer_size_, push_buffer_.end(),
                              pop_buffer_.end() - push_buffer_size_);
                    pop_buffer_size_ = push_buffer_size_;
                    push_buffer_size_ = 0;
                }
                return;
            }
            auto push_buffer_it = push_buffer_.end();
            auto top_node_it = data_.front().end();
            auto pop_buffer_it = pop_buffer_.end();
            auto const push_buffer_begin = push_buffer_.end() - push_buffer_size_;
            while (push_buffer_it != push_buffer_begin && value_compare(*(top_node_it - 1), *(push_buffer_it - 1))) {
                --push_buffer_it;
            }
            push_buffer_size_ = static_cast<std::size_t>(std::distance(push_buffer_it, push_buffer_.end()));
            while (top_node_it != data_.front().begin() && push_buffer_it != push_buffer_begin) {
                if (value_compare(*(top_node_it - 1), *(push_buffer_it - 1))) {
                    *(--pop_buffer_it) = std::move(*(--push_buffer_it));
                } else {
                    *(--pop_buffer_it) = std::move(*(--top_node_it));
                }
            }
            while (top_node_it != data_.front().begin()) {
                *(--pop_buffer_it) = std::move(*(--top_node_it));
            }
            while (push_buffer_it != push_buffer_begin) {
                *(--pop_buffer_it) = std::move(*(--push_buffer_it));
            }
            pop_buffer_size_ = static_cast<std::size_t>(std::distance(pop_buffer_it, pop_buffer_.end()));
            pop_node();
        }
    }

    template <typename Iter>
    void insert(Iter first, Iter last) {
        auto reverse_value_comparator = [this](const_reference lhs, const_reference rhs) {
            return value_compare(rhs, lhs);
        };
        auto index = data_.size();
        data_.push_back({});
        while (index > 0) {
            auto const parent = parent_index(index);
            assert(parent < index);
            if (!value_compare(*first, data_[parent].back())) {
                break;
            }
            inplace_merge(data_[parent].rbegin(), std::reverse_iterator{last}, data_[index].rbegin(), NodeSize,
                          reverse_value_comparator);
            index = parent;
        }
        std::move(first, last, data_[index].begin());
        assert(is_heap());
    }

    inline void clear() noexcept {
        pop_buffer_size_ = 0;
        push_buffer_size_ = 0;
        data_.clear();
    }
};

struct identity {
    template <typename T>
    constexpr T const &operator()(T const &t) const noexcept {
        return t;
    }
};

template <typename T, typename Comparator = std::less<T>, std::size_t NodeSize = 64,
          typename Allocator = std::allocator<T>>
using value_merge_heap = merge_heap<T, T, identity, Comparator, NodeSize, Allocator>;

}  // namespace multiqueue
