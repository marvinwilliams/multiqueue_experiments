#ifndef GALOIS_STEALINGMULTIQUEUE_H
#define GALOIS_STEALINGMULTIQUEUE_H

// adapted to standalone from https://github.com/npostnikova/mq-based-schedulers

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <optional>
#include <vector>

namespace Galois {

template <typename T, std::size_t Alignment = 128>
struct alignas(Alignment) AlignedObject {
    T data;
};

namespace WorkList {

/**
 * Class-helper, consists of a sequential heap and
 * a stealing buffer.
 *
 * @tparam T Type of the elements.
 * @tparam Compare Elements comparator.
 * @tparam STEAL_NUM Number of elements to steal at once.
 * @tparam D Arity of the heap.
 */
template <typename T, typename Compare, std::size_t STEAL_NUM, std::size_t D = 4>
class HeapWithStealBuffer {
    using value_type = T;
    using index_t = std::size_t;
    // Local priority queue.
    std::vector<T> heap;
    // Other threads steal the whole buffer at once.
    std::array<T, STEAL_NUM> stealBuffer;
    // Represents epoch & stolen flag
    // version mod 2 = 0  -- elements are stolen
    // version mod 2 = 1  -- can steal
    std::atomic<std::size_t> version;

   public:
    // Represents a flag for empty buffer cells.
    static T dummy;
    // Comparator.
    Compare compare;

    HeapWithStealBuffer() : version(0) {
        for (std::size_t i = 0; i < STEAL_NUM; i++) {
            stealBuffer[i] = dummy;
        }
    }

    //! Checks whether the element is "null".
    static bool isDummy(T const& element) {
        return element == dummy;
    }

    //! Gets current version of the stealing buffer.
    std::size_t getVersion() {
        return version.load(std::memory_order_acquire);
    }

    //! Checks whether elements in the buffer are stolen.
    std::size_t isBufferStolen() {
        return getVersion() % 2 == 0;
    }

    //! Fills stealing buffer if the current tasks are stolen.
    void fillBufferIfStolen() {
        if (isBufferStolen()) {
            fillBuffer();
        }
    }

    //! Get min among elements that can be stolen.
    //! Sets a flag to true, if operation failed because of a race.
    T getBufferMin(bool& raceHappened) {
        auto v1 = getVersion();
        if (v1 % 2 == 0) {
            return dummy;
        }
        T minVal = stealBuffer[0];
        auto v2 = getVersion();
        if (v1 == v2) {
            return minVal;
        }
        // Somebody has stolen the elements.
        raceHappened = true;
        return dummy;
    }

    //! Returns min element from the buffer, updating the buffer if empty.
    //! Can be called only by the thread-owner.
    T getMinWriter() {
        auto v1 = getVersion();
        if (v1 % 2 != 0) {
            T minVal = stealBuffer[0];
            auto v2 = getVersion();
            if (v1 == v2) {
                return minVal;
            }
        }
        return fillBuffer();
    }

    //! Fills the steal buffer.
    //! Called when the elements from the previous epoch are empty.
    T fillBuffer() {
        if (heap.empty())
            return dummy;
        std::array<T, STEAL_NUM> elements;
        elements.fill(dummy);
        for (std::size_t i = 0; i < STEAL_NUM && !heap.empty(); i++) {
            elements[i] = popLocally();
        }
        stealBuffer = elements;
        version.fetch_add(1, std::memory_order_acq_rel);
        return elements[0];
    }

    //! Tries to steal the elements from the stealing buffer.
    std::optional<std::array<T, STEAL_NUM>> trySteal(bool& raceHappened) {
        auto emptyRes = std::optional<std::array<T, STEAL_NUM>>();
        auto v1 = getVersion();
        if (v1 % 2 == 0) {
            // Already stolen.
            return emptyRes;
        }
        std::array<T, STEAL_NUM> buffer = stealBuffer;
        if (version.compare_exchange_weak(v1, v1 + 1, std::memory_order_acq_rel)) {
            return buffer;
        }
        // Another thread got ahead.
        raceHappened = true;
        return emptyRes;
    }

    //! Retrieves an element from the heap.
    T popLocally() {
        auto res = heap[0];
        heap[0] = heap.back();
        heap.pop_back();
        if (heap.size() > 0) {
            sift_down(0);
        }
        return res;
    }

    //! Extract min from the structure: both the buffer and the heap
    //! are considered. Called from the owner-thread.
  std::optional<value_type> extractMin() {
        if (heap.empty()) {
            // Only check the steal buffer.
            return tryStealLocally();
        }
        bool raceFlag = false;  // useless now
        auto bufferMin = getBufferMin(raceFlag);
        if (!isDummy(bufferMin) && compare(heap[0], bufferMin)) {
            auto stolen = tryStealLocally();
            if (stolen) {
                fillBuffer();
                return stolen;
            }
        }
        auto localMin = popLocally();
        if (isDummy(bufferMin))
            fillBuffer();
        return localMin;
    }

    //! Inserts the element into the heap.
    void pushLocally(T const& val) {
        index_t index = heap.size();
        heap.push_back({val});
        sift_up(index);
    }

   private:
    //! Tries to steal elements from local buffer.
    //! Return minimum among stolen elements.
    std::optional<T> tryStealLocally() {
        bool raceFlag = false;  // useless now
        auto stolen = trySteal(raceFlag);
        if (stolen) {
            auto elements = *stolen;
            for (std::size_t i = 1; i < STEAL_NUM && !isDummy(elements[i]); i++) {
                pushLocally(elements[i]);
            }
            return elements[0];
        }
        return std::optional<T>();
    }

    ///////////////////////// HEAP /////////////////////////
    void swap(index_t i, index_t j) {
        T t = heap[i];
        heap[i] = heap[j];
        heap[j] = t;
    }

    //! Check whether the index of the root passed.
    bool is_root(index_t index) {
        return index == 0;
    }

    //! Check whether the index is not out of bounds.
    bool is_valid_index(index_t index) {
        return index >= 0 && index < heap.size();
    }

    //! Get index of the parent.
    std::optional<index_t> get_parent(index_t index) {
        if (!is_root(index) && is_valid_index(index)) {
            return (index - 1) / D;
        }
        return {};
    }

    //! Get index of the smallest (due `Comparator`) child.
    std::optional<index_t> get_smallest_child(index_t index) {
        if (!is_valid_index(D * index + 1)) {
            return {};
        }
        index_t smallest = D * index + 1;
        for (std::size_t k = 2; k <= D; k++) {
            index_t k_child = D * index + k;
            if (!is_valid_index(k_child))
                break;
            if (compare(heap[smallest], heap[k_child]))
                smallest = k_child;
        }
        return smallest;
    }

    //! Sift down without decrease key info update.
    void sift_down(index_t index) {
        auto smallest_child = get_smallest_child(index);
        while (smallest_child && compare(heap[index], heap[*smallest_child])) {
            swap(index, *smallest_child);
            index = *smallest_child;
            smallest_child = get_smallest_child(index);
        }
    }

    //! Sift up the element with provided index.
    index_t sift_up(index_t index) {
        std::optional<index_t> parent = get_parent(index);

        while (parent && compare(heap[*parent], heap[index])) {
            swap(index, *parent);
            index = *parent;
            parent = get_parent(index);
        }
        return index;
    }
};

template <typename T, typename Compare, std::size_t STEAL_NUM, std::size_t D>
T HeapWithStealBuffer<T, Compare, STEAL_NUM, D>::dummy;

template <typename T, typename Comparer, size_t StealProb, size_t StealBatchSize, bool Concurrent = true>
class StealingMultiQueue {
   private:
    using Heap = HeapWithStealBuffer<T, Comparer, StealBatchSize, 4>;
    std::unique_ptr<AlignedObject<Heap>[]> heaps;
    std::unique_ptr<AlignedObject<std::vector<T>>[]> stealBuffers;
    Comparer compare;
    const size_t nQ;

    //! Thread local random.
    uint32_t random() {
        static thread_local uint32_t x = std::chrono::system_clock::now().time_since_epoch().count() % 16386 + 1;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        return x;
    }

    //! Index of a random heap.
    size_t rand_heap() {
        return random() % nQ;
    }

    //! Tries to steal from a random queue.
    //! Repeats if failed because of a race.
    std::optional<T> trySteal(int tId) {
        T localMin = heaps[tId].data.getMinWriter();
        bool nextIterNeeded = true;
        while (nextIterNeeded) {
            auto randId = rand_heap();
            if (randId == tId)
                continue;
            nextIterNeeded = false;
            Heap* randH = &heaps[randId].data;
            auto randMin = randH->getBufferMin(nextIterNeeded);
            if (randH->isDummy(randMin)) {
                // Nothing to steal.
                continue;
            }
            if (Heap::isDummy(localMin) || compare(localMin, randMin)) {
                auto stolen = randH->trySteal(nextIterNeeded);
                if (stolen) {
                    auto& buffer = stealBuffers[tId].data;
                    auto elements = *stolen;
                    for (size_t i = 1; i < elements.size() && !Heap::isDummy(elements[i]); i++) {
                        buffer.push_back(elements[i]);
                    }
                    std::reverse(buffer.begin(), buffer.end());
                    return elements[0];
                }
            }
        }
        return std::nullopt;
    }

   public:
    StealingMultiQueue(int num_threads) : nQ(num_threads) {
        std::memset(reinterpret_cast<void*>(&Heap::dummy), 0xff, sizeof(Heap::dummy));
        heaps = std::make_unique<AlignedObject<Heap>[]>(nQ);
        stealBuffers = std::make_unique<AlignedObject<std::vector<T>>[]>(nQ);
    }

    using value_type = T;

    //! Change the concurrency flag.
    template <bool _concurrent>
    struct rethread {
        using type = StealingMultiQueue<T, Comparer, StealProb, StealBatchSize, _concurrent>;
    };

    //! Change the type the worklist holds.
    template <typename T_>
    struct retype {
        using type = StealingMultiQueue<T_, Comparer, StealProb, StealBatchSize, Concurrent>;
    };

    template <typename RangeTy>
    unsigned int push_initial(const RangeTy& range) {
        auto rp = range.local_pair();
        return push(rp.first, rp.second);
    }

    template <typename Iter>
    unsigned int push(int tId, Iter b, Iter e) {
        if (b == e)
            return 0;
        unsigned int pushedNum = 0;
        Heap* heap = &heaps[tId].data;
        while (b != e) {
            heap->pushLocally(*b++);
            pushedNum++;
        }
        heap->fillBufferIfStolen();
        return pushedNum;
    }

    std::optional<T> pop(int tId) {
        auto& buffer = stealBuffers[tId].data;
        if (!buffer.empty()) {
            auto val = buffer.back();
            buffer.pop_back();
            heaps[tId].data.fillBufferIfStolen();
            return val;
        }
        // rand == 0 -- try to steal
        // otherwise, pop locally
        if (nQ > 1 && random() % StealProb == 0) {
            std::optional<T> stolen = trySteal(tId);
            if (stolen)
                return stolen;
        }
        auto minVal = heaps[tId].data.extractMin();
        if (minVal)
            return minVal;

        // Our heap is empty.
        return nQ == 1 ? std::nullopt : trySteal(tId);
    }
};

/* GALOIS_WLCOMPILECHECK(StealingMultiQueue) */

}  // namespace WorkList
}  // namespace Galois

#endif  // GALOIS_STEALINGMULTIQUEUE_H
