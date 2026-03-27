#ifndef RNAALIGNREFACTORED_MSPC_H
#define RNAALIGNREFACTORED_MSPC_H
#include <atomic>
#include <vector>

/**
 * Each thread buffers its reads before attempting to write.
 * We do this to avoid contention among threads trying to race to
 * write to the output, since we want keep the queue to the writer
 * thread lock-free.
 * 
 * The type `T` must implement a `reserve` subroutine that allocates
 * the needed memory from heap.
 */
template <typename T> struct QueueItem {
    T data;
    std::atomic<QueueItem*> next;

    QueueItem() {
        /**
         * Pre-allocate a chunk of memory immediately upon creation.
         * This has to be a pretty significant amount of data (say, a few MB) 
         * for the queue to be IO efficient.
         * 
         * TODO: (ARVIN) Pull this number out. It should ideally just
         *       replace the output buffer size that we had previously.
         */
        data.reserve(8ULL * 1024 * 1024); 
        next.store(nullptr, std::memory_order_relaxed);
    }
};

/**
 * We allocate buffers once at startup.
 * This pool of buffers will be used to pass hot buffers to each thread that
 * is ready to produce data.
 * In our case, this won't be entirely allocation-free while the algorithm is
 * running, because the type `T` will be `std::string` and can be extended if
 * a read alignment turns out to be really really long.
 * This should happen rarely if the buffers are large enough, and it also means
 * that we don't have to bother resizing the buffer for longer alignments 
 * manually.
 */
template <typename T> class BufferPool {
private:
    std::vector<QueueItem<T>*> pool;
    std::atomic_flag lock = ATOMIC_FLAG_INIT;

public:
    BufferPool(size_t numBuffers) {
        pool.reserve(numBuffers);
        for (size_t i = 0; i < numBuffers; ++i) {
            pool.push_back(new QueueItem<T>());
        }
    }

    QueueItem<T>* pop() {
        // Spin until we acquire the lock
        while (lock.test_and_set(std::memory_order_acquire)) {
            // CPU hint to yield execution slightly during spin
            __builtin_ia32_pause(); 
        }

        QueueItem<T>* node = nullptr;
        if (!pool.empty()) {
            node = pool.back();
            pool.pop_back();
        }

        // Release the lock
        lock.clear(std::memory_order_release);
        return node;
    }

    void push(QueueItem<T>* node) {
        while (lock.test_and_set(std::memory_order_acquire)) {
            __builtin_ia32_pause();
        }
        pool.push_back(node);
        lock.clear(std::memory_order_release);
    }
};

/**
 * A bare-bones multi-producer, single-consumer, lock-free
 * queue. The aligner threads will produce and only a single
 * thread will write to the output.
 * 
 * Notes
 * -----
 * Each producer (aligner) thread must buffer its data until
 * it reaches a significant enough size in order to make this
 * queue efficient. The consumer thread here will be concerned
 * with IO (mainly, populating the output SAM file).
 * A single thread won't be enough to saturate the write throughput
 * if writes are too small (even if completely sequential).
 * For our output SAM file, we expect chunks of a few MB of data ...
 */
template<typename T> class MPSCQueue {
private:
    std::atomic<QueueItem<T>*> head;
    QueueItem<T>* tail;
    QueueItem<T> stub; // Sentinel node to simplify queue logic

public:
    MPSCQueue() {
        stub.next.store(nullptr, std::memory_order_relaxed);
        head.store(&stub, std::memory_order_relaxed);
        tail = &stub;
    }

    void push(QueueItem<T>* node) {
        node->next.store(nullptr, std::memory_order_relaxed);
        
        // Atomically make this node the new head of the list.
        // Returns the _OLD_ head so we can link it.
        QueueItem<T>* prev = head.exchange(node, std::memory_order_acq_rel);

        // Link the old head to the new node
        prev->next.store(node, std::memory_order_release);
    }

    QueueItem<T>* pop() {
        QueueItem<T>* tail_node = tail;
        QueueItem<T>* next = tail_node->next.load(std::memory_order_acquire);

        if (tail_node == &stub) {
            if (!next) return nullptr; // Queue is empty
            tail = next;
            tail_node = next;
            next = next->next.load(std::memory_order_acquire);
        }

        if (next) {
            tail = next;
            return tail_node;
        }

        QueueItem<T>* head_node = head.load(std::memory_order_acquire);
        if (tail_node != head_node) {
            // This may happen if a producer is in the middle of a push (between exchange and store)
            return nullptr;
        }

        push(&stub);
        next = tail_node->next.load(std::memory_order_acquire);
        
        if (next) {
            tail = next;
            return tail_node;
        }

        return nullptr;
    }
};
#endif // RNAALIGNREFACTORED_MSPC_H