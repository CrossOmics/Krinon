#ifndef RNAALIGNMENT_ATOMIC_GENERATION_H
#define RNAALIGNMENT_ATOMIC_GENERATION_H
#include <mutex>
#include <cassert>

namespace rna
{
    /**
     * A latch with a fixed size that atomically resets to its original 
     * value when it hits zero.
     */
    class AtomicGenerationCounter {
    private:
        int size_;
        int value_;
        std::mutex lock_;
    public:
        AtomicGenerationCounter(int size) :
            size_(size),
            value_(size)
        {
            assert(size > 0);
        }

        bool decrement() {
            bool reset = false;
            lock_.lock();
            if (--value_ == 0) {
                value_ = size_;
                reset = true;
            }
            lock_.unlock();
            return reset;
        }
    };
} // namespace rna
#endif // RNAALIGNMENT_ATOMIC_GENERATION_H