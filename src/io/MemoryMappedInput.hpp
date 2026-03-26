#ifndef RNAALIGNREFACTORED_MEMORYMAPPEDINPUT_HPP
#define RNAALIGNREFACTORED_MEMORYMAPPEDINPUT_HPP
#include <atomic>
#include <string>
#include <thread>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <stdexcept>
#include <filesystem>

namespace RefactorProcessing {
    /**
     * Thread-safe wrapper for a memory-mapped input file. It assumes that the input
     * is read sequantially and by multiple threads.
     */
    class MemoryMappedInput {
    private:
        size_t size_;
        int fd_;
        char* mapPtr_;
        char* endPtr_;
        std::atomic<size_t> offset_;

    public:
        MemoryMappedInput(const std::string& path) :
                fd_(-1),
                mapPtr_(nullptr),
                offset_{0}
        {
            size_ = std::filesystem::file_size(path.c_str());
            
            fd_ = open(path.c_str(), O_RDONLY, S_IRUSR | S_IWUSR);
            if (fd_ == -1) {
                throw std::runtime_error("Failed get file handle for " + path);
            }

            mapPtr_ = static_cast<char*>(mmap(NULL, size_, PROT_READ, MAP_SHARED, fd_, 0));
            if (mapPtr_ == MAP_FAILED) {
                close(fd_);
                throw std::runtime_error("Failed to map " + path + " into memory");
            }
            
            madvise(mapPtr_, size_, MADV_SEQUENTIAL);
            endPtr_ = mapPtr_ + size_;
        }

        ~MemoryMappedInput() {
            memClose();
        }

        MemoryMappedInput(const MemoryMappedInput&) = delete;
        MemoryMappedInput& operator=(const MemoryMappedInput&) = delete;

        const char* getMapPtr() {
            return mapPtr_;
        }
        const char* getEndPtr() {
            return endPtr_;
        }
        size_t size() const {
            return size_;
        }

        /**
         * Atomically claim the interval `[mapPtr_ + offset_, mapPtr_ + offset_ + howMuch)`.
         * Returns pointer to the start of the interval.
         */
        const char* claim(size_t howMuch) {
            return mapPtr_ + offset_.fetch_add(howMuch, std::memory_order_relaxed);
        }

        void memClose() {
            if (mapPtr_ != nullptr && mapPtr_ != MAP_FAILED)
                munmap(mapPtr_, size_);
            if (fd_ != -1)
                close(fd_);
            mapPtr_ = nullptr;
            fd_ = -1;
        }
    };
}
#endif //RNAALIGNREFACTORED_MEMORYMAPPEDINPUT_HPP