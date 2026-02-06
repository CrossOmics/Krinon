#ifndef RNAALIGNREFACTORED_MEMORYMAPPEDFILER_HPP
#define RNAALIGNREFACTORED_MEMORYMAPPEDFILER_HPP
#include <string>
#include <thread>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <stdexcept>
namespace RefactorProcessing {
    class MemoryMappedFile {
    private:
        int64_t size_;
        int fd_;
        char* mapPtr_;

    public:
        MemoryMappedFile(const std::string& path, long long maxSize) :
                size_(maxSize),
                fd_(-1),
                mapPtr_(nullptr)
        {
            fd_ = open(path.c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
            if (fd_ == -1) {
                throw std::runtime_error("Failed get file handle for " + path);
            }

            if (ftruncate(fd_, size_) == -1) {
                close(fd_);
                throw std::runtime_error("Failed to truncate " + path);
            }

            mapPtr_ = static_cast<char*>(mmap(NULL, size_, PROT_WRITE, MAP_SHARED, fd_, 0));
            if (mapPtr_ == MAP_FAILED) {
                close(fd_);
                throw std::runtime_error("Failed to map " + path + " into memory");
            }
        }

        ~MemoryMappedFile() {
            memClose();
        }

        MemoryMappedFile(const MemoryMappedFile&) = delete;
        MemoryMappedFile& operator=(const MemoryMappedFile&) = delete;

        char* getMapPtr() {
            return mapPtr_;
        }

        int64_t size() const {
            return size_;
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
#endif //RNAALIGNREFACTORED_MEMORYMAPPEDFILER_HPP
