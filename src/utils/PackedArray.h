#ifndef RNAALIGNREFACTORED_PACKEDARRAY_H
#define RNAALIGNREFACTORED_PACKEDARRAY_H

#include <cstdint>
#include "types.h"
#include <fstream>

namespace rna {
    class PackedArray {
    private:
        uint64_t *data_;
        uint32_t wordLengthBits_;
        uint32_t wordCompLength_;
        uint64_t bitMask_;
        uint64_t length_;//number of elements
        uint64_t dataLength;//length of the array (i.e. data)
    public:
        PackedArray() : data_(nullptr), wordLengthBits_(0), length_(0), dataLength(0) {};

        ~PackedArray() {
            delete[] data_;
        }

        PackedArray(uint64_t length, uint32_t wordLengthBits)
                : length_(length), wordLengthBits_(wordLengthBits) {
            dataLength = (length * wordLengthBits + 63) / 64; // calculate the number of 32-bit integers needed
            data_ = new uint64_t[dataLength]();
            wordCompLength_ = 64 - wordLengthBits_;
            bitMask_ = (~0ULL) >> wordCompLength_;
        }

        void init(uint64_t length, uint32_t wordLengthBits);

        inline void setLength(uint64_t length) {
            length_ = length;
            dataLength = (length * wordLengthBits_ + 63) / 64; // recalculate the number of 32-bit integers needed
        }

        void buildFromArray(uint64_t *array, uint32_t wordLengthBits, uint64_t length) {
            setLength(length);
            wordLengthBits_ = wordLengthBits;
            wordCompLength_ = 64 - wordLengthBits_;
            bitMask_ = (~0ULL) >> wordCompLength_;
            data_ = array;
        }

        void set(uint64_t index, uint64_t value);

        inline uint64_t length() const {
            return length_;
        }

        uint64_t get(uint64_t index) const;

        void clear();

        void writeToFile(std::ofstream &out, std::ofstream &logOut) const;

        void loadFromFile(std::ifstream &in, std::ifstream &logIn);


    };
}

#endif //RNAALIGNREFACTORED_PACKEDARRAY_H
