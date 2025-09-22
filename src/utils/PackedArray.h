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
        uint64_t dataLength;//length of the array to store the elements (i.e. data)
        uint64_t reservedInd;
    public:
        bool isOwner{true}; // is the Owner of data_
        PackedArray() : data_(nullptr), wordLengthBits_(0), length_(0), dataLength(0), reservedInd(0) {};

        ~PackedArray() {
            if(isOwner) delete[] data_;
        }


        void init(uint64_t length, uint32_t wordLengthBits,size_t reservedLength = 0);

        inline int getWordLengthBits() const {
            return wordLengthBits_;
        }

        inline void setLength(uint64_t length) {
            length_ = length;
            dataLength = ((length+reservedInd) * wordLengthBits_ + 63) / 64; // recalculate the number of 32-bit integers needed
        }

        inline void setReservedInd(uint64_t reservedInd){
            this->reservedInd = reservedInd;
            dataLength = ((length_+reservedInd) * wordLengthBits_ + 63) / 64; // recalculate the number of 32-bit integers needed
        }

        void buildFromArray(uint64_t *array, uint32_t wordLengthBits, uint64_t length) {
            setLength(length);
            wordLengthBits_ = wordLengthBits;
            wordCompLength_ = 64 - wordLengthBits_;
            bitMask_ = (~0ULL) >> wordCompLength_;
            data_ = array;
        }

        void buildFromPackedArrayExtendForward(PackedArray &other, uint64_t forwardExtendLength);

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
