#include "PackedArray.h"
#include <fstream>

namespace RefactorProcessing {
    PackedArray::PackedArray() : data(nullptr), allocated(false) {};

    PackedArray::~PackedArray() {
        if (allocated) delete[] data;
    }

    void PackedArray::initialize(int64_t wordNum, int wordLengthBits, int reservedLength) {
        wordNum_ = wordNum;
        wordLengthBits_ = wordLengthBits;
        reservedLength_ = reservedLength;
        wordLengthBits_ = wordLengthBits;
        wordCompLength_ = 64 - wordLengthBits_;
        bitMask_ = (~0ULL) >> wordCompLength_;
        arrayLength_ = ((wordNum_ + reservedLength_) * wordLengthBits_ + 63) / 64;
        data = new uint64_t[arrayLength_];
    }

    void PackedArray::setValue(int64_t index, uint64_t value) {
        index += reservedLength_;
        uint64_t b = index * wordLengthBits_;
        uint64_t B = b / 64;
        uint64_t S = b % 64;
        auto wordPtr = data + B;
        *wordPtr = (*wordPtr & ~(bitMask_ << S)) | (value << S);
        if (wordLengthBits_ + S > 64) {
            // Handle the case where the value spans three 32-bit integers
            *(wordPtr + 1) = (*(wordPtr + 1) & ~(bitMask_ >> (64 - S))) | (value >> (64 - S));
        }
    }

    uint64_t PackedArray::getValue(int64_t index) const {
        index += reservedLength_;
        uint64_t b = index * wordLengthBits_;
        uint64_t B = b / 64;
        uint64_t S = b % 64;
        auto wordPtr = data + B;
        uint64_t value = (*wordPtr >> S) & bitMask_;
        if (wordLengthBits_ + S > 64) {
            // Handle the case where the value spans three 32-bit integers
            value |= (*(wordPtr + 1) & (bitMask_ >> (64 - S))) << (64 - S);
        }
        return value;
    }

    uint64_t PackedArray::operator[](int64_t index) const {
        return getValue(index);
    }

    bool PackedArray::setReservedLength(int64_t newReservedLength) {
        if (newReservedLength > reservedLength_) return false;
        int64_t oldReservedLength = reservedLength_;
        reservedLength_ = newReservedLength;
        int64_t shiftLength = oldReservedLength - newReservedLength;
        for (int64_t i = 0; i < wordNum_; ++i) {
            int64_t v = getValue(i + shiftLength);
            setValue(i, v);
        }
        return true;
    }

    int PackedArray::writeToFile(const std::string &fileName) const {
        std::ofstream out(fileName, std::ios::binary);
        if (!out.is_open()) return -1;
        out.write(reinterpret_cast<const char*>(&wordLengthBits_), sizeof(wordLengthBits_));
        out.write(reinterpret_cast<const char*>(&wordNum_), sizeof(wordNum_));
        out.write(reinterpret_cast<const char*>(&arrayLength_), sizeof(arrayLength_));
        out.write(reinterpret_cast<const char*>(&reservedLength_), sizeof(reservedLength_));
        out.write(reinterpret_cast<const char*>(data), arrayLength_ * sizeof(uint64_t));
        out.close();
        return 0;
    }

    int PackedArray::loadFromFile(const std::string &fileName, int reservedLength) {
        //todo optimize loading speed
        std::ifstream in(fileName, std::ios::binary);
        if (!in.is_open()) return -1;
        in.read(reinterpret_cast<char*>(&wordLengthBits_), sizeof(wordLengthBits_));
        in.read(reinterpret_cast<char*>(&wordNum_), sizeof(wordNum_));
        in.read(reinterpret_cast<char*>(&arrayLength_), sizeof(arrayLength_));
        in.read(reinterpret_cast<char*>(&reservedLength_), sizeof(reservedLength_));
        wordCompLength_ = 64 - wordLengthBits_;
        bitMask_ = (~0ULL) >> wordCompLength_;
        if (allocated) delete[] data;
        data = new uint64_t[arrayLength_];
        in.read(reinterpret_cast<char*>(data), arrayLength_ * sizeof(uint64_t));
        in.close();
        reservedLength_ = reservedLength;
        return 0;
    }

    void PackedArray::buildFromReservedArray(uint64_t *array, uint32_t wordLengthBits, uint64_t length,uint64_t reservedLength) {
        reservedLength_ = reservedLength;
        arrayLength_ = length + reservedLength;
        wordNum_ = length;
        wordLengthBits_ = wordLengthBits;
        wordCompLength_ = 64 - wordLengthBits_;
        bitMask_ = (~0ULL) >> wordCompLength_;
        data = array;
        allocated = false;
    }

    void PackedArray::setLength(int64_t length) {
        wordNum_ = length;
        arrayLength_ = ((wordNum_ + reservedLength_) * wordLengthBits_ + 63) / 64;
    }

}