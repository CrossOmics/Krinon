#include "PackedArray.h"
namespace rna {

    void PackedArray::set(uint64_t index, uint64_t value) {
            uint64_t b = index * wordLengthBits_;
            uint64_t B = b / 64;
            uint64_t S = b % 64;
            auto wordPtr = data_ + B;
            *wordPtr = (*wordPtr & ~(bitMask_ << S)) | (value << S);
            if (wordLengthBits_ + S > 64) {
                // Handle the case where the value spans three 32-bit integers
                *(wordPtr + 1) = (*(wordPtr + 1) & ~(bitMask_ >> (64 - S))) | (value >> (64 - S));
            }
    }
    void PackedArray::init(uint64_t length, uint32_t wordLengthBits) {
        length_ = length;
        wordLengthBits_ = wordLengthBits;
        dataLength = (length * wordLengthBits + 63) / 64; // calculate the number of 32-bit integers needed
        data_ = new uint64_t[dataLength]();
        wordCompLength_ = 64 - wordLengthBits_;
        bitMask_ = (~0ULL) >> wordCompLength_;
    }
    uint64_t PackedArray::get(uint64_t index) const {
        uint64_t b = index * wordLengthBits_;
        uint64_t B = b / 64;
        uint64_t S = b % 64;
        auto wordPtr = data_ + B;
        uint64_t value = (*wordPtr >> S) & bitMask_;
        if (wordLengthBits_ + S > 64) {
            // Handle the case where the value spans three 32-bit integers
            value |= (*(wordPtr + 1) & (bitMask_ >> (64 - S))) << (64 - S);
        }
        return value;
    }
    void PackedArray::clear() {
        if (data_) {
            delete[] data_;
            data_ = nullptr;
        }
        wordLengthBits_ = 0;
        wordCompLength_ = 0;
        bitMask_ = 0;
        length_ = 0;
        dataLength = 0;
    }
    void PackedArray::writeToFile(std::ofstream &out, std::ofstream &logOut) const {
        logOut<< wordLengthBits_ << " " << length_ << " " << dataLength << "\n";
        out.write(reinterpret_cast<const char*>(data_), dataLength * sizeof(uint64_t));
    }
    void PackedArray::loadFromFile(std::ifstream &in, std::ifstream &logIn) {
        logIn >> wordLengthBits_ >> length_ >> dataLength;
        wordCompLength_ = 64 - wordLengthBits_;
        bitMask_ = (~0ULL) >> wordCompLength_;
        delete[] data_;
        data_ = new uint64_t[dataLength];
        in.read(reinterpret_cast<char*>(data_), dataLength * sizeof(uint64_t));
    }
}
