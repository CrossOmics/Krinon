#ifndef RNAALIGNREFACTORED_PACKEDINDEXARRAY_H
#define RNAALIGNREFACTORED_PACKEDINDEXARRAY_H
#include <array>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <fstream>

namespace rna {
    // packed array for storing fixed-length index structures
    // no need to reserve memory
    class PackedIndexArray {
    private:
        uint64_t *data_;
        std::vector<int> elementBitLengths_;
        std::vector<uint64_t> bitMasks_;
        std::vector<int> bitOffsets_;
        int structLengthBits_;
        uint64_t length_;
        uint64_t dataLength_;
        bool lessThanLL;
    public:
        PackedIndexArray() = default;
        PackedIndexArray(std::vector<int> elementBitLengths)
                : data_(nullptr), elementBitLengths_(std::move(elementBitLengths)), structLengthBits_(0), length_(0),
                  dataLength_(0),lessThanLL(false) {

            for (int bitLength: elementBitLengths_) {
                bitOffsets_.push_back(structLengthBits_);
                structLengthBits_ += bitLength;
                bitMasks_.push_back((~0ULL) >> (64 - bitLength));
            }
            if (structLengthBits_ <= 64) {
                lessThanLL = true;
                structLengthBits_ = 64;
            }

        }

        ~PackedIndexArray() {
            delete[] data_;
        }

        int64_t size() const {
            return length_;
        }

        void init(uint64_t length){
            length_ = length;
            dataLength_ = (length * structLengthBits_ + 63) / 64; // calculate the number of 64-bit integers needed
            data_ = new uint64_t[dataLength_]();
        }

        void set(uint64_t index, int elementIndex, int64_t value){
            // handle minus value
            if (value < 0){
                uint64_t maxValue = (1ULL << elementBitLengths_[elementIndex]);
                value = maxValue + value;
            }
            if (lessThanLL){
                auto wordPtr = data_ + index;
                *wordPtr = (*wordPtr & ~(bitMasks_[elementIndex] << bitOffsets_[elementIndex])) | ((uint64_t)value << bitOffsets_[elementIndex]);
            }else{
                uint64_t b = index * structLengthBits_ + bitOffsets_[elementIndex];
                uint64_t B = b / 64;
                uint64_t S = b % 64;
                auto wordPtr = data_ + B;
                *wordPtr = (*wordPtr & ~(bitMasks_[elementIndex] << S)) | ((uint64_t)value << S);
                if (bitOffsets_[elementIndex] + elementBitLengths_[elementIndex] + S > 64) {
                    // Handle the case where the value spans two 64-bit integers
                    *(wordPtr + 1) = (*(wordPtr + 1) & ~(bitMasks_[elementIndex] >> (64 - S))) | ((uint64_t)value >> (64 - S));
                }
            }

        }

        int64_t get(uint64_t index, int elementIndex) const{
            if(lessThanLL){
                auto wordPtr = data_ + index;
                uint64_t value = (*wordPtr >> bitOffsets_[elementIndex]) & bitMasks_[elementIndex];
                if ((elementIndex != 0) && (value == (1ULL << elementBitLengths_[elementIndex]) - 1)) {
                    // -1
                    value = -1;
                }
                return (int64_t)value;
            }else{
                uint64_t b = index * structLengthBits_ + bitOffsets_[elementIndex];
                uint64_t B = b / 64;
                uint64_t S = b % 64;
                auto wordPtr = data_ + B;
                uint64_t value = (*wordPtr >> S) & bitMasks_[elementIndex];
                if (bitOffsets_[elementIndex] + elementBitLengths_[elementIndex] + S > 64) {
                    // Handle the case where the value spans two 64-bit integers
                    value |= (*(wordPtr + 1) & (bitMasks_[elementIndex] >> (64 - S))) << (64 - S);
                }

                if ((elementIndex == 1) && ((value >> (elementBitLengths_[elementIndex] - 1)) & 1)) {
                    // negative number
                    uint64_t maxValue = (1ULL << elementBitLengths_[elementIndex]);
                    value = value - maxValue;
                }

                return (int64_t)value;
            }

        }

        void writeToFile(std::ostream &out,std::ofstream &logOut) const {
            logOut << "IndexData:\n";
            logOut << structLengthBits_ << " " << length_ << " " << dataLength_ << " " << lessThanLL <<" " <<elementBitLengths_.size() <<"\n";
            for (int bitLength: elementBitLengths_) {
                logOut << bitLength << " ";
            }
            logOut << "\n";

            for (int bitOffset: bitOffsets_) {
                logOut << bitOffset << " ";
            }
            logOut << "\n";
            for (uint64_t bitMask: bitMasks_) {
                logOut << bitMask << " ";
            }
            logOut << "\n";

            out.write(reinterpret_cast<const char *>(data_), dataLength_ * sizeof(uint64_t));
        }
        void loadFromFile(std::istream &in, std::ifstream &logIn) {
            std::string n;
            int elementNum;
            logIn >> n; // Index Data:
            logIn >> structLengthBits_ >> length_ >> dataLength_ >> lessThanLL >> elementNum;

            elementBitLengths_.resize(elementNum);
            for (int i = 0; i < elementNum; ++i) {
                logIn >> elementBitLengths_[i];
            }
            bitOffsets_.resize(elementNum);
            for (int i = 0; i < elementNum; ++i) {
                logIn >> bitOffsets_[i];
            }
            bitMasks_.resize(elementNum);
            for (int i = 0; i < elementNum; ++i) {
                logIn >> bitMasks_[i];
            }
            //delete[] data_;
            data_ = new uint64_t[dataLength_];

            in.read(reinterpret_cast<char *>(data_), dataLength_ * sizeof(uint64_t));
        }
    };
}
#endif //RNAALIGNREFACTORED_PACKEDINDEXARRAY_H
