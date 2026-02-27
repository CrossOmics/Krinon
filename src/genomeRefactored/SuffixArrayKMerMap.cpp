#include <cstdlib>
#include "SuffixArrayKMerMap.h"
#include "../log/ErrorRecord.h"

namespace RefactorProcessing{
    SuffixArrayKMerMap::SuffixArrayKMerMap(int INDEX_LENGTH_BITS, int INDEX_LEFT_SA_INDEX_BITS, int INDEX_UPPER_RANGE_BITS) {
        ELEMENT_BIT_LENGTHS[INDEX_LENGTH] = INDEX_LENGTH_BITS;
        ELEMENT_BIT_LENGTHS[INDEX_LEFT_SA_INDEX] = INDEX_LEFT_SA_INDEX_BITS;
        ELEMENT_BIT_LENGTHS[INDEX_UPPER_RANGE] = INDEX_UPPER_RANGE_BITS;
        EMPTY_UPPER_RANGE = (1ULL << INDEX_UPPER_RANGE_BITS) - 1;
        EMPTY_SA_INDEX = (1ULL << INDEX_LEFT_SA_INDEX_BITS) - 1;
        BIT_OFFSETS[INDEX_LENGTH] = 0;
        BIT_OFFSETS[INDEX_LEFT_SA_INDEX] = INDEX_LENGTH_BITS;
        BIT_OFFSETS[INDEX_UPPER_RANGE] = INDEX_LENGTH_BITS + INDEX_LEFT_SA_INDEX_BITS;
        dataLengthLessThanLL_ = (ELEMENT_BIT_LENGTHS[INDEX_LEFT_SA_INDEX] + ELEMENT_BIT_LENGTHS[INDEX_UPPER_RANGE]) <= 64;
        data_ = nullptr;
    }

    SuffixArrayKMerMap::~SuffixArrayKMerMap() {
        if (data_) {
            delete[] data_;
            data_ = nullptr;
        }
    }

    void SuffixArrayKMerMap::init(uint64_t length) {
        size_t totalBits = ELEMENT_BIT_LENGTHS[INDEX_LENGTH] + ELEMENT_BIT_LENGTHS[INDEX_LEFT_SA_INDEX] + ELEMENT_BIT_LENGTHS[INDEX_UPPER_RANGE];
        size_t bitsPerElement = dataLengthLessThanLL_ ? 64 : totalBits;
        size_t totalLengthBits = length * bitsPerElement;
        size_t dataLength = (totalLengthBits + 63) / 64;
        data_ = new uint64_t[dataLength]();
        length_ = length;
    }

    int64_t SuffixArrayKMerMap::size() const {
        return length_;
    }

    void SuffixArrayKMerMap::set(uint64_t index, int elementIndex, int64_t value) {
        if (value < 0) {
            rna::ErrorRecord().reportError("value cannot be negative in compressed KMerMap");
        }
        if (dataLengthLessThanLL_) {
            auto wordPtr = data_ + index;
            *wordPtr = (*wordPtr & ~( ((1ULL << ELEMENT_BIT_LENGTHS[elementIndex]) - 1) << BIT_OFFSETS[elementIndex])) | ((uint64_t)value << BIT_OFFSETS[elementIndex]) ;
        } else {
            uint64_t b = index * (ELEMENT_BIT_LENGTHS[INDEX_LENGTH] + ELEMENT_BIT_LENGTHS[INDEX_LEFT_SA_INDEX] + ELEMENT_BIT_LENGTHS[INDEX_UPPER_RANGE]) + BIT_OFFSETS[elementIndex];
            uint64_t B = b / 64;
            uint64_t S = b % 64;
            auto wordPtr = data_ + B;
            *wordPtr = (*wordPtr & ~(((1ULL << ELEMENT_BIT_LENGTHS[elementIndex]) - 1) << S)) | ((uint64_t)value << S);
            if (ELEMENT_BIT_LENGTHS[elementIndex] + S > 64) {
                // Handle the case where the value spans two 64-bit integers
                *(wordPtr + 1) = (*(wordPtr + 1) & ~(((1ULL << ELEMENT_BIT_LENGTHS[elementIndex]) - 1) >> (64 - S))) | ((uint64_t)value >> (64 - S));
            }
        }
    }

    int64_t SuffixArrayKMerMap::get(uint64_t index, int elementIndex) const {
        if (dataLengthLessThanLL_) {
            auto wordPtr = data_ + index;
            uint64_t value = (*wordPtr >> BIT_OFFSETS[elementIndex]) & ((1ULL << ELEMENT_BIT_LENGTHS[elementIndex]) - 1);

            return (int64_t)value;
        } else {
            uint64_t b = index * (ELEMENT_BIT_LENGTHS[INDEX_LENGTH] + ELEMENT_BIT_LENGTHS[INDEX_LEFT_SA_INDEX] + ELEMENT_BIT_LENGTHS[INDEX_UPPER_RANGE]) + BIT_OFFSETS[elementIndex];
            uint64_t B = b / 64;
            uint64_t S = b % 64;
            auto wordPtr = data_ + B;
            uint64_t value = (*wordPtr >> S) & ((1ULL << ELEMENT_BIT_LENGTHS[elementIndex]) - 1);
            if (ELEMENT_BIT_LENGTHS[elementIndex] + S > 64) {
                // Handle the case where the value spans two 64-bit integers
                value |= (*(wordPtr + 1) & (((1ULL << ELEMENT_BIT_LENGTHS[elementIndex]) - 1) >> (64 - S))) << (64 - S);
            }

            return (int64_t)value;
        }
    }

    int SuffixArrayKMerMap::writeToFile(const std::string &fileName) const {
        std::ofstream ofs(fileName, std::ios::binary);
        if (!ofs.is_open()) return -1;

        // magic + version
        const char magic[4] = {'K','M','R','P'}; // KMerMap marker
        ofs.write(magic, 4);


        // write numeric members
        ofs.write(reinterpret_cast<const char*>(&ELEMENT_BIT_LENGTHS[INDEX_LENGTH]), sizeof(ELEMENT_BIT_LENGTHS[INDEX_LENGTH]));
        ofs.write(reinterpret_cast<const char*>(&ELEMENT_BIT_LENGTHS[INDEX_LEFT_SA_INDEX]), sizeof(ELEMENT_BIT_LENGTHS[INDEX_LEFT_SA_INDEX]));
        ofs.write(reinterpret_cast<const char*>(&ELEMENT_BIT_LENGTHS[INDEX_UPPER_RANGE]), sizeof(ELEMENT_BIT_LENGTHS[INDEX_UPPER_RANGE]));
        ofs.write(reinterpret_cast<const char*>(&length_), sizeof(length_));


        //write data
        size_t totalBits = ELEMENT_BIT_LENGTHS[INDEX_LENGTH] + ELEMENT_BIT_LENGTHS[INDEX_LEFT_SA_INDEX] + ELEMENT_BIT_LENGTHS[INDEX_UPPER_RANGE];
        size_t bitsPerElement = dataLengthLessThanLL_ ? 64 : totalBits;
        size_t totalLengthBits = length_ * bitsPerElement;
        size_t dataLength = (totalLengthBits + 63) / 64;
        ofs.write(reinterpret_cast<const char*>(data_), dataLength * sizeof(uint64_t));

        return 0;
    }

    int SuffixArrayKMerMap::loadFromFile(const std::string &fileName) {
        std::ifstream ifs(fileName, std::ios::binary);
        if (!ifs.is_open()) return -1;

        char magic[4];
        ifs.read(magic, 4);
        if (!ifs) return -2;
        if (!(magic[0]=='K' && magic[1]=='M' && magic[2]=='R' && magic[3]=='P')) return -3;


        int indexLengthBits, indexLeftSAIndexBits, indexUpperRangeBits;
        ifs.read(reinterpret_cast<char*>(&indexLengthBits), sizeof(indexLengthBits));
        ifs.read(reinterpret_cast<char*>(&indexLeftSAIndexBits), sizeof(indexLeftSAIndexBits));
        ifs.read(reinterpret_cast<char*>(&indexUpperRangeBits), sizeof(indexUpperRangeBits));
        if (!ifs) return -5;

        length_ = 0;
        ifs.read(reinterpret_cast<char*>(&length_), sizeof(length_));
        if (!ifs) return -6;

        // reinitialize the object with the loaded parameters
        ELEMENT_BIT_LENGTHS[INDEX_LENGTH] = indexLengthBits;
        ELEMENT_BIT_LENGTHS[INDEX_LEFT_SA_INDEX] = indexLeftSAIndexBits;
        ELEMENT_BIT_LENGTHS[INDEX_UPPER_RANGE] = indexUpperRangeBits;
        EMPTY_UPPER_RANGE = (1ULL << indexUpperRangeBits) - 1;
        EMPTY_SA_INDEX = (1ULL << indexLeftSAIndexBits) - 1;
        BIT_OFFSETS[INDEX_LENGTH] = 0;
        BIT_OFFSETS[INDEX_LEFT_SA_INDEX] = indexLengthBits;
        BIT_OFFSETS[INDEX_UPPER_RANGE] = indexLengthBits + indexLeftSAIndexBits;
        dataLengthLessThanLL_ = (ELEMENT_BIT_LENGTHS[INDEX_LEFT_SA_INDEX] + ELEMENT_BIT_LENGTHS[INDEX_UPPER_RANGE]) <= 64;
        init(length_);




        // read data
        size_t totalBits = ELEMENT_BIT_LENGTHS[INDEX_LENGTH] + ELEMENT_BIT_LENGTHS[INDEX_LEFT_SA_INDEX] + ELEMENT_BIT_LENGTHS[INDEX_UPPER_RANGE];
        size_t bitsPerElement = dataLengthLessThanLL_ ? 64 : totalBits;
        size_t totalLengthBits = length_ * bitsPerElement;
        size_t dataLength = (totalLengthBits + 63) / 64;
        ifs.read(reinterpret_cast<char*>(data_), dataLength * sizeof(uint64_t));

        return 0;

    }


}
