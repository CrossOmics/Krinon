#include "../utilsRefactored/defines.h"
#include "GenomeIndex.h"
#include "../io/Parameters.h"

namespace RefactorProcessing {
    GenomeIndex::GenomeIndex() = default;

    GenomeIndex::~GenomeIndex() = default;

    void GenomeIndex::setParam(const Parameters &P) {
        kMerSize_ = P.kMerSize;
        kMerNum_ = (1ULL << (kMerSize_ * 2)) + 1;
        extendHashTableByte_ = P.extendAlternativeByte;
        extendHashTableNum_ = extendHashTableByte_/4;
        needInsertSJ_ = P.insertSJ;
        sjdbOverhang_ = P.sjdbOverhang;
        limitSjdbInsertN_ = P.limitSjdbInsertN;
        genome_.setParam(P);
    }

    void GenomeIndex::loadGenome() {
        genome_.loadFromFasta(genome_.genomeFileName_);
    }

    void GenomeIndex::build() {
        suffixArray_.build(genome_.sequence_);

        buildLCP();

        buildKmerMap();

        buildExtendedIndexHash();
    }


    void GenomeIndex::buildLCP() {
        longestCommonPrefix_.reserve(suffixArray_.length_ + 1);
        PackedArray rk;
        rk.initialize(suffixArray_.length_, suffixArray_.wordBits_);
        for (int64_t i = 0; i < suffixArray_.length_; ++i) {
            rk.setValue(i, 0);
        }
        for (int64_t i = 0; i < suffixArray_.length_; ++i) {
            rk.setValue(suffixArray_[i], i);
        }

        //todo optimize: multi-thread
        size_t k = 0;
        for (int64_t i = 0; i < suffixArray_.length_; ++i) {
            size_t j = rk[i];
            if (j == 0) {
                k = 0;
                continue;
            }
            if (k > 0) k--;
            while (genome_.sequence_[i + k] == genome_.sequence_[suffixArray_[j - 1] + k]) k++;
            if (k < std::numeric_limits<uint8_t>::max()) longestCommonPrefix_[j] = k;
            else longestCommonPrefix_[j] = std::numeric_limits<uint8_t>::max();
        }
    }

    inline void setFlag(int32_t hash, int32_t length, std::vector<std::vector<int8_t>> &appearance_flag) {
        if (length == 1) {
            appearance_flag[length][0] |= (1 << hash);
        } else {
            int32_t pos = hash >> 3;
            int32_t bit = hash & 0x7;
            appearance_flag[length][pos] |= (1 << bit);
        }
    }

    inline bool getFlag(int32_t hash,int32_t length, const std::vector<std::vector<int8_t>> &appearance_flag) {
        if (length == 1) {
            return (appearance_flag[length][0] >> hash) & 1;
        } else {
            int32_t pos = hash >> 3;
            int32_t bit = hash & 0x7;
            return (appearance_flag[length][pos] >> bit) & 1;
        }
    }

    void GenomeIndex::buildKmerMap() {
        // calculate the appearance of each k-mer in the genome
        std::vector<std::vector<int8_t>> appearance_flag; // records the appearance of each k-mer in the genome
        appearance_flag.resize(kMerSize_ + 1); // 0 is empty, for convenience
        for (int i = 1; i <= kMerSize_; ++i) {
            appearance_flag[i].resize(1 << (std::max(2 * i - 3, 0)), 0);
        }
        std::vector<int64_t> nowWindowHash;
        nowWindowHash.resize(kMerSize_ + 1, 0); // Also, 0 is empty
        int nowAvailableLength = 0;
        const std::string &seq = genome_.sequence_;
        int64_t genomeLength = genome_.sequence_.length();
        for (size_t i = 0; i < genomeLength; ++i) {
            charToIndex(seq[i]) >= 0 ? ++nowAvailableLength : nowAvailableLength = 0;
            if (nowAvailableLength == 0) continue;
            if (nowAvailableLength <= kMerSize_) {
                nowWindowHash[nowAvailableLength] = (nowWindowHash[nowAvailableLength - 1] << 2) | charToIndex(seq[i]);
                for (int j = 1; j < nowAvailableLength; ++j) {
                    nowWindowHash[j] = (nowWindowHash[j] << 2) | charToIndex(seq[i]);
                    nowWindowHash[j] &= (1 << (2 * j)) - 1; // keep only the last 2*j bits
                }
                for (int j = 1; j <= nowAvailableLength; ++j) {
                    setFlag(nowWindowHash[j], j,appearance_flag);
                }
            } else {
                for (int j = 1; j <= kMerSize_; ++j) {
                    nowWindowHash[j] = (nowWindowHash[j] << 2) | charToIndex(seq[i]);
                    nowWindowHash[j] &= (1 << (2 * j)) - 1; // keep only the last 2*j bits
                    setFlag(nowWindowHash[j], j,appearance_flag);
                }
            }
        }

        // initialize patternMerMap_
        // In fact, length cannot be 0. ATCG all appear in the reference genome.
        // Therefore, we use 0 to show that the k-mer does not appear in the genome.
        // (For the sake of data compression)
        for (size_t i = 0; i < kMerNum_ - 1; ++i) {
            int32_t hash = i;
            int32_t length = kMerSize_;
            while (!getFlag(hash, length,appearance_flag) && length > 0) {
                --length;
                hash >>= 2; // remove the last two bits
            }
            patternMerMap_.set(i, patternMerMap_.INDEX_LENGTH, length);
        }

        int64_t prevHash = -1;
        //initialize
        for (int64_t i = 0; i < kMerNum_ - 1; ++i) {
            patternMerMap_.set(i, patternMerMap_.INDEX_LEFT_SA_INDEX, patternMerMap_.EMPTY_SA_INDEX);
            patternMerMap_.set(i, patternMerMap_.INDEX_UPPER_RANGE, patternMerMap_.EMPTY_UPPER_RANGE);
        }

        for (int64_t i = 0; i < suffixArray_.length_; ++i){
            int64_t hash = encodeKmer(seq.substr(suffixArray_[i],kMerSize_),kMerSize_);

            if (hash > prevHash){
                if (prevHash != -1 && patternMerMap_.get(prevHash,patternMerMap_.INDEX_UPPER_RANGE) == patternMerMap_.EMPTY_UPPER_RANGE)
                    patternMerMap_.set(prevHash,patternMerMap_.INDEX_UPPER_RANGE,(i - 1 - patternMerMap_.get(prevHash, patternMerMap_.INDEX_LEFT_SA_INDEX)));
                prevHash = hash;
            }

            if (hash == -1) {
                if (prevHash != -1 && patternMerMap_.get(prevHash,patternMerMap_.INDEX_UPPER_RANGE) == patternMerMap_.EMPTY_UPPER_RANGE)
                    patternMerMap_.set(prevHash,patternMerMap_.INDEX_UPPER_RANGE,(i - 1 - patternMerMap_.get(prevHash, patternMerMap_.INDEX_LEFT_SA_INDEX)));
                continue;
            }
            if (patternMerMap_.get(hash, patternMerMap_.INDEX_LEFT_SA_INDEX) == patternMerMap_.EMPTY_SA_INDEX)
                patternMerMap_.set(hash, patternMerMap_.INDEX_LEFT_SA_INDEX, i);
        }
        if (prevHash != -1 && patternMerMap_.get(prevHash,patternMerMap_.INDEX_UPPER_RANGE) == patternMerMap_.EMPTY_UPPER_RANGE)
            patternMerMap_.set(prevHash,patternMerMap_.INDEX_UPPER_RANGE,(suffixArray_.length_ - 1 - patternMerMap_.get(prevHash, patternMerMap_.INDEX_LEFT_SA_INDEX)));

        int64_t prevRightBound = 0;
        for (int64_t i = 0; i < kMerNum_ - 1; ++i) {
            if (patternMerMap_.get(i, patternMerMap_.INDEX_LENGTH) < kMerSize_) {
                patternMerMap_.set(i, patternMerMap_.INDEX_LEFT_SA_INDEX, prevRightBound + 1);
            } else {
                prevRightBound =
                        patternMerMap_.get(i, patternMerMap_.INDEX_LEFT_SA_INDEX) + patternMerMap_.get(i, patternMerMap_.INDEX_UPPER_RANGE);
            }
        }
    }



    void GenomeIndex::buildExtendedIndexHash() {
        extendedIndexHash_.resize(kMerNum_ * extendHashTableByte_ / 4, 0);
        for (size_t i = 0; i < kMerNum_ - 1; ++i) {
            buildExtendedIndexHashSingle(i);
        }
    }

    void GenomeIndex::buildExtendedIndexHashSingle(int64_t h) {
        int64_t nowLeftSAIndex = patternMerMap_.get(h,patternMerMap_.INDEX_LEFT_SA_INDEX);
        int64_t num = patternMerMap_.get(h, patternMerMap_.INDEX_UPPER_RANGE);
        int indNum = extendHashTableByte_ / 4;

        if (num > 16) {
            for (int j = 0; j < indNum; ++j) {
                int length = 16;
                uint64_t index = suffixArray_[nowLeftSAIndex +
                                              j * num * 4 / extendHashTableByte_];
                uint64_t hash = 0;

                for (int k = 0; k < length; ++k) {
                    int c = charToIndex(genome_.sequence_[index + kMerSize_ + k]);
                    if (c == -1) {
                        hash <<= 2;
                        hash |= c;
                    }else {
                        // set all remaining bps to T
                        hash <<= (2 * (length - k));
                        hash |= (1ULL << (2 * (length - k))) - 1;
                        break;
                    }
                }

                extendedIndexHash_[h * indNum + j] = (uint32_t)hash;
            }
        }
    }


}
