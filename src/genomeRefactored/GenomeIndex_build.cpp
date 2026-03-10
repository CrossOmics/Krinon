#include "../utilsRefactored/defines.h"
#include "GenomeIndex.h"
#include "../utilsRefactored/seqFunctions.h"
#include "../io/Parameters.h"
#include "../utilsRefactored/timeFunctions.h"

namespace RefactorProcessing {
    GenomeIndex::GenomeIndex(){
        //std::cout << "GenomeIndex constructor\n";
    };

    GenomeIndex::~GenomeIndex(){
        std::cout << "GenomeIndex destructor\n";
    };

    void GenomeIndex::setParam(const Parameters &P) {
        suffixArray_.setParam(P);
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
        std::cout << getCurrentTimeString() << "\n";
        std::cout << "Loading genome...\n" ;
        loadGenome();
        std::cout << getCurrentTimeString() << "\n";
        std::cout << "Building suffix array...\n" ;
        suffixArray_.build(genome_.sequence_);
        std::cout << getCurrentTimeString() << "\n";
        std::cout << "Building LCP array...\n" ;
        buildLCP();
        std::cout << getCurrentTimeString() << "\n";
        std::cout << "Building SAI...\n" ;
        buildKmerMap();
        std::cout << getCurrentTimeString() << "\n";
        std::cout << "Building extended SAI...\n" ;
        buildExtendedIndexHash();
        std::cout << getCurrentTimeString() << "\n";
        std::cout << "Finished building genome index\n" ;
    }


    void GenomeIndex::buildLCP() {

        longestCommonPrefix_.resize(suffixArray_.length_ + 1);
        PackedArray rk;
        size_t sequenceLength = genome_.sequence_.length();
        rk.initialize(sequenceLength, suffixArray_.wordBits_);
        for (size_t i = 0; i < sequenceLength; ++i) {
            rk.setValue(i, 0);
        }
        for (size_t i = 0; i < suffixArray_.length_; ++i) {
            rk.setValue(suffixArray_[i], i);
        }

        //todo optimize: multi-thread
        size_t k = 0;
        for (size_t i = 0; i < sequenceLength; ++i) {
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
        size_t genomeLength = genome_.sequence_.length();
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
        patternMerMap_.init(kMerNum_);
        for (size_t i = 0; i < kMerNum_ - 1; ++i) {
            int32_t hash = i;
            int32_t length = kMerSize_;
            while (length > 0 && !getFlag(hash, length,appearance_flag)) {
                --length;
                hash >>= 2; // remove the last two bits
            }
            patternMerMap_.set(i, patternMerMap_.INDEX_LENGTH, length);
        }

        int64_t prevHash = -1;
        //initialize
        for (size_t i = 0; i < kMerNum_ - 1; ++i) {
            patternMerMap_.set(i, patternMerMap_.INDEX_LEFT_SA_INDEX, patternMerMap_.EMPTY_SA_INDEX);
            patternMerMap_.set(i, patternMerMap_.INDEX_UPPER_RANGE, patternMerMap_.EMPTY_UPPER_RANGE);
        }

        for (size_t i = 0; i < suffixArray_.length_; ++i){
            int64_t hash = encodeKmer(seq.substr(suffixArray_[i],kMerSize_),kMerSize_);
            if (hash < 0) hash = -1;

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
        for (size_t i = 0; i < kMerNum_ - 1; ++i) {
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
                    if (c != -1) {
                        hash <<= 2;
                        hash |= c;
                    }else {
                        // non-ACGT
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



    int GenomeIndex::writeToFile(const std::string &fileName) const {
        std::ofstream ofs(fileName, std::ios::binary);
        if (!ofs.is_open()) return -1;

        // header: magic + version
        const char magic[4] = {'G','I','D','X'};
        ofs.write(magic, 4);
        uint32_t version = 1;
        ofs.write(reinterpret_cast<const char*>(&version), sizeof(version));

        // write simple members
        ofs.write(reinterpret_cast<const char*>(&kMerSize_), sizeof(kMerSize_));
        ofs.write(reinterpret_cast<const char*>(&kMerNum_), sizeof(kMerNum_));

        writeString(ofs, additionalIndexType_);

        ofs.write(reinterpret_cast<const char*>(&extendHashTableByte_), sizeof(extendHashTableByte_));
        ofs.write(reinterpret_cast<const char*>(&extendHashTableNum_), sizeof(extendHashTableNum_));

        ofs.write(reinterpret_cast<const char*>(&maxAlignNum), sizeof(maxAlignNum));

        ofs.write(reinterpret_cast<const char*>(&needInsertSJ_), sizeof(needInsertSJ_));
        ofs.write(reinterpret_cast<const char*>(&sjdbOverhang_), sizeof(sjdbOverhang_));
        ofs.write(reinterpret_cast<const char*>(&limitSjdbInsertN_), sizeof(limitSjdbInsertN_));

        // write vectors
        uint64_t lcpSize = static_cast<uint64_t>(longestCommonPrefix_.size());
        ofs.write(reinterpret_cast<const char*>(&lcpSize), sizeof(lcpSize));
        if (lcpSize) ofs.write(reinterpret_cast<const char*>(longestCommonPrefix_.data()), lcpSize * sizeof(uint8_t));

        uint64_t extSize = static_cast<uint64_t>(extendedIndexHash_.size());
        ofs.write(reinterpret_cast<const char*>(&extSize), sizeof(extSize));
        if (extSize) ofs.write(reinterpret_cast<const char*>(extendedIndexHash_.data()), extSize * sizeof(uint32_t));

        if (!ofs) return -2;

        // write complex objects into separate files for speed and modularity
        {
            std::string gfile = fileName + ".gen";
            if (genome_.writeToFile(gfile) != 0) return -3;
        }
        {
            std::string sfile = fileName + ".sa";
            if (suffixArray_.writeToFile(sfile) != 0) return -4;
        }
        {
            std::string kfile = fileName + ".saI";
            if (patternMerMap_.writeToFile(kfile) != 0) return -5;
        }

        return 0;
    }

    int GenomeIndex::loadFromFile(const std::string &fileName) {
        std::ifstream ifs(fileName, std::ios::binary);
        if (!ifs.is_open()) return -1;

        // header validation
        char magic[4];
        ifs.read(magic, 4);
        if (!ifs) return -2;
        if (!(magic[0]=='G' && magic[1]=='I' && magic[2]=='D' && magic[3]=='X')) return -3;

        uint32_t version;
        ifs.read(reinterpret_cast<char*>(&version), sizeof(version));
        if (!ifs) return -4;
        // if needed, handle different versions

        // read simple members
        ifs.read(reinterpret_cast<char*>(&kMerSize_), sizeof(kMerSize_));
        ifs.read(reinterpret_cast<char*>(&kMerNum_), sizeof(kMerNum_));

        if (!readString(ifs, additionalIndexType_)) return -5;

        ifs.read(reinterpret_cast<char*>(&extendHashTableByte_), sizeof(extendHashTableByte_));
        ifs.read(reinterpret_cast<char*>(&extendHashTableNum_), sizeof(extendHashTableNum_));

        ifs.read(reinterpret_cast<char*>(&maxAlignNum), sizeof(maxAlignNum));

        ifs.read(reinterpret_cast<char*>(&needInsertSJ_), sizeof(needInsertSJ_));
        ifs.read(reinterpret_cast<char*>(&sjdbOverhang_), sizeof(sjdbOverhang_));
        ifs.read(reinterpret_cast<char*>(&limitSjdbInsertN_), sizeof(limitSjdbInsertN_));

        // read vectors
        uint64_t lcpSize = 0;
        ifs.read(reinterpret_cast<char*>(&lcpSize), sizeof(lcpSize));
        if (!ifs) return -6;
        longestCommonPrefix_.clear();
        longestCommonPrefix_.resize(static_cast<size_t>(lcpSize));
        if (lcpSize) ifs.read(reinterpret_cast<char*>(longestCommonPrefix_.data()), lcpSize * sizeof(uint8_t));
        if (!ifs) return -7;

        uint64_t extSize = 0;
        ifs.read(reinterpret_cast<char*>(&extSize), sizeof(extSize));
        if (!ifs) return -8;
        extendedIndexHash_.clear();
        extendedIndexHash_.resize(static_cast<size_t>(extSize));
        if (extSize) ifs.read(reinterpret_cast<char*>(extendedIndexHash_.data()), extSize * sizeof(uint32_t));
        if (!ifs) return -9;

        // load complex objects from their files
        {
            std::string gfile = fileName + ".gen";
            if (genome_.loadFromFile(gfile) != 0) return -10;
        }
        {
            std::string sfile = fileName + ".sa";
            if (suffixArray_.loadFromFile(sfile) != 0) return -11;
        }
        {
            std::string kfile = fileName + ".saI";
            if (patternMerMap_.loadFromFile(kfile) != 0) return -12;
        }

        return 0;
    }





}
