#include "WindowScanningKMerBuilder.h"
#include <iostream>
#include "GenomeIndex.h"
#include <assert.h>
namespace rna{
    void WindowScanningKMerBuilder::build(GenomeIndex* index) {
        genomeIndex = index;
        buildKMerIndex();
        fillKMerIndex();
        //outputDebugData();
    }


    inline void WindowScanningKMerBuilder::setFlag(int32_t hash,int32_t length){
        if (length == 1){
            appearance_flag[length][0] |= (1 << hash);
        }else {
            int32_t pos = hash >> 3;
            int32_t bit = hash & 0x7;
            appearance_flag[length][pos] |= (1 << bit);
        }
    }
    inline bool WindowScanningKMerBuilder::getFlag(int32_t hash,int32_t length){
        if (length == 1){
            return (appearance_flag[length][0] >> hash) & 1;
        }else {
            int32_t pos = hash >> 3;
            int32_t bit = hash & 0x7;
            return (appearance_flag[length][pos] >> bit) & 1;
        }
    }
    void WindowScanningKMerBuilder::scanGenome() {
        // initialize appearance_flag
        const std::string &seq = genomeIndex->genome->sequence_;
        appearance_flag.resize(genomeIndex->MER_LENGTH+1); // 0 is empty, for convenience
        for (int i = 1; i <= genomeIndex->MER_LENGTH; ++i) {
            appearance_flag[i].resize(1 << (std::max(2*i-3,1)), 0);
        }
        nowWindowHash.resize(genomeIndex->MER_LENGTH+1, 0); // Also, 0 is empty
        int nowAvailableLength = 0;
        for(size_t i = 0; i < seq.length();++i){
            charToIndex(seq[i]) >= 0 ? ++nowAvailableLength : nowAvailableLength = 0;
            if (nowAvailableLength == 0) continue;
            if (nowAvailableLength <= genomeIndex->MER_LENGTH){
                nowWindowHash[nowAvailableLength] = (nowWindowHash[nowAvailableLength-1] << 2) | charToIndex(seq[i]);
                for (int j = 1; j < nowAvailableLength; ++j) {
                    nowWindowHash[j] = (nowWindowHash[j] << 2) | charToIndex(seq[i]);
                    nowWindowHash[j] &= (1 << (2 * j)) - 1; // keep only the last 2*j bits
                }
                for (int j = 1; j <= nowAvailableLength; ++j) {
                    setFlag(nowWindowHash[j], j);
                }
            }else{
                for (int j = 1;j <= genomeIndex->MER_LENGTH;++j){
                    nowWindowHash[j] = (nowWindowHash[j] << 2) | charToIndex(seq[i]);
                    nowWindowHash[j] &= (1 << (2 * j)) - 1; // keep only the last 2*j bits
                    setFlag(nowWindowHash[j], j);
                }
            }

        }
    }
    void WindowScanningKMerBuilder::buildKMerIndex() {
        scanGenome();
        for (size_t i = 0;i< genomeIndex->MER_NUM - 1;++i){
            int32_t hash = i;
            int32_t length = genomeIndex->MER_LENGTH;
            while (!getFlag(hash, length) && length > 0) {
                --length;
                hash &= (1<<(2*length)) - 1; // remove the first two bits
            }
            genomeIndex->patternLongMerMap_[i].length = length;
        }
        for (size_t i = 0;i< genomeIndex->SHORT_MER_NUM - 1;++i){
            int32_t hash = i;
            int32_t length = genomeIndex->SHORT_MER_LENGTH;
            while (!getFlag(hash, length) && length > 0) {
                --length;
                hash &= (1<<(2*length)) - 1; // remove the first two bits
            }
            genomeIndex->patternShortMerMap_[i].length = length;
        }
    }
    void WindowScanningKMerBuilder::fillKMerIndex() {
        auto &genome = genomeIndex->genome;
        auto &suffixArray = genomeIndex->suffixArray;
        auto &patternLongMerMap_ = genomeIndex->patternLongMerMap_;
        auto &patternShortMerMap_ = genomeIndex->patternShortMerMap_;
        int MER_LENGTH = genomeIndex->MER_LENGTH;
        int SHORT_MER_LENGTH = genomeIndex->SHORT_MER_LENGTH;
        size_t MER_NUM = genomeIndex->MER_NUM;
        size_t SHORT_MER_NUM = genomeIndex->SHORT_MER_NUM;
        const std::string &text_ = genome->sequence_;
        const PackedArray & sa_ = suffixArray.sa_;
        for (size_t i = 0; i < sa_.length();++i){
            int64_t hash = encodeKMer(text_.substr(sa_.get(i), MER_LENGTH), MER_LENGTH);
            if (hash!= -1){
                if (patternLongMerMap_[hash].leftSAIndex==-1||patternLongMerMap_[hash].leftSAIndex > i)
                    patternLongMerMap_[hash].leftSAIndex = i;
            }
            int64_t shortHash = encodeKMer(text_.substr(sa_.get(i), SHORT_MER_LENGTH),SHORT_MER_LENGTH);
            if (shortHash != -1) {
                if (patternLongMerMap_[hash].leftSAIndex==-1||patternShortMerMap_[shortHash].leftSAIndex > i)
                    patternShortMerMap_[shortHash].leftSAIndex = i;
            }
        }


        uint64_t last = sa_.length();
        int64_t last_hash = -1;
        patternLongMerMap_[MER_NUM-1].leftSAIndex = last;
        for (int64_t i = MER_NUM-2; i >= 0;--i){
            if (patternLongMerMap_[i].length < MER_LENGTH){
                patternLongMerMap_[i].leftSAIndex = last;
            }else{
                if(last - patternLongMerMap_[i].leftSAIndex > max_14_mer_match ) {
                    max_14_mer_match = last - patternLongMerMap_[i].leftSAIndex;
                    max_14_mer_hash = i;
                }
                last_hash = i;
                last = patternLongMerMap_[i].leftSAIndex;
                assert(last!= -1);
                ++matched_14_mer_num;

            }
        }

        last = sa_.length();
        patternShortMerMap_[SHORT_MER_NUM-1].leftSAIndex = last;
        for (int64_t i = SHORT_MER_NUM-1; i > 0;--i){
            if (patternShortMerMap_[i].length < SHORT_MER_LENGTH){
                patternShortMerMap_[i].leftSAIndex = last;
            }else{
                last = patternShortMerMap_[i].leftSAIndex;
            }
        }
    }

    void WindowScanningKMerBuilder::outputDebugData() const {
        std::cout << "Max 14-mer hash: " << max_14_mer_hash << std::endl;
        std::cout << "Max 14-mer match: " << max_14_mer_match << std::endl;
        std::cout << "Matched 14-mer num: " << matched_14_mer_num << std::endl;

    }
}