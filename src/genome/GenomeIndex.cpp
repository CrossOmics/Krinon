#include "GenomeIndex.h"
#include "KMerBuilder.h"
#include "SuffixAutomatonKMerBuilder.h"
#include "WindowScanningKMerBuilder.h"
#include <queue>
namespace rna {
    GenomeIndex::GenomeIndex(Genome& g) : genome(std::make_unique<const Genome>(g)) {
        suffixArray.build(genome->sequence_);
        kMerBuilder = nullptr;
        setConfig(GenomeIndexConfig());
    }

    void GenomeIndex::setConfig(const rna::GenomeIndexConfig &cfg) {
        config = cfg;
        MER_LENGTH = config.kMerSize;
        SHORT_MER_LENGTH = config.shortKMerSize;
        MER_NUM = (1<<(MER_LENGTH * 2))+1;
        SHORT_MER_NUM = (1<<(SHORT_MER_LENGTH * 2))+1;
    }

    void GenomeIndex::build() {
        patternLongMerMap_.resize(MER_NUM);
        patternShortMerMap_.resize(SHORT_MER_NUM);
        if (config.method == GenomeIndexConfig::FULL_SAM) {
            kMerBuilder = std::make_unique<SuffixAutomatonKMerBuilder>();
        } else if (config.method == GenomeIndexConfig::WINDOW_SCANNING) {
            kMerBuilder = std::make_unique<WindowScanningKMerBuilder>();
        }
        kMerBuilder->build(this);
        //genome->writeChrInfo("test");
    }


    int64_t GenomeIndex::encodeKMer(const std::string_view &pattern,const int mer_length) {
        if (pattern.length() < mer_length) return -1;
        uint32_t code = 0;
        for (char c : pattern.substr(0, mer_length)) {
            code = (code << 2) | charToIndex(c);
        }
        return code;
    }


    MatchStats GenomeIndex::find(const std::string& pattern) const {
        if (pattern.length() < MER_LENGTH) {
            return {0,0,0,0};
        }
        uint32_t code = encodeKMer(pattern, MER_LENGTH);
        if(patternLongMerMap_[code].length < MER_LENGTH) {
            return {patternLongMerMap_[code].length,0,0,0};
        }else{
            return getMatchStats(pattern,
                                 patternLongMerMap_[code].leftSAIndex,
                                 patternLongMerMap_[code+1].leftSAIndex,MER_LENGTH);
        }
    }

    MatchStats GenomeIndex::findShortMatch(const std::string& pattern) const {
        if (pattern.length() < SHORT_MER_LENGTH) {
            return {0,0,0,0};
        }

        uint32_t code = encodeKMer(pattern,SHORT_MER_LENGTH);
        if(patternShortMerMap_[code].length < MER_LENGTH) {
            return {patternShortMerMap_[code].length,0,0,0};
        }else{
            return getMatchStats(pattern,
                                 patternShortMerMap_[code].leftSAIndex,
                                 patternShortMerMap_[code+1].leftSAIndex,SHORT_MER_LENGTH);
        }
    }

    inline size_t matchLength(size_t pos, const std::string& text, const std::string& pattern,size_t matchedLength){
        size_t maxLen = 0;
        for (size_t i = matchedLength; i < pattern.length() && pos + i < text.length(); ++i) {
            if (text[pos + i] != pattern[i]) {
                break;
            }
            maxLen++;
        }
        return matchedLength + maxLen;
    }

    MatchStats GenomeIndex::getMatchStats(const std::string& pattern,
                                          size_t rangeLeft,
                                          size_t rangeRight,size_t matchedLength) const
    {
        const std::string &text_ = genome->sequence_;
        const PackedArray &sa_ = suffixArray.sa_;
        if (rangeLeft >= rangeRight) {
            return {};
        }

        // 计算第一个位置的匹配长度作为初始值
        size_t pos = sa_.get(rangeLeft);
        size_t leftLength = matchLength(pos, text_, pattern, matchedLength);

        if (rangeRight - rangeLeft > 1){
            size_t mid = (rangeLeft + rangeRight) / 2;
            pos = sa_.get(mid);
            size_t midLength = matchLength(pos, text_, pattern, matchedLength);
            size_t rightPos = sa_.get(rangeRight - 1);
            size_t rightLength = matchLength(rightPos, text_, pattern, matchedLength);
            if (leftLength > rightLength){
                if (midLength < leftLength) {
                    //此时只需找一侧，另一边可匹配的长度必定小于leftLength
                    return getMatchStats(pattern, rangeLeft, mid, midLength);;
                }
            }else if(leftLength <rightLength){
                if (midLength < rightLength) {
                    return getMatchStats(pattern, mid, rangeRight, midLength);;
                }
            }
            auto leftMatch = getMatchStats(pattern, rangeLeft, mid, leftLength);
            auto rightMatch = getMatchStats(pattern, mid, rangeRight, rightLength);
            if (leftMatch.maxLength < rightMatch.maxLength) {
                return rightMatch;
            }else if(leftMatch.maxLength > rightMatch.maxLength) {
                return leftMatch;
            }else{
                return {static_cast<size_t>(leftMatch.maxLength), leftMatch.count + rightMatch.count,
                        leftMatch.LongestMatchLeftPosInSA, rightMatch.LongestMatchRightPosInSA};
            }
        }
        return {
                leftLength,
                1,
                rangeLeft,
                rangeRight - 1
        };

    }


}