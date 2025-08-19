#include <queue>
#include <iostream>
#include "SuffixAutomatonKMerBuilder.h"
#include "GenomeIndex.h"

namespace rna {
    void SuffixAutomatonKMerBuilder::build(GenomeIndex* index) {
        genomeIndex = index;
        buildSAM();
        buildKMerIndex(genomeIndex->patternLongMerMap_,
                       genomeIndex->MER_LENGTH);
        buildKMerIndex(genomeIndex->patternShortMerMap_,
                       genomeIndex->SHORT_MER_LENGTH);
        fillMerIndex();
        outputDebugData();
    }

    void SuffixAutomatonKMerBuilder::buildSAM() {
        const std::string &text_ = genomeIndex->genome->sequence_;
        const size_t n = text_.length();
        sam_.reserve(2 * n);
        sam_.emplace_back(0);

        int64_t last = 0;
        bool former_invalid = false; // 标记前一个字符是否无效，如果是，直接跳过。也即，连续的换行符至多在SAM中保留1个
        for (size_t i = 0; i < n; ++i) {
            char c = text_[i];

            int32_t charIdx = charToIndex(c);
            if (charIdx == -1) {
                if (former_invalid) continue; // 如果前一个字符已经是无效字符，跳过当前字符
                former_invalid = true; // 标记当前字符为无效字符
                charIdx = 4; // 用于避免 AA#GG 在自动机中变成 AAGG。
            } else {
                former_invalid = false; // 有效字符，重置标记
            }

            int64_t cur = sam_.size();
            sam_.emplace_back(sam_[last].len + 1);

            int64_t p = last;
            while (p != -1 && sam_[p].next[charIdx] == -1) {
                sam_[p].next[charIdx] = cur;
                p = sam_[p].link;
            }

            if (p == -1) {
                sam_[cur].link = 0;
            } else {
                int64_t q = sam_[p].next[charIdx];
                if (sam_[p].len + 1 == sam_[q].len) {
                    sam_[cur].link = q;
                } else {
                    int64_t clone = sam_.size();
                    sam_.push_back(sam_[q]);
                    sam_[clone].len = sam_[p].len + 1;

                    sam_[q].link = clone;
                    sam_[cur].link = clone;

                    while (p != -1 && sam_[p].next[charIdx] == q) {
                        sam_[p].next[charIdx] = clone;
                        p = sam_[p].link;
                    }
                }
            }
            last = cur;
        }
    }

    void SuffixAutomatonKMerBuilder::buildKMerIndex(std::vector<MatchResultStored> &index,int merLength){
        index.clear();
        index.resize(1+(1<<(2*merLength)));
        struct queueState{
            int8_t acceptedLength;
            int8_t matchedLength;
            int32_t hash;
        };
        std::queue<std::pair<int64_t, queueState>> q;
        q.push({0, {0,0,0}});
        while (!q.empty()) {
            auto [state, cur] = q.front();
            q.pop();

            if (cur.acceptedLength == merLength) {
                uint32_t hash = cur.hash;
                index[hash].length = cur.matchedLength;
                continue;
            }

            for (int32_t c : {0, 1, 2, 3}) {
                int32_t next = sam_[state].next[c];

                int64_t now = state;
                int8_t acceptedLength = cur.acceptedLength + 1;
                int8_t matchedLength = cur.matchedLength;
                while (next == -1){
                    now = sam_[now].link;
                    next = sam_[now].next[c];
                    matchedLength = sam_[sam_[state].link].len;
                }
                matchedLength++;

                q.push({next, {acceptedLength,matchedLength, (cur.hash << 2) | c}});
            }
        }
    }





    void SuffixAutomatonKMerBuilder::fillMerIndex() {
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
        for (int64_t i = MER_NUM-2; i > 0;--i){
            if (patternLongMerMap_[i].length < MER_LENGTH){
                patternLongMerMap_[i].leftSAIndex = last;
            }else{
                if(last - patternLongMerMap_[i].leftSAIndex > max_14_mer_match ) {
                    max_14_mer_match = last - patternLongMerMap_[i].leftSAIndex;
                    max_14_mer_hash = i;
                }
                last_hash = i;
                last = patternLongMerMap_[i].leftSAIndex;
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

    void SuffixAutomatonKMerBuilder::clear() {
        sam_.clear();
        sam_.shrink_to_fit();
        genomeIndex = nullptr;
    }

    void SuffixAutomatonKMerBuilder::outputDebugData() const {
        std::cout << "Max 14-mer hash: " << max_14_mer_hash << std::endl;
        std::cout << "Max 14-mer match: " << max_14_mer_match << std::endl;
        std::cout << "Matched 14-mer num: " << matched_14_mer_num << std::endl;

    }
}
