#ifndef RNAALIGNREFACTORED_SUFFIXAUTOMATONKMERBUILDER_H
#define RNAALIGNREFACTORED_SUFFIXAUTOMATONKMERBUILDER_H
#include "KMerBuilder.h"
#include "../utils/types.h"
#include <array>


namespace rna{
    class SuffixAutomatonKMerBuilder : public KMerIndexBuilder {
    private:
        struct SAMNode {
            GenomePos len;     // 节点代表的所有字符串的长度
            int64_t link;      // 后缀链接（指向最长的等价后缀）
            std::array<int64_t, 5> next;  // ACGT 4个方向的转移和非ACGT字符的转移

            SAMNode(GenomePos l = 0) : len(l), link(-1) {
                next.fill(-1);  // -1表示没有转移
            }
        };

        static constexpr int32_t charToIndex(char c) {
            switch(c) {
                case 'A': return 0;
                case 'C': return 1;
                case 'G': return 2;
                case 'T': return 3;
                default: return -1;
            }
        }

        std::vector<SAMNode> sam_;

        //debug data
        int64_t max_14_mer_hash{-1};
        int64_t max_14_mer_match{0};
        int64_t matched_14_mer_num{0};

        void buildSAM();
        void buildKMerIndex(std::vector<MatchResultStored>& index,int merLength);
        void fillMerIndex();
        void clear();
        void outputDebugData() const;

    public:

        void build(GenomeIndex* index) override;
    };
}

#endif //RNAALIGNREFACTORED_SUFFIXAUTOMATONKMERBUILDER_H
