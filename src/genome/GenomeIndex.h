#ifndef RNAALIGNMENT_GENOMEINDEX_H
#define RNAALIGNMENT_GENOMEINDEX_H
#include "Genome.h"
#include "SuffixArray.h"
#include "KMerBuilder.h"
#include <memory>
namespace rna{
    struct GenomeIndexConfig{
        enum KMerConstructMethod {
            FULL_SAM,
            WINDOW_SCANNING
        } method = WINDOW_SCANNING;
        int kMerSize = 7;
        int shortKMerSize = 5;
    };


    class GenomeIndex {
    private:

    public:
        int MER_LENGTH{0};  // 14-mer长度
        int SHORT_MER_LENGTH{0};  // 10-mer长度
        size_t MER_NUM{0}; // 14-mer的数量
        size_t SHORT_MER_NUM{0}; // 10-mer的数量

        GenomeIndexConfig config;
        std::unique_ptr<const Genome> genome;
        std::unique_ptr<KMerIndexBuilder> kMerBuilder; // k-mer索引构建器
        SuffixArrayNew suffixArray;
        std::vector<MatchResultStored> patternLongMerMap_;  // 14-mer的查找结果
        std::vector<MatchResultStored> patternShortMerMap_;  // 10-mer的查找结果





        GenomeIndex(Genome& g);
        void setConfig(const GenomeIndexConfig& cfg);
        void build();

        MatchStats find(const std::string& pattern) const;
        MatchStats findShortMatch(const std::string& pattern) const;

        // 获取指定范围内的最大匹配长度及其出现次数
        MatchStats getMatchStats(const std::string& pattern,
                                 size_t rangeLeft,
                                 size_t rangeRight,size_t matchedLength) const;





        static int64_t encodeKMer(const std::string_view & pattern,const int MER_LENGTH) ;


        static constexpr int32_t charToIndex(char c) {
            switch(c) {
                case 'A': return 0;
                case 'C': return 1;
                case 'G': return 2;
                case 'T': return 3;
                default: return -1;
            }
        }


    };
}

#endif //RNAALIGNMENT_GENOMEINDEX_H
