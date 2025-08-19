#ifndef RNAALIGNREFACTORED_GENOMEINDEXPREFIX_H
#define RNAALIGNREFACTORED_GENOMEINDEXPREFIX_H
#include "Genome.h"
#include "SuffixArray.h"

namespace rna{
    struct GenomeIndexPrefixConfig{
        int kMerSize = 4;
        int extendLength = 2;
        int extendIndexSize = 17; // (1 << (2 * extendLength)) + 1
        int minExtendRep = 10; // minimum number of repetitions for extension
        int maxLayer = 3; // maximum number of layers for extension
        bool twoDirections {false};
        int extendAlternativeByte = 64;// 16 * 4 byte

    };



    class GenomeIndexPrefix {
    private:
        std::vector<std::vector<int8_t>> appearance_flag; // records the appearance of each k-mer in the genome
        std::vector<int32_t> nowWindowHash;
        struct extendKMerMap{
            std::vector<MatchResultStored> extendedKMerMap; // extended k-mer map
            std::vector<int64_t> extendedKMerIndex; // index of the extended k-mer map, -1 means no extension
        };


    public:
        int MER_LENGTH{0};  // k-mer length
        size_t MER_NUM{0}; // k-mer num
        int64_t genomeLength{0}; // genome length, one strand

        GenomeIndexPrefixConfig config;
        std::unique_ptr<Genome> genome;
        SuffixArrayNew suffixArray;
        std::vector<MatchResultStored> patternMerMap_;  // stored result of k-mer search
        std::vector<uint16_t> longestCommonPrefix_; // actually, the height array of the suffix array, LCP (i,i-1)
        std::vector<int64_t> extendedKMerIndex_;
        int extendedKMerNum{};
        std::vector<int64_t> nextLayerIndex_;
        std::vector<MatchResultStored> extendedKMerMap_;

        // another way to accelerate the search
        std::vector<uint32_t> extendIndexHash;



        GenomeIndexPrefix() = default;
        GenomeIndexPrefix(Genome& g);
        void setConfig(const GenomeIndexPrefixConfig& cfg);
        void build();
        void buildKMerIndex();
        void buildLCP();
        void fillKMerIndex();
        void scanGenome();
        void buildExtendedKMerMap();



        int64_t buildExtendedKMerMapSingleIndex(int32_t length, int64_t leftSAIndex,int64_t rightSAIndex,int depth = 0);
        inline void setFlag(int32_t hash,int32_t length);
        inline bool getFlag(int32_t hash,int32_t length);
        void write(const std::string& outDir) const;
        void load(const std::string& inDir);

        MatchStats find(const std::string_view & pattern) const;
        std::vector<Align> alignRead(const std::string_view & readSeq,const Split& split) const;


        void buildAlternativeKMerMap();
        MatchStats findAlternative(const std::string_view & pattern) const;


        // given a pattern, find the longest match and its rep in the genome
        MatchStats getMatchStatsLCP(const std::string_view & pattern,
                                 size_t rangeLeft,
                                 size_t rangeRight,size_t matchedLength) const;
        inline int64_t findSameMatchingLength(size_t targetLength,int64_t lBound, int64_t rBound,size_t matchedLength, const std::string_view &pattern) const;

        //return the longest match length and whether the pattern is larger than the subsequence of the genome
        inline std::pair<size_t, bool> matchGenomeSeq(size_t matchedLength, size_t pos, const std::string_view &pattern) const;


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
#endif //RNAALIGNREFACTORED_GENOMEINDEXPREFIX_H
