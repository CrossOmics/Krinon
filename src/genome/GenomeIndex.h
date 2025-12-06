#ifndef RNAALIGNREFACTORED_GENOMEINDEX_H
#define RNAALIGNREFACTORED_GENOMEINDEX_H
#include "Genome.h"
#include "SuffixArray.h"
#include "SJDB.h"
#include "../utils/PackedIndexArray.h"

#define KMER_INDEX_LENGTH_BITS 4  // max support k-mer length is 15
#define KMER_INDEX_UPPER_RANGE_BITS 24 // max support upper range is 2^24-1 = 16777215


#define KMER_INDEX_LENGTH 0
#define KMER_INDEX_UPPER_RANGE 1
#define KMER_INDEX_LEFT_SA_INDEX 2


namespace rna{
    struct GenomeIndexConfig{
        int kMerSize = 10;
        bool twoDirections {true};
        int extendAlternativeByte = 64;// 16 * 4 byte


        bool insertSJ{false};
        bool insertSecondPass{false};

        int sjdbOverhang{100};
        int sjdbLength{220};
        int limitSjdbInsertN{20000};// the default in STAR is 1000000, 20000 is for smaller test.

    };



    class GenomeIndex {
    private:
        std::vector<std::vector<int8_t>> appearance_flag; // records the appearance of each k-mer in the genome
        std::vector<int32_t> nowWindowHash;



    public:
        int MER_LENGTH{0};  // k-mer length
        size_t MER_NUM{0}; // k-mer num
        int64_t genomeLength{0}; // genome length, one strand

        GenomeIndexConfig config;
        std::unique_ptr<Genome> genome;


        SuffixArray suffixArray;
        SuffixArray suffixArrayFirstPass;
        //SuffixArray suffixArraySecondPass;
        //SuffixArray& suffixArray ;

        GTF *gtf{nullptr};

        PackedIndexArray patternMerMap_;

        //std::vector<MatchResultStored> patternMerMap_;  // stored result of k-mer search
        std::vector<uint8_t> longestCommonPrefix_; // actually, the height array of the suffix array, LCP (i,i-1)




        //accelerate the search
        std::vector<uint32_t> extendIndexHash;



        GenomeIndex() = default;
        GenomeIndex(Genome& g);
        void setConfig(const GenomeIndexConfig& cfg);
        void build();
        void buildKMerIndex();
        void buildLCP();
        void fillKMerIndex();
        void scanGenome();

        inline void setFlag(int32_t hash,int32_t length);
        inline bool getFlag(int32_t hash,int32_t length);
        void write(const std::string& outDir) const;
        void load(const std::string& inDir);
        //void transformSJDB(const std::string& gtfFileName);


        int alignRead(const std::string* readSeq,const Split& split,Align* aligns,int &alignNum,int maxAlignPerRead) const;


        void buildKMerMap();
        void buildKMerMapSingle(int64_t hash);
        MatchStats find(const std::string_view & pattern) const;


        // given a pattern, find the longest match and its rep in the genome
        MatchStats getMatchStatsLCP(const std::string_view & pattern,
                                 size_t rangeLeft,
                                 size_t rangeRight,size_t matchedLength) const;

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

        friend GTF;


    };
}
#endif //RNAALIGNREFACTORED_GENOMEINDEX_H
