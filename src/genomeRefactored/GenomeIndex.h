#ifndef RNAALIGNREFACTORED_GENOMEINDEXR_H
#define RNAALIGNREFACTORED_GENOMEINDEXR_H

#include "Genome.h"
#include "../io/Parameters.h"
#include <vector>
#include "SuffixArray.h"
#include "SuffixArrayKMerMap.h"
#include "SJDB.h"
#include "../utilsRefactored/defines.h"
#include "../utilsRefactored/seqFunctions.h"


namespace RefactorProcessing {
    // seed search results
    struct Align {
        int64_t readPos;
        int64_t leftSAIndex;
        int64_t rightSAIndex;
        int64_t length;
        int64_t rep;
        int direction; // 0 for forward, 1 for reverse
        int iFragment; // 0 or 1, for paired-end reads
        Align() { length = 0; rep = 0;}
    };

    struct Split {
        int splitStart;//position in read, forward strand
        int length;
        int readLength;
        std::string_view forward;
        std::string_view reverse;
    };


    class GenomeIndex {
    private:
        // configs

        int kMerSize_;
        uint64_t kMerNum_;
        std::string additionalIndexType_; // Now only LCP and pre-calculated search results. Reserved for future expansion.
        int extendHashTableByte_; // Bytes of each hash's extend table
        int extendHashTableNum_; // Number of elements in each hash's extend table

        //todo add to Parameters
        int maxAlignNum = 1000; // max alignments to return for one search



        bool needInsertSJ_;
        int sjdbOverhang_;
        int limitSjdbInsertN_;
        // data

        Genome genome_;// owned genome

        SuffixArray suffixArray_;


        //todo replace with PackedIndexArray
        SuffixArrayKMerMap patternMerMap_; // stored result of k-mer search

        std::vector<uint8_t> longestCommonPrefix_; // LCP
        std::vector<uint32_t> extendedIndexHash_; // pre-calculated search results to accelerate search



        inline bool insertAlign(std::vector<Align> &results, const Align &a) const;

    public:


        GenomeIndex();

        ~GenomeIndex();

        // parameters
        void setParam(const Parameters &P);

        // genome Generate Process
        void loadGenome();

        void build();

        void buildKmerMap();

        void buildLCP();

        void buildExtendedIndexHashSingle(int64_t h);

        void buildExtendedIndexHash();

        // for additional sequences like SJDB
        void modify();

        // search
        void find(const Split pattern, std::vector<Align> &results) const;

        Align findMMP(const std::string_view &seq) const;

        Align findMMP_GetRange(const std::string_view &seq, int64_t rangeLeft, int64_t rangeRight, size_t matchedLength) const;

        inline std::pair<size_t, bool>
        matchGenomeSeq(const std::string_view &pattern, size_t matchedLength, size_t pos) const;


        // input/output
        int writeToFile(const std::string &fileName) const;

        int loadFromFile(const std::string &fileName);

    };
}

#endif //RNAALIGNREFACTORED_GENOMEINDEXR_H
