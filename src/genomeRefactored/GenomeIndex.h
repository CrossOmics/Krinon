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
        size_t readPos;
        size_t leftSAIndex;
        size_t rightSAIndex;
        size_t length;
        int64_t rep;
        int direction; // 0 for forward, 1 for reverse
        int iFragment; // 0 or 1, for paired-end reads
        bool isAnchor; // whether this align is an anchor
        Align() {
            length = 0;
            rep = 0;
        }
    };

    struct Split {
        int splitStart;//position in read, forward strand
        int length;
        int readLength;
        std::string forward;
        std::string reverse;
        int iFragment; // 0 or 1, for paired-end reads
    };


    class GenomeIndex {
    private:
        // configs

        uint32_t kMerSize_;
        uint64_t kMerNum_;
        std::string additionalIndexType_; // Now only LCP and pre-calculated search results. Reserved for future expansion.
        int extendHashTableByte_; // Bytes of each hash's extend table
        int extendHashTableNum_; // Number of elements in each hash's extend table

        //todo add to Parameters
        uint64_t maxAlignNum = 1000; // max alignments to return for one search



        bool needInsertSJ_;
        int sjdbOverhang_;
        int limitSjdbInsertN_;


        //todo replace with PackedIndexArray
        SuffixArrayKMerMap patternMerMap_; // stored result of k-mer search

        std::vector<uint8_t> longestCommonPrefix_; // LCP
        std::vector<uint32_t> extendedIndexHash_; // pre-calculated search results to accelerate search



        inline bool insertAlignResults(std::vector<Align> &results, const Align &a) const;


        inline size_t sjFindInsertPosition(const std::string_view& s) const{

            int64_t left = 0;
            int64_t right = suffixArray_.length_;
            int64_t mid;
            std::string_view gs(genome_.sequence_);
            while (right - left > 1) {
                mid = (left + right) / 2;

                if (compSeq(s, gs.substr(suffixArray_[mid]))) right = mid;
                else left = mid;
            }
            return right;

        }

    public:
        // data

        Genome genome_;// owned genome

        SuffixArray suffixArray_;

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
        void modify(SJDB &sjdb);

        // search
        void find(const Split pattern, std::vector<Align> &results) const;

        Align findMMP(const std::string_view &seq) const;

        Align findMMP_GetRange(const std::string_view &seq, int64_t rangeLeft, int64_t rangeRight,
                               size_t matchedLength) const;

        inline std::pair<size_t, bool>
        matchGenomeSeq(const std::string_view &pattern, size_t matchedLength, size_t pos) const;


        // input/output
        int writeToFile(const std::string &fileName) const;

        int loadFromFile(const std::string &fileName);

    };
}

#endif //RNAALIGNREFACTORED_GENOMEINDEXR_H
