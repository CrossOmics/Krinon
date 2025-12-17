#ifndef RNAALIGNREFACTORED_GENOMEINDEXR_H
#define RNAALIGNREFACTORED_GENOMEINDEXR_H

#include "Genome.h"
#include "../io/Parameters.h"
#include <vector>
#include "SuffixArray.h"
#include "SJDB.h"
#include "../utilsRefactored/defines.h"
#include "../utilsRefactored/seqFunctions.h"


namespace RefactorProcessing {
    // seed search results
    struct Align{
        int64_t readPos;
        int64_t leftSAIndex;
        int64_t rightSAIndex;
        int64_t length;
        int64_t rep;
        int direction; // 0 for forward, 1 for reverse
        int iFragment; // 0 or 1, for paired-end reads
    };



    class GenomeIndex {
    private:
        // configs

        int kMerSize_;
        std::string additionalIndexType_; // Now only LCP and pre-calculated search results. Reserved for future expansion.

        bool needInsertSJ_;

        int extendAlternativeByte_;
        int sjdbOverhang_;
        int limitSjdbInsertN_;
        // data

        Genome genome_;// owned genome

        SuffixArray suffixArray_;


        //todo replace with PackedIndexArray
        std::vector<uint32_t> patternMerMap_; // stored result of k-mer search

        std::vector<uint8_t> longestCommonPrefix_; // LCP
        std::vector<uint32_t> extendedIndexHash_; // pre-calculated search results to accelerate search


    public:
        GenomeIndex();

        ~GenomeIndex();

        // parameters
        void setParam(const Parameters &P);

        // genome Generate Process
        void loadGenome();

        void build(const Genome &genome);

        void buildKmerMap(const Genome &genome);

        void buildLCP(const Genome &genome);

        // for additional sequences like SJDB
        void modify();

        // search
        void find(const std::string &seq,std::vector<Align> &results) const;

        void findMMP(const std::string &seq) const;

        void findMMP_GetRange(const std::string &seq, int64_t rangeLeft,int64_t rangeRight, size_t matchedLength) const;

        inline std::pair<size_t,bool> matchGenomeSeq(const std::string &pattern ,size_t matchedLength, size_t pos) const;


        // input/output
        int writeToFile(const std::string &fileName) const;

        int loadFromFile(const std::string &fileName);

    };
}

#endif //RNAALIGNREFACTORED_GENOMEINDEXR_H
