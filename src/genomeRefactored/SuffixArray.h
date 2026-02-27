#ifndef RNAALIGNREFACTORED_SUFFIXARRAYR_H
#define RNAALIGNREFACTORED_SUFFIXARRAYR_H
#include <string>
#include <vector>
#include "../utilsRefactored/PackedArray.h"
#include "../io/Parameters.h"
namespace RefactorProcessing {
    class SuffixArray{
    private:
        //configs
        std::string buildMethod; // SA-IS or parallel O(n2) or other methods will be added later
        int64_t reservedLength;
        int threadNum;

        // reusable scratch buffers to minimize allocations across calls
        std::vector<int> scratchRank;
        std::vector<size_t> scratchSA;
        std::vector<int> scratchTmpRank;
        std::vector<size_t> scratchTmpSA;


    public:
        //data
        PackedArray suffixArray_; // Suffix Array

        int64_t length_;

        int64_t wordBits_;

        SuffixArray(){
            buildMethod = "SAIS";
            reservedLength = 0;
            threadNum = 0;
        };
        ~SuffixArray(){};
        void build(const std::string &seq);

        void buildFromOtherInit(const SuffixArray &other,int newLength);

        void setParam(const Parameters &P);

        PackedArray& getSuffixArray();

        uint64_t operator[](int index) const { return suffixArray_.getValue(index); }

        // different building methods
        void buildTraditional(const std::string &seq);
        void buildSAIS(const std::string &seq);

        int writeToFile(const std::string &fileName) const;

        int loadFromFile(const std::string &fileName);

    };
};
#endif //RNAALIGNREFACTORED_SUFFIXARRAYR_H
