#ifndef RNAALIGNREFACTORED_SUFFIXARRAYR_H
#define RNAALIGNREFACTORED_SUFFIXARRAYR_H
#include <string>
#include "../utilsRefactored/PackedArray.h"
#include "../io/Parameters.h"
namespace RefactorProcessing {
    class SuffixArray{
    private:
        //configs
        std::string buildMethod; // SA-IS or parallel O(n2) or other methods will be added later

        //data
        PackedArray suffixArray_; // Suffix Array

    public:
        SuffixArray(){};
        ~SuffixArray(){};
        void build(const std::string &seq);

        // different building methods
        void buildTraditional(const std::string &seq);
        void buildSAIS(const std::string &seq);

    };
};
#endif //RNAALIGNREFACTORED_SUFFIXARRAYR_H
