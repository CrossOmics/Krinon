#ifndef RNAALIGNREFACTORED_SUFFIXARRAY_H
#define RNAALIGNREFACTORED_SUFFIXARRAY_H
#include <string>
namespace RefactorProcessing {
    class SuffixArray{
    private:
        //configs
        std::string method; // SA-IS or parallel O(n2) or other methods will be added later

        //data
    public:
        SuffixArray(){};
        ~SuffixArray(){};
        void build(const std::string &seq);

        // different building methods
        void buildTraditional(const std::string &seq);
        void buildSAIS(const std::string &seq);

    };
};
#endif //RNAALIGNREFACTORED_SUFFIXARRAY_H
