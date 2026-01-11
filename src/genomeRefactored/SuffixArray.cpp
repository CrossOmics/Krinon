#include "../utilsRefactored/defines.h"
#include "SuffixArray.h"
namespace RefactorProcessing{
    void SuffixArray::build(const std::string &seq) {
        buildSAIS(seq);
    }

    void SuffixArray::setParam(const Parameters &P){
        reservedLength = P.limitSjdbInsertN;
    }

    PackedArray& SuffixArray::getSuffixArray(){
        return suffixArray_;
    }


}
