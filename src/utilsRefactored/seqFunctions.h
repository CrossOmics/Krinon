#ifndef RNAALIGNREFACTORED_SEQFUNCTIONS_H
#define RNAALIGNREFACTORED_SEQFUNCTIONS_H

#include <string>
namespace RefactorProcessing {
    int64_t encodeKmer(const std::string_view& seq, int kMerSize);

    int32_t charToIndex(char c);
}

#endif //RNAALIGNREFACTORED_SEQFUNCTIONS_H
