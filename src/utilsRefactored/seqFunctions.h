#ifndef RNAALIGNREFACTORED_SEQFUNCTIONS_H
#define RNAALIGNREFACTORED_SEQFUNCTIONS_H

#include <string>
#include <fstream>
namespace RefactorProcessing {
    int64_t encodeKmer(const std::string_view& seq, int kMerSize);

    int32_t charToIndex(char c);

    bool writeString(std::ofstream &ofs, const std::string &s);

    bool readString(std::ifstream &ifs, std::string &s);
}

#endif //RNAALIGNREFACTORED_SEQFUNCTIONS_H
