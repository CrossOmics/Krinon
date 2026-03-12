#ifndef RNAALIGNREFACTORED_SEQFUNCTIONS_H
#define RNAALIGNREFACTORED_SEQFUNCTIONS_H

#include <string>
#include <fstream>
namespace RefactorProcessing {
    int64_t encodeKmer(const std::string_view &seq, unsigned int kMerSize);

    inline int32_t charToIndex(char c){

            switch (c) {
                case 'A':
                    return 0;
                case 'C':
                    return 1;
                case 'G':
                    return 2;
                case 'T':
                    return 3;
                default:
                    return -1; // For 'N' or any other character
            }

    }

    bool writeString(std::ofstream &ofs, const std::string &s);

    bool readString(std::ifstream &ifs, std::string &s);

    inline bool compSeq(const std::string_view& seq1, const std::string_view& seq2){
        // <=
        size_t l1 = seq1.length();
        size_t l2 = seq2.length();
        size_t l = std::min(l1, l2);
        for (size_t i = 0; i < l; ++i) {
            if(seq1[i] == seq2[i]) continue;
            if (seq1[i] == 'N' || seq1[i] == '#') return false;
            if (seq2[i] == 'N' || seq2[i] == '#') return true;
            if (seq1[i] < seq2[i]) return true;
            if (seq1[i] > seq2[i]) return false;
        }
        return l1 > l2;
    }

}
#endif //RNAALIGNREFACTORED_SEQFUNCTIONS_H
