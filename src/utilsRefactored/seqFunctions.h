#ifndef RNAALIGNREFACTORED_SEQFUNCTIONS_H
#define RNAALIGNREFACTORED_SEQFUNCTIONS_H

#include <string>
#include <fstream>

namespace RefactorProcessing {
    int64_t encodeKmer(const std::string_view &seq, unsigned int kMerSize);


    // fast lookup table for nucleotide ranking used in comparisons
    static inline const int8_t *getCharRankTable() {
        static int8_t table[256];
        static bool inited = []() {
            for (int i = 0; i < 256; ++i) table[i] = -1;
            table[(unsigned char) 'A'] = 0;
            table[(unsigned char) 'C'] = 1;
            table[(unsigned char) 'G'] = 2;
            table[(unsigned char) 'T'] = 3;
            table[(unsigned char) 'a'] = 0;
            table[(unsigned char) 'c'] = 1;
            table[(unsigned char) 'g'] = 2;
            table[(unsigned char) 't'] = 3;
            // leave other characters (N, #, etc.) as -1
            return true;
        }();
        (void) inited;
        return table;
    }

    inline int32_t charToIndex(char c) {

        return getCharRankTable()[(unsigned char) c];

    }

    bool writeString(std::ofstream &ofs, const std::string &s);

    bool readString(std::ifstream &ifs, std::string &s);

    inline bool compSeq(const std::string_view &seq1, const std::string_view &seq2) {
        // <=
        size_t l1 = seq1.length();
        size_t l2 = seq2.length();
        size_t l = std::min(l1, l2);
        for (size_t i = 0; i < l; ++i) {
            if (seq1[i] == seq2[i]) continue;
            if (seq1[i] == 'N' || seq1[i] == '#') return false;
            if (seq2[i] == 'N' || seq2[i] == '#') return true;
            if (seq1[i] < seq2[i]) return true;
            if (seq1[i] > seq2[i]) return false;
        }
        return l1 > l2;
    }

}
#endif //RNAALIGNREFACTORED_SEQFUNCTIONS_H
