#include "defines.h"
#include "seqFunctions.h"

namespace RefactorProcessing {

    // encode a k-mer into an integer hash, return -1 if the seq length is less than kMerSize, or -i-2 if the i-th character is invalid
    int64_t encodeKmer(const std::string_view &seq, int kMerSize) {
        if (seq.length() < kMerSize) {
            return -1;
        }
        int64_t hash = 0;
        for (int i = 0; i < kMerSize; ++i) {
            int32_t idx = charToIndex(seq[i]);
            if (idx < 0) return -i - 2;
            hash = (hash << 2) | idx;
        }
        return hash;
    }

    int32_t charToIndex(char c) {
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
}
