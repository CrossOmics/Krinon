#ifndef RNAALIGNREFACTORED_SUFFIXARRAY_H
#define RNAALIGNREFACTORED_SUFFIXARRAY_H
#pragma once

#include "../utils/types.h"
#include "../utils/PackedArray.h"
#include <vector>
#include <cstdint>
#include <array>
#include <unordered_map>
#include <variant>
#include <string>
#include <utility>
#include <limits>
#include <cstring>

namespace rna {


    class SuffixArray {
        // a less memory consuming version of SuffixArray
    public:
        SuffixArray() = default;

        void build(const std::string &seq, bool skipInvalidSuffixes = true, size_t reservedLength = 0);

        int64_t findInsertPosition(const std::string_view &pattern,const std::string_view& genomeSeq) const;

        inline GenomePos operator[](size_t index) const { return sa_.get(index); }

        size_t length() const { return sa_.length(); }

        PackedArray sa_; //compressed suffix array
        int64_t fullLength_; // full length of the original text, including invalid characters
    private:
        static constexpr int32_t charToIndex(char c) {
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
                    return 4;
            }
        }

        std::vector<unsigned char> text_;         // original text
        static constexpr size_t EMPTY = std::numeric_limits<int64_t>::max() / 2 + 1;
        static constexpr size_t UNIQUE = EMPTY + 1;
        static constexpr size_t MULTI = EMPTY + 2;

        inline bool patCharType(size_t cur, size_t prev, bool lastScanned) {
            return (cur < prev) || (cur == prev && lastScanned);
        }

// 比较两个 LMS 子串
        int lmsStrCmp(const size_t *s, size_t p1, size_t p2, size_t len);

        void renamePat(size_t *pat, size_t *sa, size_t patLen, size_t saLen);

        size_t sortLmsChar(size_t *pat, size_t *sa, size_t patLen, size_t saLen);

        void sortLmsSubstr(size_t *pat, size_t *sa, size_t patLen, size_t saLen);

        bool constructPat1(size_t *pat, size_t *sa, size_t lmsCnt, size_t patLen, size_t saLen);

        void sortLmsSuf(size_t *pat, size_t *sa, size_t patLen, size_t saLen, size_t lmsCnt, bool dup);

        void inducedSort(size_t *pat, size_t *sa, size_t patLen, size_t saLen);

        void computeSuffixArray(size_t *pat, size_t *sa, size_t patLen, size_t saLen);

        size_t *suffixArray(const std::vector<uint8_t> &pat,size_t reservedLength);


    };

} // namespace rna
#endif //RNAALIGNREFACTORED_SUFFIXARRAY_H
