#include "SuffixArray.h"
#include <omp.h>
#include <algorithm>
#include <vector>
#include <cstring>
#include <cstdint>
#include "../utilsRefactored/seqFunctions.h"

namespace RefactorProcessing {

    void SuffixArray::buildTraditional(const std::string &seq) {
        const int64_t n = seq.length();
        size_t patLen = (size_t)n + 1; // include sentinel


        scratchRank.resize(patLen);
        scratchSA.resize(patLen);
        scratchTmpRank.resize(patLen);


        std::vector<int> &rank = scratchRank;
        std::vector<size_t> &sa = scratchSA;
        std::vector<int> &tmp_rank = scratchTmpRank;

        for (int64_t i = 0; i < n; ++i) {
            int c = charToIndex(seq[(size_t)i]);
            rank[(size_t)i] = (c < 0) ? 5 : c + 1;
        }
        rank[(size_t)n] = 0;


        for (size_t i = 0; i < patLen; ++i) sa[i] = i;


        for (size_t k = 1; k <= (size_t)n; k <<= 1) {
            std::cout<< k << "\n";
            // comparator that uses current rank and k
            auto comp = [&](size_t a, size_t b)->bool {
                if (rank[a] != rank[b]) return rank[a] < rank[b];
                int ra = (a + k <= (size_t)n) ? rank[a + k] : 0;
                int rb = (b + k <= (size_t)n) ? rank[b + k] : 0;
                if (ra != rb) return ra < rb;
                return a < b;
            };


            int numThreads = (threadNum > 0) ? threadNum : 1;


            size_t numBlocks = (size_t)numThreads;
            size_t blockSize = (patLen + numBlocks - 1) / numBlocks;


            std::vector<size_t> &buffer = scratchTmpSA;
            buffer.resize(patLen);


            #pragma omp parallel for schedule(static) num_threads(numThreads)
            for (size_t t = 0; t < numBlocks; ++t) {
                size_t left = t * blockSize;
                size_t right = std::min(left + blockSize, patLen);
                if (left < right) std::sort(sa.begin() + left, sa.begin() + right, comp);
            }


            std::vector<size_t> *src = &sa;
            std::vector<size_t> *dst = &buffer;
            size_t mergeSize = blockSize;
            while (mergeSize < patLen) {
                size_t step = 2 * mergeSize;
                #pragma omp parallel for schedule(static) num_threads(numThreads)
                for (size_t left = 0; left < patLen; left += step) {
                    size_t mid = std::min(left + mergeSize, patLen);
                    size_t right = std::min(left + step, patLen);
                    size_t i = left, j = mid, p = left;
                    while (i < mid && j < right) {
                        if (comp((*src)[i], (*src)[j])) (*dst)[p++] = (*src)[i++];
                        else (*dst)[p++] = (*src)[j++];
                    }
                    while (i < mid) (*dst)[p++] = (*src)[i++];
                    while (j < right) (*dst)[p++] = (*src)[j++];
                }
                // swap src/dst
                std::swap(src, dst);
                mergeSize = step;
            }


            if (src != &sa) {
                std::copy(src->begin(), src->begin() + patLen, sa.begin());
            }


            tmp_rank[sa[0]] = 0;
            for (size_t i = 1; i < patLen; ++i) {
                size_t prev = sa[i - 1];
                size_t cur = sa[i];
                if (rank[prev] != rank[cur]) {
                    tmp_rank[cur] = tmp_rank[prev] + 1;
                } else {
                    int rprev = (prev + k <= (size_t)n) ? rank[prev + k] : 0;
                    int rcur = (cur + k <= (size_t)n) ? rank[cur + k] : 0;
                    tmp_rank[cur] = tmp_rank[prev] + (rprev != rcur ? 1 : 0);
                }
            }


            for (size_t i = 0; i < patLen; ++i) rank[i] = tmp_rank[i];

            // if all ranks are unique, done
            if ((size_t)rank[sa[patLen - 1]] == patLen - 1) break;
        }


        size_t saLen = std::max<int64_t>(n, 256);
        auto* rawSA = new size_t[saLen + reservedLength];
        std::memset(rawSA, 0, (saLen + reservedLength) * sizeof(size_t));
        size_t* raw = rawSA + reservedLength;

        for (size_t i = 0; i < patLen; ++i) raw[i] = sa[i];


        size_t maxValue = (n > 0) ? (size_t)(n - 1) : 0;
        uint32_t bitsNeeded = 0;
        while (maxValue > 0) {
            bitsNeeded++;
            maxValue >>= 1;
        }

        // build PackedArray from raw
        suffixArray_.buildFromReservedArray(rawSA, bitsNeeded, n, reservedLength);

        size_t validCount = 0;
        for (size_t i = 1; i <= (size_t)n; ++i) {
            if (charToIndex(seq[raw[i]]) >= 0) {
                suffixArray_.setValue(validCount, raw[i]);
                ++validCount;
            }
        }
        suffixArray_.setLength(validCount);
        length_ = (int64_t)validCount;
        wordBits_ = bitsNeeded;

        delete[] rawSA;
    }
}
