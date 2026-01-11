#include <cstring>
#include "SuffixArray.h"
#include "../utilsRefactored/defines.h"
#include "../utilsRefactored/seqFunctions.h"


namespace RefactorProcessing{
    constexpr size_t EMPTY = std::numeric_limits<int64_t>::max() / 2 + 1;
    constexpr size_t UNIQUE = EMPTY + 1;
    constexpr size_t MULTI = EMPTY + 2;
    bool patCharType(size_t cur, size_t prev, bool lastScanned) {
        return (cur < prev) || (cur == prev && lastScanned);
    }

    int lmsStrCmp(const size_t *s, size_t p1, size_t p2, size_t len){
        for (size_t i = 0; i < len; ++i) {
            if (s[p1 + i] != s[p2 + i]) return s[p1 + i] < s[p2 + i] ? -1 : 1;
        }
        return 0;
    }

    void renamePat(size_t *pat, size_t *sa, size_t patLen, size_t saLen){
        for (size_t i = 0; i < saLen; ++i) sa[i] = 0;
        for (size_t i = 0; i < patLen; ++i) ++sa[pat[i]];
        for (size_t i = 1; i < saLen; ++i) sa[i] += sa[i - 1];

        for (size_t i = 0; i < patLen - 1; ++i) pat[i] = sa[pat[i]] - 1;
        for (size_t i = 0; i < saLen; ++i) sa[i] = 0;
        for (size_t i = 0; i < patLen; ++i) ++sa[pat[i]];

        bool lastType = true;
        pat[patLen - 1] = 0;
        for (int64_t i = patLen - 2; i >= 0; --i) {
            if (patCharType(pat[i], pat[i + 1], lastType)) {
                lastType = true;
            } else {
                pat[i] -= sa[pat[i]] - 1;
                lastType = false;
            }
        }
    }

    size_t sortLmsChar(size_t *pat, size_t *sa, size_t patLen, size_t saLen){
        for (size_t i = 0; i < saLen; ++i) sa[i] = EMPTY;

        bool lastType = true;
        for (int64_t i = patLen - 2; i >= 0; --i) {
            if (patCharType(pat[i], pat[i + 1], lastType)) {
                lastType = true;
            } else {
                if (lastType) ++sa[pat[i + 1]];
                lastType = false;
            }
        }

        size_t lmsCnt = 0;
        lastType = true;
        for (int64_t i = patLen - 2; i >= 0; --i) {
            if (patCharType(pat[i], pat[i + 1], lastType)) {
                lastType = true;
            } else {
                size_t e = pat[i + 1];
                if (lastType) {
                    ++lmsCnt;
                    if (sa[e] == UNIQUE) {
                        sa[e] = i + 1;
                    } else if (sa[e] >= MULTI) {
                        if (sa[e-1] == EMPTY) {
                            if (sa[e-2] == EMPTY) {
                                sa[e-2] = i + 1;
                                sa[e-1] = 1;
                            }else {
                                sa[e] = i + 1;
                                sa[e-1] = EMPTY;
                            }
                        } else {
                            size_t c = sa[e-1];
                            if (sa[e-2-c] == EMPTY) {
                                sa[e-2-c] = i + 1;
                                sa[e-1]++;
                            } else {
                                for (int64_t j = c; j >= 1;--j){
                                    sa[e-c+j] = sa[e-c+j-2];
                                }
                                sa[e-c] = i + 1;
                                sa[e-c-1] = EMPTY;
                            }
                        }
                    } else if (sa[e] < EMPTY) {
                        for (int64_t j = e-1; j>=0; --j) {
                            if(sa[j] == EMPTY){
                                sa[j] = i+1;
                                break;
                            }
                        }
                    }
                }
                lastType = false;
            }
        }

        for (int64_t i = patLen - 1; i >= 0; --i) {
            if (sa[i] >= MULTI) {
                int64_t c = sa[i - 1];
                for (int64_t j = c; j >= 1; --j) sa[i - c + j] = sa[i- c + j - 2];
                sa[i - c - 1] = EMPTY;
                sa[i - c] = EMPTY;
            }
        }
        return lmsCnt;
    }

    void inducedSort(size_t *pat, size_t *sa, size_t patLen, size_t saLen);
    void sortLmsSubstr(size_t *pat, size_t *sa, size_t patLen, size_t saLen){
        inducedSort(pat, sa,patLen,saLen);

        size_t lmsCnt = 0;
        size_t tailPtr = patLen;
        size_t bucket = EMPTY;
        size_t numS = 0;

        for (int64_t i = patLen - 1; i > 0; --i) {
            size_t v = pat[sa[i]];
            if (v != bucket) {
                numS = 0;

                size_t l = 0;
                while (pat[sa[i - l]] == pat[sa[i]]) {
                    size_t pos = sa[i - l];
                    if (pat[pos] < pat[pos + 1]) {
                        size_t k = pos;
                        while (k > 0 && pat[k - 1] == pat[pos]) --k;
                        numS += pos - k + 1;
                    } else {
                        break;
                    }
                    ++l;
                }
                tailPtr = i;
                bucket = pat[sa[tailPtr]];
            }

            if (numS > 0 && i > tailPtr - numS && sa[i] > 0 && pat[sa[i]] < pat[sa[i] - 1]) {
                sa[patLen - 1 - lmsCnt] = sa[i];
                ++lmsCnt;
            }
        }

        sa[patLen - 1 - lmsCnt] = sa[0];
        ++lmsCnt;
        for (size_t i = 0; i < patLen - lmsCnt; ++i ) sa[i] = EMPTY;
    }

    bool constructPat1(size_t *pat, size_t *sa, size_t lmsCnt, size_t patLen, size_t saLen){
        size_t rank = 0;
        size_t prevLen = 1;
        bool dup = false;
        sa[(patLen - 1) / 2] = 0;

        for (size_t i = patLen - lmsCnt + 1; i < patLen; ++i) {
            size_t j = sa[i];
            while (pat[j] <= pat[j + 1]) ++j;
            size_t k = j;
            while (k + 1 < patLen && pat[k] >= pat[k + 1]) ++k;
            size_t curLen = k + 1 - sa[i];

            int cmp =  lmsStrCmp(pat, sa[i], sa[i - 1], std::min(curLen, prevLen));
            if (cmp != 0) ++rank;

            if (rank == sa[sa[i - 1] / 2]) dup = true;
            sa[sa[i] / 2] = rank;
            prevLen = curLen;
        }

        size_t w = 0;
        for (size_t i = 0; i < patLen - lmsCnt; ++i) {
            if (sa[i] != EMPTY){
                sa[w] = sa[i];
                if (i > w) {
                    sa[i] = EMPTY;
                }
                ++w;
            }
        }
        for (size_t i =lmsCnt; i < patLen; ++i) sa[i] = EMPTY;

        return dup;
    }

    void computeSuffixArray(size_t *pat, size_t *sa, size_t patLen, size_t saLen);
    void sortLmsSuf(size_t *pat, size_t *sa, size_t patLen, size_t saLen, size_t lmsCnt, bool dup){
        size_t* pat1 = sa;
        size_t pat1Len = lmsCnt;
        size_t* sa1 = sa + patLen - lmsCnt;
        size_t sa1Len = saLen - patLen + lmsCnt;

        if (dup) {
            computeSuffixArray(pat1, sa1, pat1Len, sa1Len);
        } else {
            for (size_t i = 0; i < lmsCnt; ++i) sa1[pat1[i]] = i;
        }

        for (size_t i = 0; i < lmsCnt; ++i) sa[i] = sa[patLen - lmsCnt + i];

        size_t j = 0;
        bool lastType = true;
        for (int64_t i = patLen - 2; i >= 0; --i) {
            if (patCharType(pat[i], pat[i + 1], lastType)) {
                lastType = true;
            } else {
                if (lastType) {
                    sa[patLen - 1 - j] = i + 1;
                    ++j;
                }
                lastType = false;
            }
        }

        for (size_t i = 0; i < lmsCnt; ++i) {
            size_t r = sa[i];
            sa[i] = sa[patLen - lmsCnt + r];
            sa[patLen - lmsCnt + r] = EMPTY;
        }

        size_t tail = EMPTY, rfp = EMPTY;
        for (int64_t i = lmsCnt - 1; i >= 1; --i) {
            if (pat[sa[i]] != tail) {
                tail = pat[sa[i]];
                rfp = tail;
            }
            sa[rfp] = sa[i];
            if (rfp != (size_t)i) sa[i] = EMPTY;
            --rfp;
        }
    }

    void inducedSort(size_t *pat, size_t *sa, size_t patLen, size_t saLen){
        bool lastType = true;
        for (int64_t j = patLen - 2; j >= 0; --j) {
            if (patCharType(pat[j], pat[j + 1], lastType)) {
                lastType = true;
            } else {
                sa[pat[j]]++;
                lastType = false;
            }
        }

        int64_t i = 0;
        while (i < patLen) {
            if (sa[i] < EMPTY && sa[i] > 0) {
                size_t j = sa[i] - 1;
                bool isLType = false;

                if (pat[j] > pat[j + 1]) {
                    isLType = true;
                } else if (pat[j] == pat[j + 1]) {
                    size_t next_i = sa[pat[sa[i]]];
                    if (next_i >= MULTI) {
                        isLType = true;
                    } else if (next_i < EMPTY && pat[sa[i]] + 1 < patLen) {
                        size_t t = sa[pat[sa[i]] + 1];
                        if (t == EMPTY) isLType = true;
                        else if (t < EMPTY) {
                            if (pat[t] == pat[sa[i]]) isLType = true;
                        }
                    }
                }

                if (isLType) {
                    if (sa[pat[j]] == UNIQUE) {
                        sa[pat[j]] = j;
                    } else if (sa[pat[j]] >= MULTI && sa[pat[j] + 1] == EMPTY) {
                        if (sa[pat[j]] - EMPTY > 2) {
                            sa[pat[j] + 2] = j;
                            sa[pat[j] + 1] = 1;
                        } else {
                            sa[pat[j]] = j;
                        }
                    } else if (sa[pat[j]] >= MULTI && sa[pat[j] + 1] != EMPTY) {
                        size_t e = pat[j];
                        size_t c = sa[e + 1];
                        size_t lfp = e + c + 2;
                        if (c + 2 < sa[pat[j]] - EMPTY) {
                            sa[lfp] = j;
                            ++sa[e + 1];
                        } else {
                            for (size_t k = 1; k <= c; ++k) sa[e + k - 1] = sa[e + k + 1];
                            sa[e + c] = j;
                            sa[e + c + 1] = EMPTY;
                            if (i >= e + 2 && i <= e + c + 1) i -= 2;
                        }
                    } else if (sa[pat[j]] < EMPTY) {
                        for (size_t k = pat[j]; k < patLen; ++k) {
                            if (sa[k] == EMPTY) {
                                sa[k] = j;
                                break;
                            }
                        }
                    }
                }

            } else if (sa[i] >= MULTI) ++i;

            ++i;
        }

        lastType = true;
        for (int64_t j = patLen - 2; j >= 0; --j) {
            if (patCharType(pat[j], pat[j + 1], lastType)) {
                lastType = true;
            } else {
                if (lastType) {
                    if (sa[pat[j + 1]] <= EMPTY) {
                        sa[pat[j + 1]] = UNIQUE;
                    } else {
                        ++sa[pat[j + 1]];
                    }
                }
                lastType = false;
            }
        }
        i = patLen - 1;
        while (i > 0) {
            if (sa[i] > EMPTY) {
                size_t c = sa[i] - EMPTY;
                for (size_t k = 0; k < c; ++k) sa[i - k] = EMPTY;
                i -= c - 1;
            }
            i--;
        }
        sa[0] = patLen - 1;

        lastType = true;
        for (int64_t j = patLen - 2; j >= 0; --j) {
            if (patCharType(pat[j], pat[j + 1], lastType)) {
                if (sa[pat[j]] >= EMPTY) {
                    ++sa[pat[j]];
                } else {
                    sa[pat[j]] = UNIQUE;
                }
                lastType = true;
            } else {
                lastType = false;
            }
        }

        i = patLen - 1;
        while (i > 0) {
            if (sa[i] < EMPTY && sa[i] > 0) {
                size_t j = sa[i] - 1;
                bool isSType = false;
                if (pat[j] < pat[j + 1]) {
                    isSType = true;
                } else if (pat[j] == pat[j + 1]) {
                    size_t next_i = sa[pat[sa[i]]];
                    if (next_i >= MULTI) isSType = true;
                    else if (next_i < EMPTY && pat[sa[i]] - 1 > 0) {
                        size_t t = sa[pat[sa[i]] - 1];
                        if (t == EMPTY) isSType = true;
                        else if (t < EMPTY) {
                            if (pat[t] == pat[sa[i]]) isSType = true;
                        }
                    }
                }

                if (isSType) {
                    if (sa[pat[j]] == UNIQUE) sa[pat[j]] = j;
                    else if (sa[pat[j]] >= MULTI && sa[pat[j] - 1] == EMPTY) {
                        if (sa[pat[j]] - EMPTY > 2) {
                            sa[pat[j] - 2] = j;
                            sa[pat[j] - 1] = 1;
                        } else sa[pat[j]] = j;
                    } else if (sa[pat[j]] >= MULTI && sa[pat[j] - 1] != EMPTY) {
                        size_t e = pat[j];
                        size_t c = sa[e - 1];
                        size_t num = sa[pat[j]] - EMPTY;
                        if (c + 2 < num) {
                            size_t rfp = e - c - 2;
                            sa[rfp] = j;
                            ++sa[e - 1];
                        } else {
                            for (size_t k = 1; k <= c; ++k) sa[e - k + 1] = sa[e - k - 1];
                            sa[e - c] = j;
                            sa[e - c - 1] = EMPTY;
                            if (i >= e - num + 1 && i <= e - 2) i += 2;
                        }
                    } else if (sa[pat[j]] < EMPTY) {
                        for (size_t k = pat[j]; k > 0; --k) {
                            if (sa[k - 1] == EMPTY) {
                                sa[k - 1] = j;
                                break;
                            }
                        }
                    }
                }
            } else if (sa[i] >= MULTI) {
                --i;
            }
            --i;
        }
    }

    void computeSuffixArray(size_t *pat, size_t *sa, size_t patLen, size_t saLen){
        renamePat(pat, sa, patLen, saLen);
        size_t lmsCnt = sortLmsChar(pat, sa, patLen, saLen);
        sortLmsSubstr(pat, sa, patLen, saLen);
        bool dup = constructPat1(pat, sa, lmsCnt, patLen, saLen);
        sortLmsSuf(pat, sa, patLen, saLen, lmsCnt, dup);
        inducedSort(pat, sa, patLen, saLen);
    }




    void SuffixArray::buildSAIS(const std::string &seq) {
        //todo: optimize memory allocations: text_ is useless
        std::vector<unsigned char> text_;
        const int64_t n = seq.length();
        text_.resize(n , -1);
        for (size_t i = 0; i < n;++i){
            int c = charToIndex(seq[i]);
            text_[i] = c < 0 ? 5 : c + 1; // ensure that invalid characters are at the end
        }
        size_t patLen = n + 1;
        size_t* p = new size_t[patLen];
        for (size_t i = 0; i < n; ++i) {
            p[i] = text_[i];
        }
        p[n] = 0;
        size_t saLen = std::max(n, (int64_t) 256);
        auto* rawSA = new size_t[saLen + reservedLength];
        memset(rawSA,0, (saLen + reservedLength) * sizeof(size_t));
        size_t* sa = rawSA + reservedLength;
        computeSuffixArray(p, sa, patLen, saLen);
        delete[] p;

        size_t maxValue = n - 1;
        uint32_t bitsNeeded = 0;
        while (maxValue > 0) {
            bitsNeeded++;
            maxValue >>= 1;
        }
        suffixArray_.buildFromReservedArray(rawSA, bitsNeeded, n, reservedLength);
        size_t validCount = 0;
        for (size_t i = 1; i <= n; ++i) {
            if (text_[sa[i]] != 5) {
                suffixArray_.setValue(validCount, sa[i]);
                validCount++;
            }
        }
        suffixArray_.setLength(validCount);
        length_ = validCount;
        wordBits_ = bitsNeeded;
    }
}
