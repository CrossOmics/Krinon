#include "SuffixArray.h"
#include "../utilsRefactored/seqFunctions.h"
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <immintrin.h>
#include "omp.h"
namespace RefactorProcessing {
    constexpr uint64_t DEFAULT_SUBPROBLEM_COUNT = 8192;
    constexpr uint64_t NESTED_PAR_GRAIN = 1ul << 13;  // 8192,1ul<<13


    inline uint64_t lcp_avx(const char* a, const char* b, uint64_t n) {
        uint64_t i = 0;
        while (i + 32 <= n) {
            __m256i va = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a + i));
            __m256i vb = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(b + i));
            uint32_t mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(va, vb));
            if (mask != 0xFFFFFFFF) {
                i += __builtin_ctz(~mask);
                return i;
            }
            i += 32;
        }
        while (i < n && a[i] == b[i]) ++i;
        return i;
    }


    void merge(const uint64_t* X, uint64_t len_x,
               const uint64_t* Y, uint64_t len_y,
               const uint64_t* LCP_x, const uint64_t* LCP_y,
               uint64_t* Z, uint64_t* LCP_z,
               const char* T, uint64_t total_n, uint64_t max_context) {
        uint64_t m = 0;          // LCP of last compared pair
        uint64_t i = 0, j = 0, k = 0;

        while (i < len_x && j < len_y) {
            uint64_t lx = LCP_x[i];

            if (lx > m) {
                Z[k] = X[i];
                LCP_z[k] = lx;

            } else if (lx < m) {
                Z[k] = Y[j];
                LCP_z[k] = m;
                m = lx;

            } else {  // lx == m
                uint64_t max_n = total_n - std::max(X[i], Y[j]); // shorter suffix length
                uint64_t context = std::min(max_context, max_n);
                uint64_t ncommon = m + lcp_avx(T + (X[i] + m), T + (Y[j] + m), context - m);

                Z[k] = (ncommon == max_n ? std::max(X[i],Y[j]):(T[X[i] + ncommon] < T[Y[j] +ncommon]? X[i]:Y[j]));

                LCP_z[k] = (Z[k] == X[i] ? lx : m);
                m = ncommon;
            }
            if (Z[k] == X[i]) i++;
            else {
                j++;
                std::swap(X, Y);
                std::swap(len_x, len_y);
                std::swap(LCP_x, LCP_y);
                std::swap(i, j);
            }
            k++;


        }

        std::memcpy(Z + k, X + i, (len_x - i) * sizeof(uint64_t));
        std::memcpy(LCP_z + k, LCP_x + i, (len_x - i) * sizeof(uint64_t));
        std::memcpy(Z + k, Y + j, (len_y - j) * sizeof(uint64_t));
        std::memcpy(LCP_z + k, LCP_y + j, (len_y - j) * sizeof(uint64_t));

        if (k < len_x + len_y) LCP_z[k] = m;


    }


    void merge_sort(uint64_t* X, uint64_t* Y, uint64_t n,
                    uint64_t* LCP, uint64_t* W,
                    const char* T, uint64_t total_n, uint64_t max_context) {
        assert(std::memcmp(X, Y, n * sizeof(uint64_t)) == 0);

        if (n == 1) {
            LCP[0] = 0;
            return;
        }

        uint64_t m = n / 2;

        if (m < NESTED_PAR_GRAIN || n - m < NESTED_PAR_GRAIN) {
            // 串行递归
            merge_sort(Y, X, m, W, LCP, T, total_n, max_context);
            merge_sort(Y + m, X + m, n - m, W + m, LCP + m, T, total_n, max_context);
        } else {
#pragma omp task
            merge_sort(Y, X, m, W, LCP, T, total_n, max_context);
#pragma omp task
            merge_sort(Y + m, X + m, n - m, W + m, LCP + m, T, total_n, max_context);
#pragma omp taskwait
        }

        merge(X, m, X + m, n - m, W, W + m, Y, LCP, T, total_n, max_context);
    }


    void sample_pivots(const uint64_t* X, uint64_t n, uint64_t m, uint64_t* P) {
        assert(m <= n);
        uint64_t gap = n / m;
        for (uint64_t i = 0; i < m; ++i)
            P[i] = X[(i + 1) * gap - 1];
    }


    uint64_t upper_bound(const uint64_t* X, uint64_t n, const char* P, uint64_t P_len,
                         const char* T, uint64_t total_n, uint64_t max_context) {
        int64_t l = -1, r = n;
        uint64_t soln = n;
        uint64_t lcp_l = 0, lcp_r = 0;
        constexpr uint64_t cutoff = 65536;

        while (r - l > 1) {
            uint64_t c = (l + r) / 2;
            const char* suf = T + X[c];
            uint64_t suf_len = total_n - X[c];

            uint64_t lcp_c = std::min(lcp_l, lcp_r);
            lcp_c = std::min(lcp_c, cutoff);
            uint64_t max_lcp = std::min({suf_len, P_len, max_context, cutoff});
            lcp_c += lcp_avx(suf + lcp_c, P + lcp_c, max_lcp - lcp_c);

            if (lcp_c == max_lcp) {
                if (lcp_c == P_len) {
                    if (P_len == suf_len)
                        return c + 1;
                    else
                        r = c, lcp_r = lcp_c, soln = c;
                } else {
                    l = c, lcp_l = lcp_c;
                }
            } else {
                if (suf[lcp_c] < P[lcp_c]) {
                    l = c, lcp_l = lcp_c;
                } else {
                    r = c, lcp_r = lcp_c, soln = c;
                }
            }
        }
        return soln;
    }


    void sort_partition(uint64_t* X, uint64_t* Y, uint64_t n, const uint64_t* S,
                        uint64_t* LCP_x, uint64_t* LCP_y,
                        const char* T, uint64_t total_n, uint64_t max_context) {
        if (n == 1) return;

        uint64_t m = n / 2;
        uint64_t left_size = S[m] - S[0];
        uint64_t right_size = S[n] - S[m];

        if (left_size < NESTED_PAR_GRAIN || right_size < NESTED_PAR_GRAIN) {
            sort_partition(Y, X, m, S, LCP_y, LCP_x, T, total_n, max_context);
            sort_partition(Y + left_size, X + left_size, n - m, S + m,
                           LCP_y + left_size, LCP_x + left_size, T, total_n, max_context);
        } else {
#pragma omp task
            sort_partition(Y, X, m, S, LCP_y, LCP_x, T, total_n, max_context);
#pragma omp task
            sort_partition(Y + left_size, X + left_size, n - m, S + m,
                           LCP_y + left_size, LCP_x + left_size, T, total_n, max_context);
#pragma omp taskwait
        }

        merge(X, left_size, X + left_size, right_size,
              LCP_x, LCP_x + left_size, Y, LCP_y,
              T, total_n, max_context);
    }




    void caps_sa(const char* T, uint64_t n,
                     uint64_t subproblem_count, uint64_t max_context,
                     uint64_t* SA, uint64_t* LCP) {

        if (max_context == 0) max_context = n;
        uint64_t p = std::min(subproblem_count > 0 ? subproblem_count : DEFAULT_SUBPROBLEM_COUNT, n / 16);
        if (p > n) {
            std::cerr << "Error: subproblem count too large.\n";
            std::exit(EXIT_FAILURE);
        }
        if (n < 16) {
            std::cerr << "Error: input too short.\n";
            std::exit(EXIT_FAILURE);
        }

        uint64_t pivot_per_part = std::min(static_cast<uint64_t>(std::ceil(32.0 * std::log(n))), n / p - 1);
        uint64_t subarr_size = n / p;
        uint64_t last_subarr_size = subarr_size + n % p;


        uint64_t* SA_w = static_cast<uint64_t*>(std::malloc(n * sizeof(uint64_t)));
        uint64_t* LCP_w = static_cast<uint64_t*>(std::malloc(n * sizeof(uint64_t)));
        uint64_t* pivot = static_cast<uint64_t*>(std::malloc(p * pivot_per_part * sizeof(uint64_t)));
        uint64_t* part_size_scan = static_cast<uint64_t*>(std::malloc((p + 1) * sizeof(uint64_t)));
        uint64_t* part_ruler = static_cast<uint64_t*>(std::malloc(p * (p + 1) * sizeof(uint64_t)));

        if (!SA_w || !LCP_w || !pivot || !part_size_scan || !part_ruler) {
            std::cerr << "Memory allocation failed.\n";
            std::exit(EXIT_FAILURE);
        }


#pragma omp parallel for
        for (uint64_t i = 0; i < n; ++i) {
            SA[i] = i;
            SA_w[i] = i;
        }


#pragma omp parallel
        {
#pragma omp single
            {
                for (uint64_t i = 0; i < p; ++i) {
                    uint64_t start = i * subarr_size;
                    uint64_t len = (i == p-1) ? last_subarr_size : subarr_size;
#pragma omp task
                    merge_sort(SA_w + start, SA + start, len,
                               LCP + start, LCP_w + start,
                               T, n, max_context);
                }
            }
        }


#pragma omp parallel for
        for (uint64_t i = 0; i < p; ++i) {
            uint64_t start = i * subarr_size;
            uint64_t len = (i == p-1) ? last_subarr_size : subarr_size;
            sample_pivots(SA + start, len, pivot_per_part, pivot + i * pivot_per_part);
        }


        uint64_t sample_count = p * pivot_per_part;
        uint64_t* pivot_w = static_cast<uint64_t*>(std::malloc(sample_count * sizeof(uint64_t)));
        uint64_t* temp1 = static_cast<uint64_t*>(std::malloc(sample_count * sizeof(uint64_t)));
        uint64_t* temp2 = static_cast<uint64_t*>(std::malloc(sample_count * sizeof(uint64_t)));
        if (!pivot_w || !temp1 || !temp2) {
            std::cerr << "Memory allocation failed.\n";
            std::exit(EXIT_FAILURE);
        }
        std::memcpy(pivot_w, pivot, sample_count * sizeof(uint64_t));

#pragma omp parallel
        {
#pragma omp single
            merge_sort(pivot, pivot_w, sample_count, temp1, temp2, T, n, max_context);
        }


        sample_pivots(pivot_w, sample_count, p - 1, pivot);
        std::free(pivot_w);
        std::free(temp1);
        std::free(temp2);


        uint64_t* P = static_cast<uint64_t*>(std::malloc(p * (p + 1) * sizeof(uint64_t)));
        if (!P) {
            std::cerr << "Memory allocation failed.\n";
            std::exit(EXIT_FAILURE);
        }

#pragma omp parallel for
        for (uint64_t i = 0; i < p; ++i) {
            uint64_t start = i * subarr_size;
            uint64_t len = (i == p-1) ? last_subarr_size : subarr_size;
            uint64_t* Pi = P + i * (p + 1);
            Pi[0] = 0;
            Pi[p] = len;
            for (uint64_t j = 0; j < p - 1; ++j) {
                Pi[j + 1] = upper_bound(SA + start, len, T + pivot[j], n - pivot[j],
                                        T, n, max_context);
            }
        }


#pragma omp parallel for
        for (uint64_t j = 0; j < p; ++j) {
            uint64_t total = 0;
            for (uint64_t i = 0; i < p; ++i) {
                const uint64_t* Pi = P + i * (p + 1);
                total += (Pi[j + 1] - Pi[j]);
            }
            part_size_scan[j] = total;
        }


        uint64_t cur = 0;
        for (uint64_t j = 0; j < p; ++j) {
            uint64_t sz = part_size_scan[j];
            part_size_scan[j] = cur;
            cur += sz;
        }
        part_size_scan[p] = cur;
        assert(part_size_scan[p] == n);


#pragma omp parallel for
        for (uint64_t j = 0; j < p; ++j) {
            uint64_t dest_start = part_size_scan[j];
            uint64_t* dest_SA = SA_w + dest_start;
            uint64_t* dest_LCP = LCP_w + dest_start;
            uint64_t* ruler = part_ruler + j * (p + 1);
            uint64_t offset = 0;
            for (uint64_t i = 0; i < p; ++i) {
                const uint64_t* Pi = P + i * (p + 1);
                uint64_t src_start = i * subarr_size + Pi[j];
                uint64_t len = Pi[j + 1] - Pi[j];
                ruler[i] = offset;
                if (len > 0) {
                    std::memcpy(dest_SA + offset, SA + src_start, len * sizeof(uint64_t));
                    std::memcpy(dest_LCP + offset, LCP + src_start, len * sizeof(uint64_t));
                    dest_LCP[ruler[i]] = 0;
                    offset += len;
                }
            }
            ruler[p] = offset;
            assert(offset == part_size_scan[j + 1] - part_size_scan[j]);
        }

        std::free(P);


#pragma omp parallel for
        for (uint64_t j = 0; j < p; ++j) {
            uint64_t start = part_size_scan[j];
            uint64_t len = part_size_scan[j + 1] - start;
            std::memcpy(SA + start, SA_w + start, len * sizeof(uint64_t));
            std::memcpy(LCP + start, LCP_w + start, len * sizeof(uint64_t));
        }


#pragma omp parallel
        {
#pragma omp single
            {
                for (uint64_t j = 0; j < p; ++j) {
                    uint64_t start = part_size_scan[j];
                    uint64_t len = part_size_scan[j + 1] - start;
                    const uint64_t* ruler = part_ruler + j * (p + 1);
#pragma omp task
                    sort_partition(SA_w + start, SA + start, p, ruler,
                                   LCP_w + start, LCP + start,
                                   T, n, max_context);
                }
            }
        }


#pragma omp parallel for
        for (uint64_t j = 1; j < p; ++j) {
            uint64_t idx = part_size_scan[j];
            LCP[idx] = lcp_avx(T + SA[idx - 1], T + SA[idx],
                               n - std::max(SA[idx - 1], SA[idx]));
        }


        std::free(SA_w);
        std::free(LCP_w);
        std::free(pivot);
        std::free(part_size_scan);
        std::free(part_ruler);

        std::cerr << "Suffix array construction completed.\n";
    }

    void SuffixArray::buildCAPS(const std::string &seq) {
        if (threadNum > 1) omp_set_num_threads(threadNum);
        const size_t n = seq.length();

        size_t patLen = n + 1;
        char *p = new char[patLen];

        for (size_t i = 0; i < n; ++i) {
            int c = charToIndex(seq[i]);
            p[i] = c < 0 ? 5 : c + 1;
        }
        p[n] = 0;

        size_t saLen = std::max(n, (size_t) 256);
        auto* rawSA = new size_t[saLen + reservedLength];
        memset(rawSA,0, (saLen + reservedLength) * sizeof(size_t));
        size_t* sa = rawSA + reservedLength;
        size_t* rawLCP = new size_t [patLen];

        caps_sa(p,n,DEFAULT_SUBPROBLEM_COUNT,0,sa,rawLCP);

        delete[] p;
        delete[] rawLCP;

        size_t maxValue = n - 1;
        uint32_t bitsNeeded = 0;
        while (maxValue > 0) {
            bitsNeeded++;
            maxValue >>= 1;
        }
        suffixArray_.buildFromReservedArray(rawSA, bitsNeeded, n, reservedLength);
        size_t validCount = 0;
        for (size_t i = 0; i < n; ++i) {
            if (charToIndex(seq[sa[i]]) >= 0) {
                suffixArray_.setValue(validCount, sa[i]);
                validCount++;
            }
        }
        suffixArray_.setLength(validCount);
        length_ = validCount;
        wordBits_ = bitsNeeded;

    }
}