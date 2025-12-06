#ifndef RNAALIGNMENT_TYPES_H
#define RNAALIGNMENT_TYPES_H

#include <cstdint>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>



#define READ_PAIRED 0x1
#define READ_MAPPED_IN_PROPER_PAIR 0x2
#define READ_UNMAPPED 0x4
#define READ_MATE_UNMAPPED 0x8
#define READ_REVERSE_STRAND 0x10
#define READ_MATE_REVERSE_STRAND 0x20
#define READ_FIRST_IN_PAIR 0x40
#define READ_SECOND_IN_PAIR 0x80

namespace rna {
    using GenomePos = int64_t;
    using ReadPos = int64_t;
    using Score = int64_t;
    using Length = int64_t;
    struct Read {
        std::string sequence[2]; // 0 for forward, 1 for complementary
        std::string quality;
        std::string name;
        int64_t length;
        int64_t mate1Length; // for paired-end reads
        int64_t mate2Length;
    };

    struct MatchResultStored {
        uint8_t length;
        int32_t upperRange;
        int64_t leftSAIndex;

        MatchResultStored() : length(0), upperRange(-1), leftSAIndex(-1) {};
    };

    struct MatchStats {
        int64_t maxLength;
        size_t count;
        size_t LongestMatchLeftPosInSA;
        size_t LongestMatchRightPosInSA;

        MatchStats(int64_t len = 0, size_t cnt = 0, size_t l = 0, size_t r = 0) : maxLength(len), count(cnt),
                                                                                  LongestMatchLeftPosInSA(l),
                                                                                  LongestMatchRightPosInSA(r) {}
    };

    struct Split {
        ReadPos readStart{0};
        size_t length{0};
        int iFragment{0};
    };

    struct Align {
        ReadPos readStart{0};
        int64_t length{0};
        size_t leftSAIndex{0};
        size_t rightSAIndex{0};
        size_t rep{0};
        int direction{0};
        bool isAnchor{false};
        int iFragment{0};

        bool operator<(const Align &other) const {
            return readStart < other.readStart;
        }
    };











    using ReadPtr = std::shared_ptr<Read>;


}
#endif //RNAALIGNMENT_TYPES_H
