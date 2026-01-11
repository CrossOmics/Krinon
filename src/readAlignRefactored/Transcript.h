#ifndef RNAALIGNREFACTORED_TRANSCRIPT_H
#define RNAALIGNREFACTORED_TRANSCRIPT_H
#include <string>
#include <vector>
namespace RefactorProcessing {
    struct Exon {

    };

    struct SpliceJunction{

    };

    struct WindowAlign{
        int64_t readStart{0};
        int64_t genomeStart{0};
        int64_t length{0};
        int64_t score{0};
        bool isAnchor{false};
        int64_t isj{-1};// annotation index
        int iFragment{0};

        bool operator<(const WindowAlign &other) const {
            if (genomeStart != other.genomeStart) return genomeStart < other.genomeStart;
            return length < other.length;
        }
    };

    struct Transcript {
        std::string readName;
        std::string chr;
        std::string CIGAR;
        int strand{0};
        int64_t matched{0};
        int64_t unmatched{0};
        int64_t nIns{0};
        int64_t nDel{0};
        std::vector<Exon> exons;
        std::vector<SpliceJunction> sj;
        std::vector<WindowAlign> aligns;//Actually no need to store this, for debugging only
        int64_t readStart{0};
        int64_t genomeStart{0};
        int64_t posInChr{0};
        int64_t score{0};
        int64_t readLength{0};


    };
}
#endif //RNAALIGNREFACTORED_TRANSCRIPT_H
