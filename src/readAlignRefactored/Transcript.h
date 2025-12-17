#ifndef RNAALIGNREFACTORED_TRANSCRIPT_H
#define RNAALIGNREFACTORED_TRANSCRIPT_H
#include <string>
#include <vector>
namespace RefactorProcessing {
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
