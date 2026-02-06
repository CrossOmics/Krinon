#ifndef RNAALIGNREFACTORED_TRANSCRIPT_H
#define RNAALIGNREFACTORED_TRANSCRIPT_H
#include <string>
#include <vector>
namespace RefactorProcessing {
    struct Exon {
        int64_t genomeStart;
        int64_t length;
        int64_t readStart;
    };

    struct SpliceJunction{
        int type{0};//canon
        int64_t length{0};
        bool isAnnotated{false};
        int64_t shiftLeft{0};
        int64_t shiftRight{0};
        int64_t sjStrand{0};
    };



    struct Transcript {
        std::string readName;
        std::string chr;
        std::string CIGAR;
        int strand{0};
        int64_t matched{0};
        int64_t mismatched{0};
        int64_t nIns{0};
        int64_t nDel{0};
        std::vector<Exon> exons;
        std::vector<SpliceJunction> sj;
        //std::vector<WindowAlign> aligns;//Actually no need to store this, for debugging only
        int64_t readStart{0};
        int64_t genomeStart{0};
        int64_t posInChr{0};
        int64_t score{0};
        int64_t readLength{0};

        std::string getCIGAR() const;


    };
}
#endif //RNAALIGNREFACTORED_TRANSCRIPT_H
