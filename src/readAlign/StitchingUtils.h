#ifndef RNAALIGNREFACTORED_STITCHINGUTILS_H
#define RNAALIGNREFACTORED_STITCHINGUTILS_H

#include <cstdint>
namespace rna{
    struct SpliceJunction {

        int64_t type{0};//canon
        int64_t length{0};
        bool isAnnotated{false};
        int64_t shiftLeft{0};
        int64_t shiftRight{0};
        int64_t sjStrand{0};
    };


// the change of the score, matches, mismatches, exons length of one single stitching step (stitching WA[j] to WA[i])
    struct StitchingRecord {
        enum Type {
            CANNOT_STITCH,
            PERFECT_MATCH,
            SAME_GAP,
            DELETION,
            INSERTION,
            SPLICE_JUNCTION,
            CROSS_FRAGMENTS
        } type{CANNOT_STITCH};
        int64_t score{0}; // including the gap penalty
        int matches{0};
        int mismatches{0};
        int64_t formerExonLengthShift{0};
        int64_t latterExonLengthShift{0};

        SpliceJunction stitchingType{0, 0, false};
    };

    struct ExtensionRecord {
        //score can be calculated by matches * MATCH_SCORE - i * MATCH_SCORE
        struct singleExtensionRecord {
            int length{0}; // extension length
            int matched{0}; // number of matched bases, in case of the appearance of 'N'
            int mismatches{
                    0}; // number of mismatches, since the best extension may not reach the max mismatch tolerance
        };
        // array of length [max_mismatch], to store the max extension length when allowing a certain number of mismatches
        // Note: must free the memory after use
        singleExtensionRecord *maxExtensionLengthWithMismatch{nullptr};
        int maxMismatch{0};
        int maxExtensionScore{0};

    };



// selected pieces and corresponding stitching records
// can be used to reconstruct the transcript
// must ensure that the score is right (same as the transcript's final score)
    struct RawTranscript {
        int previousTranscriptId{-1};//id in the allWindowAligns_ array
        int newAlignId{-1}; //
        int64_t score{0}; // true score
        int mismatches{0}; // number of mismatches in the transcript
        int matches{0};
        int exonCount{0}; // number of exons in the transcript
        int extendedLengthForward{0}; //extension 5' end
        int extendedLengthBackward{0}; // extension 3' end
        int startAlignId{-1};
        bool crossFragments{false};
        int extendedLengthFragmentFormer{0};
        int extendedLengthFragmentLatter{0};
        int64_t firstFragmentMatchEnd{0};
    };

    struct FragmentStitchingRecord {
        int64_t score{0};
        int mismatches{0};
        int matches{0};
        RawTranscript* firstFragmentRecord{nullptr};
        RawTranscript* secondFragmentRecord{nullptr};
    };

}

#endif //RNAALIGNREFACTORED_STITCHINGUTILS_H
