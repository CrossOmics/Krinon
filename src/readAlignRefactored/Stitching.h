#ifndef RNAALIGNREFACTORED_STITCHINGR_H
#define RNAALIGNREFACTORED_STITCHINGR_H
//todo
#include <vector>
#include "../genomeRefactored/GenomeIndex.h"
#include "../io/Parameters.h"
#include "Transcript.h"
#include "../utilsRefactored/seqFunctions.h"

namespace RefactorProcessing {
    struct WindowAlign{
        int64_t readStart{0};
        int64_t genomeStart{0};
        int64_t length{0};
        bool isAnchor{false};
        int64_t isj{-1};// annotation index
        int iFragment{0};

        bool operator<(const WindowAlign &other) const {
                if (genomeStart != other.genomeStart) return genomeStart < other.genomeStart;
                return length < other.length;
        }
    };

    struct StitchRecord{
        //pre-calculated stitching record between two window-alignments
        enum Type {
            CANNOT_STITCH,
            PERFECT_MATCH,
            SAME_GAP,
            DELETION,
            INSERTION,
            SPLICE_JUNCTION,
            CROSS_FRAGMENTS
        } type{CANNOT_STITCH};
        int spliceJunctionType{0}; //todo 0: none, 1: GT-AG, 2: GC-AG, 3: AT-AC, 4: non-canonical, to be merged with type?
        int64_t score{0}; // including the gap penalty
        int matches{0}; // ADDITIONAL matches, matches contributed by window aligns not included
        int mismatches{0}; //ADDITIONAL mismatches. same as above
        int64_t formerExonLengthShift{0}; // length extension of the former exon
        int64_t latterExonLengthShift{0}; // length extension of the latter exon

    };

    struct ExtendRecord{
        struct singleExtendRecord{
            int length{0}; // extension length
            int matched{0}; // number of matched bases, in case of the appearance of 'N'
            int mismatches{0}; // number of mismatches, since the best extension may not reach the max mismatch tolerance
        };
        // array of length [max_mismatch], to store the max extension length when allowing a certain number of mismatches
        singleExtendRecord *maxExtensionLengthWithMismatch{nullptr};
        int maxMismatch{0};
        int maxExtensionScore{0};
    };



    class Window {

    };

    class RawTranscript {

    public:
        Transcript convertToTranscript();
    };

    class Stitching {
        // stitch alignment fragments
    private:
        // configs
        bool isPaired_;
        int maxAnchorRep_;
        int winBinSizeLog_;
        int winAnchorDistBins_;
        int flankSize_;
        int maxWindows_;
        int maxSeedPerWindows_;
        int maxRep_;
        int maxExons_;
        int transcriptStoredMax_;

        //output filter
        int outFilterMultimapMax_;
        int maxMismatch_;
        int multimapScoreRange_;
        double outFilterScoreMinOverLRead_;
        double outFilterMatchMinOverLRead_;


        // scoring
        int MATCH_SCORE_;
        int MISMATCH_PENALTY_;
        int GAP_OPEN_PENALTY_;
        int DEL_OPEN_PENALTY_;
        int DEL_EXTEND_PENALTY_;
        int INS_OPEN_PENALTY_;
        int INS_EXTEND_PENALTY_;
        int SCORE_STITCH_SJ_SHIFT_;
        int SCORE_GAP_GCAG_;
        int SCORE_GAP_ATAC_;
        int SCORE_GAP_NON_CANONICAL_;
        int SCORE_ANNOTATED_SJ_;
        int MAX_SJ_REPEAT_SEARCH_;
        int MIN_INTRON_LENGTH_;
        int MAX_INTRON_LENGTH_;
        int MAX_MISMATCH_FOR_SJ_[4];// max mismatch allowed for different SJ types




    public:

        // data
        std::vector<Align> alignments_;//  alignments to be stitched

        std::vector<WindowAlign> windowAlignments_;

        std::vector<Window> windows_;

        std::vector<int32_t> winBinMap_[2];

        std::vector<StitchRecord> stitchRecords_;// stitching records between window-alignments

        std::vector<ExtendRecord> extendRecords_[2];// extension records, forward & backward

        Stitching(){};
        ~Stitching(){};

        // parameters
        void setParam(const Parameters &P);

        // stitching process
        void init();

        void stitchAlignments();

        void createWindows(const std::vector<Align>& alignments);

        void assignAlignmentsToWindows();

        void convertAlignment();// convert alignment to window-alignment

        void assignAlignment();

        void stitchWindow_SingleEnd();

        void stitchWindow_PairedEnd();

    };
};

#endif //RNAALIGNREFACTORED_STITCHINGR_H
