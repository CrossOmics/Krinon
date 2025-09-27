#ifndef RNAALIGNREFACTORED_WINDOWMANAGEMENT_H
#define RNAALIGNREFACTORED_WINDOWMANAGEMENT_H

#include "../utils/types.h"
#include "../genome/GenomeIndex.h"
#include <vector>
#include <unordered_map>

namespace rna {
    struct StitchingConfig {
        int maxAnchorRep{50};
        int winBinSizeLog{16};                      // the size of each window is k* 2^winBinSizeLog
        int winAnchorDistBins{9};
        int flankSize{4};
        int maxWindows{10000};
        int maxSeedPerWindows{50};
        //size_t maxTranscriptsPerWindow{20};//only use in recursive stitching

        //int minAlignLength{10};
        int maxRep{10000};
        int maxExons{20};                           // max exons number in a transcript
        int transcriptStoredMax{100};               // maximum number of transcripts to store

        int outFilterMultimapMax{10};
        int maxMismatch{10};
        int multimapScoreRange{1};
        double outFilterScoreMinOverLRead{0.66};
        double outFilterMatchMinOverLRead{0.66};
    };

    struct StitchingScoreConfig {
        int MATCH_SCORE = 1;
        int MISMATCH_PENALTY = -1;
        int GAP_OPEN_PENALTY = 0;
        int DEL_OPEN_PENALTY = -2;
        int DEL_EXTEND_PENALTY = -2;
        int INS_OPEN_PENALTY = -2;
        int INS_EXTEND_PENALTY = -2;
        int SCORE_STITCH_SJ_SHIFT = 1;              // maximum score reduction while searching for SJ boundaries
        int SCORE_GAP_GCAG = -4;                    // GC/AG and CT/GC junction penalty
        int SCORE_GAP_ATAC = -8;                    // AT/AC and GT/AT junction penalty
        int SCORE_GAP_NON_CANONICAL = -8;
        int SCORE_ANNOTATED_SJ = 2;
        int MAX_SJ_REPEAT_SEARCH = 255;
        int MIN_INTRON_LENGTH = 21;
        int MAX_INTRON_LENGTH = 2147483647;
        int MAX_MISMATCH_FOR_SJ[4] = {0, -1, 0, 0}; // max mismatch allowed for different SJ types: non-canonical, 
                                                    // GT-AG, GC-AG, AT-AC, -1 means no limit
    };

    class StitchingManagement {
    public:
        explicit StitchingManagement(const StitchingConfig &config, const GenomeIndex &);

        ~StitchingManagement();

        void processAlignments(Align *alignments, int alignNum, ReadPtr read);

        TranscriptPtr getBestTranscript() const;

        std::vector<SJDBOutput> getSJDB() const;

        void clear();

        void setConfig(const StitchingConfig &config,const StitchingScoreConfig &scoreConfig){
            config_ = config;
            scoreConfig_ = scoreConfig;
        }

        enum StitchingStatus {
            SUCCESS,
            FAILED_NO_ALIGNMENTS,
            FAILED_NO_GOOD_TRANSCRIPT,
            FAILED_TOO_MANY_TRANSCRIPTS,
        } status = SUCCESS;

        Transcript *goodTranscripts_;
        int numGoodTranscripts_ = 0;
    private:
        struct PositiveStrandAlign {

            int64_t genomeStart;
            int64_t readStart;
            int64_t length;
            int direction;
            bool isAnchor;
            int iFragment{0};
        };

        struct sjSplitRecord{
            int64_t isj;
            int64_t donorLength;
            int64_t acceptorLength;
            int64_t donorStart;
            int64_t acceptorStart;
            int direction;
        };

        //Windows
        inline PositiveStrandAlign convertAlignToPositiveStrand(const Align &align, size_t SAi);

        void createWindows(const Align *alignments, int alignNum);

        void createRawWindowsSingle(const PositiveStrandAlign &a);

        std::pair<bool,sjSplitRecord> trySplitSJ(const Align &align, int64_t location);

        void assignAlignmentsToWindows(const Align *alignments, int alignNum);

        inline void assignSingleAlignment(Window &win, const PositiveStrandAlign &a,int64_t isj = -1) const;


        void identifyAnchors(Align *alignments, int alignNum);


        //DP stitching
        void stitchWindowsAlignNew(Window &window);

        void generateTranscriptsNew();

        void stitchingBetweenWindowAligns(const WindowAlign &a1, const WindowAlign &a2, int windowDir,
                                          StitchingRecord &record);

        void extendWindowAlign(const WindowAlign &a, int windowDir, int extendDir, ExtensionRecord &res);

        void refreshWinBinMap();




        static inline std::pair<int, int64_t>
        checkJunctionMotif(const std::string &genomeSeq, int64_t leftPos, int64_t rightPos);


        StitchingConfig config_;
        StitchingScoreConfig scoreConfig_;
        std::vector<Window> windows_;
        std::vector<Transcript> transcripts_;
        std::vector<SJDBOutput> sjdb;
        int64_t outFilterScoreMin_;
        int outFilterMatchMin_;

        WindowAlign *allWindowAligns_;              // window alignments stored for all windows.
        int32_t *winBinMap_[2];

        StitchingRecord *nowStitchingRecord_;       // stitching records for current window
        ExtensionRecord *nowExtensionRecord_[2];    // extension records for current window, forward & backward
        RawTranscript *nowRawTranscript_;           // raw transcripts for current window, used to store stitching records
        ExtensionRecord::singleExtensionRecord *allSingleExtensionRecord_;

        int64_t maxTranscriptScore_;


        const GenomeIndex &genomeIndex_;
        ReadPtr read_;
        int64_t readLength{0};
        TranscriptPtr bestTranscript_;


        uint32_t totalAnchors_ = 0;
        // todo move to config
        static constexpr int MATCH_SCORE = 1;
        static constexpr int MISMATCH_PENALTY = -1;
        static constexpr int GAP_OPEN_PENALTY = 0;
        static constexpr int DEL_OPEN_PENALTY = -2;
        static constexpr int DEL_EXTEND_PENALTY = -2;
        static constexpr int INS_OPEN_PENALTY = -2;
        static constexpr int INS_EXTEND_PENALTY = -2;
        static constexpr int SCORE_STITCH_SJ_SHIFT = 1;// maximum score reduction while searching for SJ boundaries
        static constexpr int SCORE_GAP_GCAG = -4; // GC/AG and CT/GC junction penalty
        static constexpr int SCORE_GAP_ATAC = -8; // AT/AC and GT/AT junction penalty
        static constexpr int SCORE_GAP_NON_CANONICAL = -8;
        static constexpr int SCORE_ANNOTATED_SJ = 2;
        static constexpr int MAX_SJ_REPEAT_SEARCH = 255;
        static constexpr uint32_t MIN_INTRON_LENGTH = 20;
        static constexpr uint32_t MAX_INTRON_LENGTH = 2147483647;
        static constexpr int MAX_MISMATCH_FOR_SJ[4] = {0, -1, 0,0}; // max mismatch allowed for different SJ types: non-canonical, GT-AG, GC-AG, AT-AC, -1 means no limit
        static constexpr int ALIGN_ENDS_PROTRUDE = 0;

        static constexpr inline int32_t charToIndex(char c) {
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
    };
} // namespace rna
#endif //RNAALIGNREFACTORED_WINDOWMANAGEMENT_H
