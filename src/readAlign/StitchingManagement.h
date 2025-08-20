#ifndef RNAALIGNREFACTORED_WINDOWMANAGEMENT_H
#define RNAALIGNREFACTORED_WINDOWMANAGEMENT_H
#include "../utils/types.h"
#include "../genome/GenomeIndex.h"
#include "../genome/GenomeIndexPrefix.h"
#include <vector>
#include <unordered_map>
namespace rna {
    struct stitchingConfig {
        size_t maxAnchorRep{50};
        size_t winBinSizeLog{16}; //the size of each window is k* 2^winBinSizeLog
        int winAnchorDistBins{9};
        int maxWindows{10000};
        int maxSeedPerWindows{50};
        size_t maxTranscriptsPerWindow{20};//only use in recursive stitching
        int flankSize{4};
        int minAlignLength{10};
        int maxRep{10000};
        int maxMismatch{10};
        int multimapScoreRange{1};
        int transcriptStoredMax{100}; // maximum number of transcripts to store
        int maxExons{20}; // max exons number in a transcript
        double outFilterScoreMinOverLRead{0.66};
        double outFilterMatchMinOverLRead{0.66};
    };
    class StitchingManagement {
    public:
        explicit StitchingManagement(const stitchingConfig& config,const GenomeIndexPrefix&);
        ~StitchingManagement();
        void processAlignments(std::vector<Align>& alignments,ReadPtr read);
        TranscriptPtr getBestTranscript() const;
        std::vector<SJDBOutput> getSJDB() const;
        void clear();


    private:
        struct PositiveStrandAlign{

            int64_t genomeStart;
            int64_t readStart;
            int64_t length;
            int direction;
            bool isAnchor;
        };
        //Windows
        inline PositiveStrandAlign convertAlignToPositiveStrand(const Align &align, size_t SAi);
        void createWindows(const std::vector<Align>& alignments);
        void assignAlignmentsToWindows(const std::vector<Align>& alignments);
        void generateTranscripts();

        inline void assignSingleAlignment(Window& win, const PositiveStrandAlign& a) const;


        void identifyAnchors(std::vector<Align>& alignments);


        //DP stitching
        void stitchWindowAligns(Window &window);

        void stitchWindowsAlignNew(Window &window);

        StitchingRecord stitchingBetweenWindowAligns(const WindowAlign& a1,const WindowAlign& a2,int windowDir);

        ExtensionRecord extendWindowAlign(const WindowAlign& a, int windowDir,int extendDir);

        void refreshWinBinMap();

        void finalizeTranscript(Transcript &t, const int chrId);



        int64_t stitchAlignmentToTranscript(
                WindowAlign& alignment,       // alignment to stitch
                Transcript& transcript           // current transcript
        );


        static inline std::pair<int, int64_t> checkJunctionMotif(const std::string &genomeSeq, int64_t leftPos, int64_t rightPos);




        stitchingConfig config_;
        std::vector<Window> windows_;
        std::vector<Transcript> transcripts_;
        std::vector<SJDBOutput> sjdb;
        int64_t outFilterScoreMin_;
        int outFilterMatchMin_;

        WindowAlign* allWindowAligns_; // window alignments stored for all windows.
        int32_t* winBinMap_[2];

        StitchingRecord* nowStitchingRecord_; // stitching records for current window
        ExtensionRecord* nowExtensionRecord_[2]; // extension records for current window, forward & backward
        RawTranscript* nowRawTranscript_; // raw transcripts for current window, used to store stitching records

        int64_t maxTranscriptScore_ ;

        Transcript* goodTranscripts_;
        int numGoodTranscripts_ = 0;


        const GenomeIndexPrefix &genomeIndex_;
        ReadPtr read_;
        TranscriptPtr bestTranscript_;


        uint32_t totalAnchors_ = 0;
        // todo move to config
        static constexpr int MATCH_SCORE = 1;
        static constexpr int MISMATCH_PENALTY = -1;
        static constexpr int GAP_OPEN_PENALTY = -4;
        static constexpr int DEL_OPEN_PENALTY = -4;
        static constexpr int DEL_EXTEND_PENALTY = -1;
        static constexpr int INS_OPEN_PENALTY = -4;
        static constexpr int INS_EXTEND_PENALTY = -1;
        static constexpr int SCORE_STITCH_SJ_SHIFT = 1;// maximum score reduction while searching for SJ boundaries
        static constexpr int SCORE_GAP_GCAG = -4; // GC/AG and CT/GC junction penalty
        static constexpr int SCORE_GAP_ATAC = -8; // AT/AC and GT/AT junction penalty
        static constexpr int SCORE_GAP_NON_CANONICAL = -8;
        static constexpr int MAX_SJ_REPEAT_SEARCH = 255;
        static constexpr uint32_t MIN_INTRON_LENGTH = 20;  // 最小内含子长度
        static constexpr uint32_t MAX_INTRON_LENGTH = 100000;  // 最大内含子长度

    };
}
#endif //RNAALIGNREFACTORED_WINDOWMANAGEMENT_H
