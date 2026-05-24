#ifndef RNAALIGNREFACTORED_STITCHINGR_H
#define RNAALIGNREFACTORED_STITCHINGR_H
//todo
#include <vector>
#include "../genomeRefactored/GenomeIndex.h"
#include "../io/Parameters.h"
#include "Transcript.h"
#include "../utilsRefactored/seqFunctions.h"
#include "Read.h"

#define transcriptMinScore -100000

namespace RefactorProcessing {
    struct WindowAlign{
        int64_t readStart{0};
        int64_t genomeStart{0};
        int64_t length{0};
        bool isAnchor{false};
        bool isAcceptor{false};
        int direction{0};
        int isj{-1};// annotation index
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
        int matches{0}; // ADDITIONAL matches, including the latter alignment's matches
        int mismatches{0}; //ADDITIONAL mismatches. same as above
        int64_t formerExonLengthShift{0}; // length extension of the former exon
        int64_t latterExonLengthShift{0}; // length extension of the latter exon
        int shiftLeft{0};
        int shiftRight{0};
        int64_t isj{-1}; 

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
    private:

    public:
        WindowAlign* aligns{nullptr}; // todo try vector
        int chrIndex{0};
        int direction{0};
        int numAnchors{0};
        int numAligns{0};
        int startBin{0};
        int endBin{0};
        int minLengthWhenFull{0};
        int numSjdbIds{0};
        int *sjdbIds;
        int numFirstFragAligns{0};
        int numSecondFragAligns{0};
        bool assignAlignment(const WindowAlign& a,int maxSeedPerWindows);


    };

    class RawTranscript {
    public:
        int previousTranscriptId{-1};
        int newAlignId{0};
        int score{0};
        int mismatches{0};
        int matches{0};
        int numExon{0};
        int extendedLengthForward{0}; //extension 5' end
        int extendedLengthBackward{0}; // extension 3' end

        int lastExonLength{0}; // to handle short exon
        // todo maybe no need
        int StartAlignId{0}; // to calc max extension

        int sjStrand{-1}; // handle sj strand consistency, -1 for unknown, 0 for +, 1 for -

        void init(const Window& win,int i,int matchScore);

    };

    class RawTranscriptPaired{
    public:
        int previousTranscriptId{-1};
        int newAlignId{0};
        int score{0};
        int mismatches{0};
        int matches{0};
        int numExon{0};
        int StartAlignId{0}; // to calc max extension
        int iFragment{0};
        void init(const Window& win,int i,int matchScore);


    };



    struct fragmentMatchRecord{
        int fragId0{0};
        int fragId1{0};
        // L0 -> F1 -> F0 -> L1
        int extendLengthFormer0{0};
        int extendLengthLatter0{0}; // extend inside
        int extendLengthFormer1{0}; // extend inside
        int extendLengthLatter1{0};
        int score{0};
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

        int outFilterScoreMin_;// calculated based on read length
        int outFilterMatchMin_;

        //sjdb
        int sjdbOverhang_;
        int sjdbLength_;


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
        const GenomeIndex& genomeIndex_;

        std::vector<Align>& alignments_;//  alignments to be stitched

        std::vector<WindowAlign> windowAlignments_;

        std::vector<Window> windows_;

        std::vector<int> windowSjdbIds_;

        std::vector<int> winBinMap_[2];

        std::vector<StitchRecord> stitchRecords_;// stitching records between window-alignments

        std::vector<ExtendRecord> extendRecords_[2];// extension records, forward & backward

        std::vector<ExtendRecord::singleExtendRecord> allSingleExtensionRecord_;// memory pool for single extension record

        Read* read_;

        std::vector<RawTranscript> rawTranscripts_;

        std::vector<RawTranscriptPaired> rawTranscriptsPaired_;



        std::vector<Transcript> transcripts_;

        std::vector<fragmentMatchRecord> fragmentMatchRecords_;

        int maxTranscriptScore_{0};

        int numGoodTranscripts_{0};

        Stitching(const GenomeIndex& g, std::vector<Align>& a):genomeIndex_(g),alignments_(a){};

        // parameters
        void setParam(const Parameters &P);

        // stitching process
        void init();

        inline bool convertAlignToPositiveStrandWindowAlign(const Align& a, size_t ind,WindowAlign &wa1, WindowAlign &wa2) const;

        void process(Read* read, std::string& outputBuffer);

        void identifyAnchors();

        void createWindowFromAnchor(const WindowAlign& anchor);

        void createWindows();

        Transcript convertRawTranscriptToTranscript(const RawTranscript& rt,const Window& win,int startAlignId);



        inline std::pair<int, int64_t>
        checkJunctionMotif(const std::string &genomeSeq, int64_t leftPos, int64_t rightPos);

        void stitchWindowAlign(const Window& window, const WindowAlign& a1, const WindowAlign& a2, StitchRecord& record);

        void extendWindowAlign(const Window& window,const WindowAlign& a,ExtendRecord& res,int extendDir);

        void stitchPaired(const Window& window);

        void stitchSingle(const Window& window);

        void assignAlignment();

        void convertToResult(std::string& outputBuffer);

        // clear data
        void clear();

        void refreshWinBinMap();

        inline int sjdbCheck(const Window& win,const int64_t start,const int64_t end) const;



    };
}

#endif //RNAALIGNREFACTORED_STITCHINGR_H
