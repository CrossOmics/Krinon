#include "Stitching.h"
#include "Read.h"
#include "Transcript.h"
namespace RefactorProcessing{
    void Stitching::setParam(const Parameters &P) {

        isPaired_ = P.isPaired;
        maxAnchorRep_ = P.maxAnchorRep;
        winBinSizeLog_ = P.winBinSizeLog;
        winAnchorDistBins_ = P.winAnchorDistBins;
        flankSize_ = P.flankSize;
        maxWindows_ = P.maxWindows;
        maxSeedPerWindows_ = P.maxSeedPerWindows;
        maxRep_ = P.maxRep;
        maxExons_ = P.maxExons;
        transcriptStoredMax_ = P.transcriptStoredMax;

        outFilterMultimapMax_ = P.outFilterMultimapMax;;
        maxMismatch_ = P.maxMismatch;
        multimapScoreRange_ = P.multimapScoreRange;
        outFilterScoreMinOverLRead_ = P.outFilterScoreMinOverLRead;
        outFilterMatchMinOverLRead_ = P.outFilterMatchMinOverLRead;

        MATCH_SCORE_ = P.MATCH_SCORE;
        MISMATCH_PENALTY_ = P.MISMATCH_PENALTY;
        GAP_OPEN_PENALTY_ = P.GAP_OPEN_PENALTY;
        DEL_OPEN_PENALTY_ = P.DEL_OPEN_PENALTY;
        DEL_EXTEND_PENALTY_ = P.DEL_EXTEND_PENALTY;
        INS_OPEN_PENALTY_ = P.INS_OPEN_PENALTY;
        INS_EXTEND_PENALTY_ = P.INS_EXTEND_PENALTY;
        SCORE_STITCH_SJ_SHIFT_ = 1;
        SCORE_GAP_GCAG_ = P.SCORE_GAP_GCAG;
        SCORE_GAP_ATAC_ = P.SCORE_GAP_ATAC;
        SCORE_GAP_NON_CANONICAL_ = P.SCORE_GAP_NON_CANONICAL;
        SCORE_ANNOTATED_SJ_ = P.SCORE_ANNOTATED_SJ;
        MAX_SJ_REPEAT_SEARCH_ = P.MAX_SJ_REPEAT_SEARCH;
        MIN_INTRON_LENGTH_ = P.MIN_INTRON_LENGTH;
        MAX_INTRON_LENGTH_ = P.MAX_INTRON_LENGTH;
        for (int i = 0; i < 4; ++i) {
            MAX_MISMATCH_FOR_SJ_[i] = P.MAX_MISMATCH_FOR_SJ[i];
        }

    }

    void Stitching::init() {
        // allocate memory for data structures based on parameters
    }


}