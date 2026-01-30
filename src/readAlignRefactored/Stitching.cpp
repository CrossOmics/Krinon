#include <cstring>
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
        windows_.reserve(maxWindows_);
        transcripts_.reserve(transcriptStoredMax_);
    }

    void Stitching::process(std::vector<RefactorProcessing::Align>& aligns, RefactorProcessing::Read *read) {
        if (aligns.size() == 0){
            return;
        }
        read_ = read;
        alignments_ = &aligns;
        identifyAnchors();
        createWindows();
        assignAlignment();

        if (isPaired_) stitchWindow_PairedEnd();
        else stitchWindow_SingleEnd();


    }

    std::pair<WindowAlign,WindowAlign> Stitching::convertAlignToPositiveStrandWindowAlign(const Align& a, int ind) const{
        // todo Add sjdb check
        WindowAlign wa;
        int64_t loc = genomeIndex_.suffixArray_[ind];
        if (loc > genomeIndex_.genome_.genomeLength_){
            //reverse strand
            wa.genomeStart = genomeIndex_.genome_.genomeLength_ *2 - loc - a.length;
            wa.readStart = read_->length - a.readPos - a.length;
            wa.length = a.length;
            wa.direction = 1 - a.direction;
            wa.isAnchor = a.isAnchor;
            wa.iFragment = a.iFragment;
        }else {
            //forward strand
            wa.genomeStart = loc;
            wa.readStart = a.readPos;
            wa.length = a.length;
            wa.direction = a.direction;
            wa.isAnchor = a.isAnchor;
            wa.iFragment = a.iFragment;
        }

        return {wa,WindowAlign()};

    }

    void Stitching::identifyAnchors() {
        // select anchors based on rep
        for (auto &align : *alignments_) {
            if (align.rep <= maxAnchorRep_) {
                align.isAnchor = true;
            } else {
                align.isAnchor = false;
            }
        }
    }

    void Stitching::createWindowFromAnchor(const RefactorProcessing::WindowAlign &anchor) {
        auto location = anchor.genomeStart;
        auto chrId = genomeIndex_.genome_.getPosChrIndex(location);
        Chromosome chr = genomeIndex_.genome_.chromosomes_[chrId];
        int chrStartBin = chr.start >> winBinSizeLog_;
        int chrEndBin = (chr.start + chr.length - 1) >> winBinSizeLog_;

        Window window;
        int dir = anchor.direction;
        window.chrIndex = chrId;
        window.direction = dir;

        int anchorBin = location >> winBinSizeLog_;
        if (winBinMap_[dir][anchorBin] != -1) {
            // window already exists
            return;
        }

        // extend left: try merge with existing windows
        int leftOverlap = -1;
        int nowBin;
        int leftBin = std::max(chrStartBin, anchorBin - winAnchorDistBins_);
        for (nowBin = anchorBin -1; nowBin >= leftBin; --nowBin) {
            int winId = winBinMap_[dir][nowBin];
            if (winId != -1) {
                // found existing window to merge
                leftOverlap = winId;
                break;
            }
        }
        int32_t nowWindowIndex;
        if (leftOverlap != -1) {
            // merge with existing window
            nowWindowIndex = leftOverlap;
            for (nowBin = nowBin + 1; nowBin <= anchorBin; nowBin++) {
                winBinMap_[dir][nowBin] = nowWindowIndex;
            }

        } else {
            //create a new window
            nowWindowIndex = windows_.size();
        }
        winBinMap_[dir][nowBin] = nowWindowIndex;

        // try merging with right existing windows
        int64_t rightBound = std::min(chrEndBin, nowBin + winAnchorDistBins_);
        int32_t rightOverlap = -1;
        for (nowBin = nowBin + 1; nowBin <= rightBound; nowBin++) {
            if (winBinMap_[dir][nowBin] != -1) {
                rightOverlap = winBinMap_[dir][nowBin];
            }
        }
        if (rightOverlap != -1) {
            nowBin = anchorBin + 1;
            //extend to reach the overlapping window
            while (winBinMap_[dir][nowBin] != rightOverlap) {
                winBinMap_[dir][nowBin] = nowWindowIndex;
                nowBin++;
            }
            // merge with right window
            while (winBinMap_[dir][nowBin] == rightOverlap) {
                winBinMap_[dir][nowBin] = nowWindowIndex;
                nowBin++;
            }

            window.endBin = nowBin;
            // kill right window
            windows_[rightOverlap].startBin = 1;
            windows_[rightOverlap].endBin = 0;
        } else {
            window.endBin = anchorBin;
        }
        window.startBin = anchorBin;
        // update window information
        if (leftOverlap == -1) {
            //new window
            windows_.emplace_back(window);
        } else {
            windows_[leftOverlap].endBin = window.endBin;
        }



    }

    void Stitching::createWindows() {
        // create and try to extend windows

        for (const auto &align : *alignments_) {
            if (!align.isAnchor) continue;

            for (size_t i = align.leftSAIndex; i <= align.rightSAIndex; ++i) {
                // get positive window align
                // todo handle sjdb
                auto [anchor, dummy] = convertAlignToPositiveStrandWindowAlign(align, i);
                createWindowFromAnchor(anchor);
            }
        }

        // flank existing windows
        for (int i = 0; i < windows_.size();++i) {
            auto &win = windows_[i];
            if (win.startBin > win.endBin) continue;
            // flank the window
            auto &chr = genomeIndex_.genome_.chromosomes_[win.chrIndex];
            int chrStartBin = chr.start >> winBinSizeLog_;
            int chrEndBin = (chr.start + chr.length - 1) >> winBinSizeLog_;
            int64_t leftBin = std::max(chrStartBin, win.startBin - flankSize_);
            int64_t rightBin = std::min(chrEndBin, win.endBin + flankSize_);
            int64_t leftKeyBound = win.startBin;
            int dir = win.direction;
            for (int64_t binKey = leftBin; binKey < leftKeyBound; binKey++) {
                winBinMap_[dir][binKey] = i;
            }
            int64_t rightKeyBound = rightBin;
            for (int64_t binKey = win.endBin; binKey <= rightKeyBound; binKey++) {
                winBinMap_[dir][binKey] = i;
            }
            win.startBin = leftBin;
            win.endBin = rightBin;
            // reserve space for alignments
            win.aligns = windowAlignments_.data() + i * maxSeedPerWindows_;
            memset(win.aligns, 0, sizeof(WindowAlign) * maxSeedPerWindows_);
            win.numAligns = 0;
        }
    }

    bool Window::assignAlignment(const WindowAlign &a,int maxSeedPerWindows) {
        // when window is full, check if it needs to be replaced
        if (a.length <= minLengthWhenFull && !a.isAnchor) return false; // ignore too short no anchor alignment

        //detect overlap
        for (int i = 0; i < numAligns; ++i) {
            if (aligns[i].isj == a.isj && a.genomeStart + aligns[i].readStart == aligns[i].genomeStart + a.readStart \
 && (((a.readStart >= aligns[i].readStart) && a.readStart < aligns[i].readStart + aligns[i].length)\
 || (a.readStart + a.length > aligns[i].readStart &&
     a.readStart + a.length <= aligns[i].readStart + aligns[i].length))) {
                // overlap

                // same
                if (a.genomeStart == aligns[i].genomeStart && a.length == aligns[i].length) return false;


                if (a.length > aligns[i].length) {
                    //delete the old alignment and insert the new one
                    //keep the array sorted by readStart

                    //find the position to insert
                    //only need to search one side

                    if (a.readStart > aligns[i].readStart) {
                        for (int j = i + 1; j < numAligns; ++j) {
                            if (aligns[j].readStart > a.readStart) {
                                aligns[j - 1] = a;
                                return true;
                            } else aligns[j - 1] = aligns[j];
                        }
                    } else {
                        for (int j = i - 1; j >= 0; --j) {
                            if (aligns[j].readStart <= a.readStart) {
                                aligns[j + 1] = a;
                                return true;
                            } else aligns[j + 1] = aligns[j];
                        }

                        // reaching here means that a.readStart is the smallest, insert at the beginning
                        aligns[0] = a;
                        return true;


                    }
                }
                // if not, do nothing
                return false;
            }
        }

        if (a.isAnchor) ++numAnchors; // anchor must be added to the window

        // handle the case that there are too many seeds in the window
        if (numAligns == maxSeedPerWindows) {
            // calculate minLengthWhenFull and the alignment to remove
            minLengthWhenFull = std::numeric_limits<int>::max();
            int removePos = -1;
            for (int i = 0; i < numAligns; ++i) {
                if (!aligns[i].isAnchor && aligns[i].length < minLengthWhenFull) {
                    minLengthWhenFull = aligns[i].length;
                    removePos = i;
                }
            }
            if (removePos == -1) return false; // all are anchors, cannot add new alignment


            //remove all shortest non-anchor alignments and insert the new one

            int newAlignNum = 0;
            int posToInsert = -1;
            for (int i = 0; i < numAligns; ++i) {
                if (aligns[i].length > minLengthWhenFull || aligns[i].isAnchor) {
                    if (aligns[i].readStart >= a.readStart && posToInsert == -1) {
                        posToInsert = newAlignNum++;
                        aligns[newAlignNum] = aligns[i];
                    } else {
                        aligns[newAlignNum++] = aligns[i];
                    }
                }
            }
            numAligns = newAlignNum;
            if (minLengthWhenFull >= a.length && !a.isAnchor) return false; // new alignment is too short, ignore
            if (posToInsert == -1) posToInsert = newAlignNum++;
            // reaching here means that a.readStart is the largest



            for (int i = newAlignNum - 1; i > posToInsert; --i) {
                aligns[i] = aligns[i - 1];
            }
            aligns[posToInsert] = a;
            numAligns = newAlignNum;

            return true;
        }


        // simple case
        // find the position to insert
        ++numAligns;
        for (int i = numAligns; i > 0; --i) {
            if (aligns[i - 1].readStart <= a.readStart) {
                aligns[i] = a;
                return true;
            } else aligns[i] = aligns[i - 1];
        }
        // reaching here means that a.readStart is the smallest, insert at the beginning
        aligns[0] = a;
        return true;
    }

    void Stitching::assignAlignment() {
        for (const auto &align : *alignments_) {
            for (size_t i = align.leftSAIndex; i <= align.rightSAIndex; ++i) {
                // get positive window align
                // todo handle sjdb
                auto [wa, dummy] = convertAlignToPositiveStrandWindowAlign(align, i);
                size_t bin = (wa.genomeStart >> winBinSizeLog_);
                auto winId = winBinMap_[wa.direction][bin];
                if (winId == -1) continue;
                Window &win = windows_[winId];
                win.assignAlignment(wa, maxSeedPerWindows_);
            }
        }
    }

    void Stitching::refreshWinBinMap() {
        for (const auto &win: windows_) {
            if (win.startBin > win.endBin) continue;
             auto& winBinArray = winBinMap_[win.direction];
            for (int64_t winBin = win.startBin; winBin <= win.endBin; ++winBin) {
                winBinArray[winBin] = -1; // reset bin map
            }

        }
    }


}