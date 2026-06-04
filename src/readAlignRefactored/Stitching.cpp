#include <cstring>
#include <cmath>
#include "Stitching.h"
#include "Read.h"
#include "Transcript.h"
#include "../log/ErrorRecord.h"

#define SIMPLE_DELETION -1;
#define SIMPLE_INSERTION -2;
namespace RefactorProcessing {


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

        sjdbOverhang_ = P.sjdbOverhang;
        sjdbLength_ = P.sjdbLength;

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
        windowAlignments_.resize(maxWindows_ * maxSeedPerWindows_);
        windowSjdbIds_.resize(maxWindows_ * maxSeedPerWindows_/2);
        transcripts_.resize(transcriptStoredMax_);
        size_t winBinNum = (genomeIndex_.genome_.genomeLength_ >> (winBinSizeLog_ - 1)) + 2;
        winBinMap_[0].resize(winBinNum, -1);
        winBinMap_[1].resize(winBinNum, -1);
        size_t maxStitchRecordNum = maxSeedPerWindows_ * maxSeedPerWindows_;
        stitchRecords_.resize(maxStitchRecordNum);
        extendRecords_[0].resize(maxSeedPerWindows_);
        extendRecords_[1].resize(maxSeedPerWindows_);
        allSingleExtensionRecord_.resize(maxSeedPerWindows_ * 2 * (1 + maxMismatch_));
        for (int i = 0; i < maxSeedPerWindows_; ++i) {
            extendRecords_[0][i].maxExtensionLengthWithMismatch =
                    allSingleExtensionRecord_.data() + (i * (1 + maxMismatch_));
            extendRecords_[1][i].maxExtensionLengthWithMismatch =
                    allSingleExtensionRecord_.data() + ((i + maxSeedPerWindows_) * (1 + maxMismatch_));
        }
        size_t maxRawTranscriptNum = maxSeedPerWindows_ * maxSeedPerWindows_;
        if (!isPaired_) rawTranscripts_.resize(maxRawTranscriptNum);
        else rawTranscriptsPaired_.resize(maxRawTranscriptNum);
        transcripts_.resize(100);
        if (isPaired_) size_t maxFragmentMatchRecordNum = maxSeedPerWindows_ * maxSeedPerWindows_;
    }

    void Stitching::clear() {
        refreshWinBinMap();
        windows_.clear();
        numGoodTranscripts_ = 0;
    }

    void Stitching::process(RefactorProcessing::Read *read, std::string &outputBuffer) {
        if (alignments_.size() == 0) {
            return;
        }
        maxTranscriptScore_ = 0;
        numGoodTranscripts_ = 0;
        read_ = read;
        outFilterScoreMin_ = int(double(read_->length - 1) * outFilterScoreMinOverLRead_);
        outFilterMatchMin_ = int(double(read_->length - 1) * outFilterMatchMinOverLRead_);

        identifyAnchors();
        createWindows();
        assignAlignment();

        if (isPaired_) {
            for (auto &window: windows_) {
                if (window.startBin > window.endBin) continue;
                stitchPaired(window);
            }
        } else {
            for (auto &window: windows_) {
                if (window.startBin > window.endBin) continue;
                stitchSingle(window);
            }
        }


        int trueNumGoodTranscripts = numGoodTranscripts_;

        for (int i = 0; i < numGoodTranscripts_; ++i) {
            if (transcripts_[i].score < maxTranscriptScore_ - multimapScoreRange_) {
                trueNumGoodTranscripts = i;
                break;
            }
        }

        if (trueNumGoodTranscripts > outFilterMultimapMax_) {
            trueNumGoodTranscripts = 0;
            clear();
            return;
        }
        numGoodTranscripts_ = trueNumGoodTranscripts;

        convertToResult(outputBuffer);

        clear();
    }

    void Stitching::convertToResult(std::string &outputBuffer) {
        for (int i = 0; i < numGoodTranscripts_; ++i) {
            transcripts_[i].convertToSAM(*read_, isPaired_, numGoodTranscripts_, outputBuffer);
        }
    }

    inline bool
    Stitching::convertAlignToPositiveStrandWindowAlign(const Align &a, size_t ind, WindowAlign &wa,
                                                       WindowAlign &wa2) const {
        /**
         * TODO: If `loc` can take negative values, comparing it with `sjdbSeqLength_`
         *       is only valid of casting `sjdbSeqLength_` to signed value does not
         *       overflow.
         *       I will put a runtime check for that and abort if it fails.
         */
        /**
         *  loc is non-negative
         *  Also, sjdbSeqLength_ is at most 2 * sjdbNum_ * sjdbLength_, which is very unlikely to overflow 64 bits
         */
        /*uint64_t max_int64 = static_cast<uint64_t>(std::numeric_limits<int64_t>::max());
        if (genomeIndex_.genome_.sjdbSeqLength_ > max_int64) {
            std::cerr << "FATAL: `sjdbSeqLength_` overflow!";
            std::abort();
        }*/
        //sjdb check
        wa2.direction = 2;// 2 means dummy
        size_t loc = genomeIndex_.suffixArray_[ind];
        if (loc >= genomeIndex_.genome_.sjdbStart_ && genomeIndex_.genome_.sjdbNum_ > 0) {
            // maybe a cross-sjdb alignment
            loc -= genomeIndex_.genome_.sjdbStart_;
            int dir = a.direction;
            int64_t readPos = a.readPos;
            if (loc > genomeIndex_.genome_.sjdbSeqLength_) {
                dir = 1 - a.direction;
                loc = 2 * genomeIndex_.genome_.sjdbSeqLength_ - loc - a.length;
                readPos = read_->length - a.readPos - a.length;
            }

            size_t startInSj = loc % sjdbLength_;
            if (startInSj < (size_t) sjdbOverhang_ && startInSj + a.length > (size_t) sjdbOverhang_) {
                // crossing sjdb
                int sjIndex = loc / sjdbLength_;
                int64_t donorStart = genomeIndex_.genome_.sjDonorStart_[sjIndex] + startInSj;
                int64_t donorLength = sjdbOverhang_ - startInSj;
                int64_t acceptorStart = genomeIndex_.genome_.sjAcceptorStart_[sjIndex];
                int64_t acceptorLength = a.length - donorLength;
                wa.genomeStart = donorStart;
                wa.readStart = readPos;
                wa.length = donorLength;
                wa.direction = dir;
                wa.isAnchor = a.isAnchor;
                wa.iFragment = a.iFragment;
                wa.isj = sjIndex;

                wa2.genomeStart = acceptorStart;
                wa2.readStart = readPos + donorLength;
                wa2.length = acceptorLength;
                wa2.direction = dir;
                wa2.isAnchor = a.isAnchor;
                wa2.iFragment = a.iFragment;
                wa2.isj = sjIndex;
                wa2.isAcceptor = true;
                return true;
            } else {
                // not crossing sjdb
                wa.direction = 2;
                return false;
            }
        }

        if (loc > genomeIndex_.genome_.genomeLength_) {
            //reverse strand
            wa.genomeStart = genomeIndex_.genome_.genomeLength_ * 2 - loc - a.length;
            wa.readStart = read_->length - a.readPos - a.length;
            wa.length = a.length;
            wa.direction = 1 - a.direction;
            wa.isAnchor = a.isAnchor;
            wa.iFragment = a.iFragment;
            wa.isj = -1;
        } else {
            //forward strand
            wa.genomeStart = loc;
            wa.readStart = a.readPos;
            wa.length = a.length;
            wa.direction = a.direction;
            wa.isAnchor = a.isAnchor;
            wa.iFragment = a.iFragment;
            wa.isj = -1;
        }

        return true;

    }

    void Stitching::identifyAnchors() {
        // select anchors based on rep
        for (auto &align: alignments_) {
            if (align.rep <= maxAnchorRep_) {
                align.isAnchor = true;
            } else {
                align.isAnchor = false;
            }
        }
    }

    void Stitching::createWindowFromAnchor(const RefactorProcessing::WindowAlign &anchor) {
        if (windows_.size() >= (size_t) maxWindows_) return;
        auto location = anchor.genomeStart;
        auto chrId = genomeIndex_.genome_.getPosChrIndex(location);
        const Chromosome &chr = genomeIndex_.genome_.chromosomes_[chrId];
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
        for (nowBin = anchorBin - 1; nowBin >= leftBin; --nowBin) {
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
        winBinMap_[dir][anchorBin] = nowWindowIndex;

        // try merging with right existing windows
        int64_t rightBound = std::min(chrEndBin, anchorBin + winAnchorDistBins_);
        int32_t rightOverlap = -1;
        for (nowBin = anchorBin + 1; nowBin <= rightBound; nowBin++) {
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
        size_t numAlignments = alignments_.size();
        WindowAlign anchor, dummy;

        for (size_t iA = 0; iA < numAlignments; ++iA) {
            auto &align = alignments_[iA];
            if (!align.isAnchor) continue;

            for (size_t i = align.leftSAIndex; i <= align.rightSAIndex; ++i) {
                // get positive window align
                // handle sjdb
                bool valid = convertAlignToPositiveStrandWindowAlign(align, i, anchor, dummy);

                if (!valid) continue; // not a valid alignment
                createWindowFromAnchor(anchor);
                // todo ? no need to create anchor2 for a pair of cross-sjdb aligns
            }
        }

        // flank existing windows
        for (size_t i = 0; i < windows_.size(); ++i) {
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
            // TODO: Are you sure this is correct? these pointer arithmetics look risky ...
            /**
             * fixed, I forgot to limit the number of windows created in createWindowFromAnchor
             */
            win.aligns = windowAlignments_.data() + i * maxSeedPerWindows_;
            win.sjdbIds = windowSjdbIds_.data() + i * maxSeedPerWindows_/2;
            memset((void *) win.aligns, 0, sizeof(WindowAlign) * maxSeedPerWindows_);
            win.numAligns = 0;
            win.numAnchors = 0;
            win.numSjdbIds = 0;
            win.numFirstFragAligns = 0;
            win.numSecondFragAligns = 0;
            
        }
    }

    inline int Stitching::sjdbCheck(const Window &win,const int64_t start, const int64_t end) const{
        // check if the del is in the GTF file, and return the sjdb index, otherwise return -1
        for (int i = 0; i < win.numSjdbIds; ++i) {
            int sjdbId = win.sjdbIds[i];
            if ((size_t) start == genomeIndex_.genome_.sjDataBase_[sjdbId].start && (size_t) end == genomeIndex_.genome_.sjDataBase_[sjdbId].end) {
                return sjdbId;
            }
        }
        return -1;
    }

    bool Window::assignAlignment(const WindowAlign &a, int maxSeedPerWindows) {
        // when window is full, check if it needs to be replaced
        if (a.length <= minLengthWhenFull && !a.isAnchor) return false; // ignore too short no anchor alignment

        //bool overlapInserted = false;
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

        if (a.isj != -1) {
            bool sjdbExists = false;
            for (int i  = 0; i < numSjdbIds; ++i) {
                if (sjdbIds[i] == a.isj){
                    sjdbExists = true;
                    break;
                }
            }
            if (!sjdbExists) {
                ++numSjdbIds;
                sjdbIds[numSjdbIds - 1] = a.isj;
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
        for (int i = numAligns - 1; i > 0; --i) {
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
        size_t numAlignments = alignments_.size();
        WindowAlign wa, wa2;
        for (size_t iA = 0; iA < numAlignments; ++iA) {
            auto &align = alignments_[iA];
            // ignore too repetitive alignments
            if (align.rep > maxRep_) continue;
            for (size_t i = align.leftSAIndex; i <= align.rightSAIndex; ++i) {
                // get positive window align
                // handle sjdb
                bool valid = convertAlignToPositiveStrandWindowAlign(align, i, wa, wa2);
                if (!valid) continue;
                size_t bin = (wa.genomeStart >> winBinSizeLog_);
                auto winId = winBinMap_[wa.direction][bin];
                if (winId == -1) continue;
                Window &win = windows_[winId];
                win.assignAlignment(wa, maxSeedPerWindows_);
                if (wa2.direction != 2) {
                    size_t bin2 = (wa2.genomeStart >> winBinSizeLog_);
                    auto winId2 = winBinMap_[wa2.direction][bin2];
                    if (winId2 == -1) continue;
                    Window &win2 = windows_[winId2];
                    win2.assignAlignment(wa2, maxSeedPerWindows_);
                }
            }
        }
    }

    void Stitching::refreshWinBinMap() {
        for (const auto &win: windows_) {
            if (win.startBin > win.endBin) continue;
            auto &winBinArray = winBinMap_[win.direction];
            for (int64_t winBin = win.startBin; winBin <= win.endBin; ++winBin) {
                winBinArray[winBin] = -1; // reset bin map
            }
        }
    }

    inline std::pair<int, int64_t> Stitching::checkJunctionMotif(const std::string &genomeSeq,
                                                                 int64_t leftPos,
                                                                 int64_t rightPos) {
        // may optimize this function later, use hash and array
        int left = charToIndex(genomeSeq[leftPos + 1]) * 5 + charToIndex(genomeSeq[leftPos + 2]);
        int right = charToIndex(genomeSeq[rightPos - 1]) * 5 + charToIndex(genomeSeq[rightPos]);
        // check the left motif
        if (left == 13) { // GT
            if (right == 2) // AG
                return {1, 0};
            else if (right == 3) // AT
                return {6, SCORE_GAP_ATAC_};
        } else if (left == 8) { //CT
            if (right == 1) //AC
                return {2, 0};
            else if (right == 11) //GC
                return {4, SCORE_GAP_GCAG_};
        } else if (left == 3) { //AT
            if (right == 1) // AC
                return {5, SCORE_GAP_ATAC_};
        } else if (left == 11) { //GC
            if (right == 2)
                return {3, SCORE_GAP_GCAG_};
        }
        return {0, SCORE_GAP_NON_CANONICAL_};
    }

    void Stitching::stitchWindowAlign(const Window &window, const WindowAlign &a1, const WindowAlign &a2,
                                      StitchRecord &record) {
        //calculate the stitching score between two alignments
        //stitch a2 to a1

        const int64_t a2ReadEnd = a2.readStart + a2.length - 1;
        const int64_t a2GenomeEnd = a2.genomeStart + a2.length - 1;
        const int64_t a1ReadEnd = a1.readStart + a1.length - 1;
        const int64_t a1GenomeEnd = a1.genomeStart + a1.length - 1;

        record.isj = -1;

        if (a1.iFragment != a2.iFragment) {
            // stitching between paired reads
            // handle mates extension while DP
            record.type = StitchRecord::CROSS_FRAGMENTS;
            record.score = a2.length * MATCH_SCORE_;
            record.matches = a2.length;
            record.mismatches = 0;
            record.spliceJunctionType = -3;
            record.formerExonLengthShift = 0;
            record.latterExonLengthShift = 0;
            return;
        }


        //sjdb check
        if (a1.isj != -1 && a1.isj == a2.isj && a2.readStart == a1ReadEnd + 1 && a1GenomeEnd + 1 < a2.genomeStart) {
            // annotated splice junction
            if (genomeIndex_.genome_.sjDataBase_[a2.isj].motif == 0
                && (a2.length <= genomeIndex_.genome_.sjDataBase_[a2.isj].shiftRight ||
                    a1.length <= genomeIndex_.genome_.sjDataBase_[a2.isj].shiftLeft)) {
                // too large repeats around non-canonical
                record.type = StitchRecord::CANNOT_STITCH;
                return;
            }
            record.type = StitchRecord::SPLICE_JUNCTION;
            record.score = SCORE_ANNOTATED_SJ_ + a2.length * MATCH_SCORE_;
            record.formerExonLengthShift = 0;
            record.latterExonLengthShift = 0;
            record.matches = a2.length;
            record.mismatches = 0;
            record.spliceJunctionType = genomeIndex_.genome_.sjDataBase_[a1.isj].motif;
            record.shiftLeft = genomeIndex_.genome_.sjDataBase_[a1.isj].shiftLeft;
            record.shiftRight = genomeIndex_.genome_.sjDataBase_[a1.isj].shiftRight;
            record.isj = a1.isj;
            //record.stitchingType = {genomeIndex_.genome->sjdb[a2.isj].motif, a2.genomeStart - (a1GenomeEnd + 1), true};
            return;
        }


        if (a2ReadEnd <= a1ReadEnd || a2GenomeEnd <= a1GenomeEnd) {
            // r or g fully overlap, cannot stitch
            record.type = StitchRecord::CANNOT_STITCH;
            return;
        }

        int64_t a2ReadStart = a2.readStart;
        int64_t a2GenomeStart = a2.genomeStart;
        int64_t a2Length = a2.length;
        const std::string &genomeSeq = genomeIndex_.genome_.sequence_;
        const std::string &readSeq = read_->sequence[window.direction];

        if (a2.readStart <= a1ReadEnd) {
            // a2 5' overlaps a1 3'
            a2GenomeStart += (a1ReadEnd - a2.readStart + 1);
            a2ReadStart = a1ReadEnd + 1;
            a2Length = a2ReadEnd - a2ReadStart + 1;
        }

        int64_t readGap = a2ReadStart - (a1ReadEnd + 1);
        int64_t genomeGap = a2GenomeStart - (a1GenomeEnd + 1);

        int nMatch = 0, nMismatch = 0;
        int64_t Del = 0, Ins = 0; // deletion and insertion length
        int64_t junctionReadPos = 0;
        int64_t junctionType = 0;
        int64_t lastIntronBase = a2GenomeStart - readGap - 1; // if all gap belongs to acceptor

        if (genomeGap == 0 && readGap == 0) {
            // perfect match, no gap
            record.type = StitchRecord::PERFECT_MATCH;
            record.score = (a2Length - a2.length) * MATCH_SCORE_;
            record.formerExonLengthShift = -a1.length;
            record.latterExonLengthShift = a1.length + a2Length - a2.length;
            record.matches = a2Length;
            record.mismatches = 0;
        } else if (genomeGap > 0 && readGap > 0 && genomeGap == readGap) {
            // calc the matches & mismatches in the gap
            for (int64_t i = 1; i <= readGap; ++i) {
                if (genomeSeq[a1GenomeEnd + i] != 'N' && readSeq[a1ReadEnd + i] != 'N') {
                    if (genomeSeq[a1GenomeEnd + i] == readSeq[a1ReadEnd + i]) {
                        nMatch++;
                    } else {
                        nMismatch++;
                    }
                }
            }

            record.type = StitchRecord::SAME_GAP;
            record.score = (nMatch + a2Length) * MATCH_SCORE_ + nMismatch * MISMATCH_PENALTY_;
            record.formerExonLengthShift = -a1.length;
            record.latterExonLengthShift = a1.length + readGap;
            record.matches = nMatch + a2Length;
            record.mismatches = nMismatch;
        } else if (genomeGap > readGap) {
            // means there is a deletion or an intron (long deletion)
            Del = genomeGap - readGap;
            record.score = MATCH_SCORE_ * a2Length;
            if (Del > MAX_INTRON_LENGTH_) {
                record.type = StitchRecord::CANNOT_STITCH;// too large deletion
                return;
            }

            // Now, find the best junction pos
            int64_t junctionReadPosTemp = 1;
            int64_t scoreTemp = 0;
            do {
                junctionReadPosTemp--;
                if (readSeq[a1ReadEnd + junctionReadPosTemp] !=
                    genomeSeq[lastIntronBase + junctionReadPosTemp] &&
                    readSeq[a1ReadEnd + junctionReadPosTemp] ==
                    genomeSeq[a1GenomeEnd + junctionReadPosTemp]) {
                    scoreTemp -= MATCH_SCORE_;
                }
            } while (scoreTemp + SCORE_STITCH_SJ_SHIFT_ >= 0 && a1.length + junctionReadPosTemp > 1);

            int64_t maxDelScore = -1000000;
            int64_t junctionPenalty = 0;
            scoreTemp = 0;

            do {
                if (readSeq[a1ReadEnd + junctionReadPosTemp] !=
                    genomeSeq[lastIntronBase + junctionReadPosTemp] &&
                    readSeq[a1ReadEnd + junctionReadPosTemp] ==
                    genomeSeq[a1GenomeEnd + junctionReadPosTemp]) {
                    scoreTemp += MATCH_SCORE_;
                }
                if (readSeq[a1ReadEnd + junctionReadPosTemp] ==
                    genomeSeq[lastIntronBase + junctionReadPosTemp] &&
                    readSeq[a1ReadEnd + junctionReadPosTemp] !=
                    genomeSeq[a1GenomeEnd + junctionReadPosTemp]) {
                    scoreTemp -= MATCH_SCORE_;
                }

                int junctionTypeTemp = SIMPLE_DELETION;
                int64_t junctionPenaltyTemp = 0;
                int64_t maxDelScoreTemp = scoreTemp;

                if (Del >= MIN_INTRON_LENGTH_) {
                    auto res = checkJunctionMotif(
                            genomeSeq,
                            a1GenomeEnd + junctionReadPosTemp,
                            lastIntronBase + junctionReadPosTemp
                    );
                    junctionTypeTemp = res.first;
                    junctionPenaltyTemp = res.second;
                    maxDelScoreTemp += junctionPenaltyTemp;
                }

                if (maxDelScoreTemp > maxDelScore) {
                    maxDelScore = maxDelScoreTemp;
                    junctionReadPos = junctionReadPosTemp;
                    junctionType = junctionTypeTemp;
                    junctionPenalty = junctionPenaltyTemp;
                }

                junctionReadPosTemp++;
            } while (junctionReadPosTemp < a2ReadEnd - a1ReadEnd);

            //repeat length for junctions is needed to be considered for short exon filtering
            int64_t jjL = 0, jjR = 0;
            while (a1GenomeEnd + junctionReadPos >= jjL &&
                   genomeSeq[a1GenomeEnd - jjL + junctionReadPos] == genomeSeq[lastIntronBase - jjL + junctionReadPos] &&
                   genomeSeq[a1GenomeEnd - jjL + junctionReadPos] != 'N' &&
                   genomeSeq[a1GenomeEnd - jjL + junctionReadPos] != '#' &&
                   jjL <= 255) {//go back
                jjL++;
            }

            while ((size_t) a1GenomeEnd + jjR + junctionReadPos + 1 < genomeSeq.length() &&
                   genomeSeq[a1GenomeEnd + jjR + junctionReadPos + 1] ==
                   genomeSeq[lastIntronBase + jjR + junctionReadPos + 1] &&
                   genomeSeq[a1GenomeEnd + jjR + junctionReadPos + 1] != 'N' &&
                   genomeSeq[a1GenomeEnd + jjR + junctionReadPos + 1] != '#' &&
                   jjR <= 255) {//go forward
                jjR++;
            }

            if (junctionType <= 0){
                // non canonical or deletion
                junctionReadPos-=jjL;
                // todo notice the filtering
                jjR+=jjL;
                jjL=0;
            }



            // score donor and acceptor
            for (int64_t i = std::min((int64_t) 1, junctionReadPos + 1); i <= std::max(readGap, junctionReadPos); ++i) {
                size_t pos = (i <= junctionReadPos) ? (a1GenomeEnd + i) : (lastIntronBase + i);
                if (genomeSeq[pos] != 'N' && readSeq[a1ReadEnd + i] != 'N') {
                    if (genomeSeq[pos] != readSeq[a1ReadEnd + i]) {
                        nMismatch++;
                        record.score -= MATCH_SCORE_;
                        if (i < 1 || i > readGap) {
                            record.score -= MATCH_SCORE_;
                            nMatch--;
                        }
                    } else {
                        if (i >= 1 && i <= readGap) {
                            nMatch++;
                            record.score += MATCH_SCORE_;
                        }
                    }
                }
            }


            if (junctionType < 0 || (MAX_MISMATCH_FOR_SJ_[(junctionType + 1) / 2] == -1 ||
                                     nMismatch <= MAX_MISMATCH_FOR_SJ_[(junctionType + 1) / 2])) {
                // valid junction
            } else {
                // not a valid junction
                record.type = StitchRecord::CANNOT_STITCH;
                return;
            }

            // score the gap
            // check if the junction is annotated
            if (window.numSjdbIds > 0) {
                int isj = sjdbCheck(window,a1GenomeEnd+junctionReadPos +1, lastIntronBase + junctionReadPos);
                if (isj != -1) {
                    // annotated junction
                    if (genomeIndex_.genome_.sjDataBase_[isj].motif == 0){
                        junctionReadPos += genomeIndex_.genome_.sjDataBase_[isj].shiftLeft;
                        jjL = genomeIndex_.genome_.sjDataBase_[isj].shiftLeft;
                        jjR = genomeIndex_.genome_.sjDataBase_[isj].shiftRight;
                    }
                    record.isj = isj;
                    record.type = StitchRecord::SPLICE_JUNCTION;
                    record.matches = nMatch + a2Length;
                    record.mismatches = nMismatch;
                    record.score += SCORE_ANNOTATED_SJ_;
                    record.formerExonLengthShift = junctionReadPos;
                    record.latterExonLengthShift = (a2ReadEnd - a1ReadEnd - junctionReadPos) - a2.length;
                    record.spliceJunctionType = genomeIndex_.genome_.sjDataBase_[isj].motif;
                    record.shiftLeft = jjL;
                    record.shiftRight = jjR;
                    return;
                }
            }

            if (Del > MIN_INTRON_LENGTH_) {
                record.score += GAP_OPEN_PENALTY_ + junctionPenalty;
                record.type = StitchRecord::SPLICE_JUNCTION;
                //int sjStrand = 0;
                //if (junctionType > 0) sjStrand = 2 - (junctionType % 2);
                record.spliceJunctionType = junctionType;

            } else {
                record.score += DEL_OPEN_PENALTY_ + Del * DEL_EXTEND_PENALTY_;
                record.type = StitchRecord::DELETION;
                record.spliceJunctionType = SIMPLE_DELETION;
            }
            record.formerExonLengthShift = junctionReadPos;
            record.latterExonLengthShift = (a2ReadEnd - a1ReadEnd - junctionReadPos) - a2.length;
            record.matches = nMatch + a2Length;
            record.mismatches = nMismatch;
            record.shiftLeft = jjL;
            record.shiftRight = jjR;


        } else if (readGap > genomeGap) {
            //insert
            record.score = MATCH_SCORE_ * a2Length;
            Ins = readGap - genomeGap;
            if (genomeGap == 0) {
                junctionReadPos = 0;
            } else if (genomeGap < 0) {
                // overlapping
                junctionReadPos = 0;
                record.score -= (-genomeGap) * MATCH_SCORE_; // reduce score
            } else {

                int64_t maxInsScore = 0;
                int64_t scoreTemp = 0;
                int64_t junctionReadPosTemp = 1;
                // find the best junction pos
                for (junctionReadPosTemp = 1; junctionReadPosTemp <= genomeGap; junctionReadPosTemp++) {
                    if (genomeSeq[a1GenomeEnd + junctionReadPosTemp] != 'N') {
                        scoreTemp += (readSeq[a1ReadEnd + junctionReadPosTemp] ==
                                      genomeSeq[a1GenomeEnd + junctionReadPosTemp]) ? MATCH_SCORE_ : -MATCH_SCORE_;
                        scoreTemp += (readSeq[a1ReadEnd + Ins + junctionReadPosTemp] ==
                                      genomeSeq[a1GenomeEnd + junctionReadPosTemp]) ? -MATCH_SCORE_ : MATCH_SCORE_;
                    }

                    if (scoreTemp > maxInsScore) {
                        maxInsScore = scoreTemp;
                        junctionReadPos = junctionReadPosTemp;
                    }
                }
                // score donor and acceptor
                for (int64_t i = 1; i <= genomeGap; i++) {
                    int64_t rPos = a1ReadEnd + i + (i <= junctionReadPos ? 0 : Ins);
                    if (genomeSeq[a1GenomeEnd + i] != 'N' && readSeq[rPos] != 'N') {
                        if (genomeSeq[a1GenomeEnd + i] != readSeq[rPos]) {
                            nMismatch++;
                            record.score -= MATCH_SCORE_;
                        } else {
                            nMatch++;
                            record.score += MATCH_SCORE_;
                        }
                    }
                }

                // todo implement the alignInsertionFlush parameter

            }

            record.score += Ins * INS_EXTEND_PENALTY_ + INS_OPEN_PENALTY_;
            junctionType = SIMPLE_INSERTION; // mark insertion
            record.type = StitchRecord::INSERTION;
            record.formerExonLengthShift = junctionReadPos;
            record.latterExonLengthShift = (a2ReadEnd - a1ReadEnd - junctionReadPos - Ins) - a2.length;
            record.matches = nMatch + a2Length;
            record.mismatches = nMismatch;
            record.spliceJunctionType = junctionType;
        }
    }

    void Stitching::extendWindowAlign(const Window &window, const WindowAlign &a, ExtendRecord &res, int extendDir) {
        const std::string &genomeSeq = genomeIndex_.genome_.sequence_;
        std::string &readSeq = read_->sequence[window.direction];


        int mismatchCount = 0;
        int matchCount = 0;
        int maxExtendScore = 0;
        int maxScoreMismatch = 0;
        int maxScoreLength = 0;
        int maxScoreMatch = 0;
        int extendScore = 0;
        if (extendDir == 0) {
            // extend forward
            int64_t extendLength = std::min(a.readStart, a.genomeStart);
            int l = 1;
            for (l = 1; l <= extendLength; ++l) {
                size_t genomePos = a.genomeStart - l;
                size_t readPos = a.readStart - l;


                if (genomeSeq[genomePos] == 'N' || readSeq[readPos] == 'N') continue;
                if (genomeSeq[genomePos] == '#') break;
                if (readSeq[readPos] == '#') break; // spacer between paired reads
                if (genomeSeq[genomePos] != readSeq[readPos]) {
                    if (extendScore > maxExtendScore) {
                        maxExtendScore = extendScore;
                        maxScoreLength = l - 1;
                        maxScoreMatch = matchCount;
                        maxScoreMismatch = mismatchCount;
                    }
                    res.maxExtensionLengthWithMismatch[mismatchCount] = {
                            .length = maxScoreLength,
                            .matched = maxScoreMatch,
                            .mismatches = maxScoreMismatch
                    };
                    ++mismatchCount;
                    --extendScore;
                    if (mismatchCount > maxMismatch_) break; // stop extending if too many mismatches
                } else {
                    ++matchCount;
                    ++extendScore;
                }
            }

            if (mismatchCount <= maxMismatch_) {
                if (extendScore > maxExtendScore) {
                    maxExtendScore = extendScore;
                    maxScoreLength = (int) l - 1;
                    maxScoreMatch = matchCount;
                    maxScoreMismatch = mismatchCount;
                }
                res.maxMismatch = mismatchCount;
                res.maxExtensionLengthWithMismatch[mismatchCount] = {
                        .length = maxScoreLength,
                        .matched = maxScoreMatch,
                        .mismatches = maxScoreMismatch
                };
            } else {
                res.maxMismatch = maxMismatch_;
            }
            res.maxExtensionScore = maxExtendScore;

        } else {
            // extend backward
            int64_t extendLength = read_->length - a.readStart - a.length;
            int64_t pieceEndReadPos = a.readStart + a.length - 1;
            int64_t pieceEndGenomePos = a.genomeStart + a.length - 1;
            int l;
            for (l = 1; l <= extendLength; ++l) {
                int64_t readPos = pieceEndReadPos + l;
                int64_t genomePos = pieceEndGenomePos + l;
                if (genomeSeq[genomePos] == 'N' || readSeq[readPos] == 'N') continue;
                if (genomeSeq[genomePos] == '#') break;
                if (readSeq[readPos] == '#') break; // spacer between paired reads
                if (genomeSeq[genomePos] != readSeq[readPos]) {
                    if (extendScore > maxExtendScore) {
                        maxExtendScore = extendScore;
                        maxScoreLength = l - 1;
                        maxScoreMatch = matchCount;
                        maxScoreMismatch = mismatchCount;
                    }
                    res.maxExtensionLengthWithMismatch[mismatchCount] = {
                            .length = maxScoreLength,
                            .matched = maxScoreMatch,
                            .mismatches = maxScoreMismatch
                    };
                    ++mismatchCount;
                    --extendScore;
                    if (mismatchCount > maxMismatch_) break; // stop extending if too many mismatches
                } else {
                    ++matchCount;
                    ++extendScore;
                }
            }

            if (mismatchCount <= maxMismatch_) {
                if (extendScore > maxExtendScore) {
                    maxExtendScore = extendScore;
                    maxScoreLength = (int) l - 1;
                    maxScoreMatch = matchCount;
                    maxScoreMismatch = mismatchCount;
                }
                res.maxMismatch = mismatchCount;
                res.maxExtensionLengthWithMismatch[mismatchCount] = {
                        .length = maxScoreLength,
                        .matched = maxScoreMatch,
                        .mismatches = maxScoreMismatch
                };
            } else {
                res.maxMismatch = maxMismatch_;
            }

            res.maxExtensionScore = maxExtendScore;

        }
    }

    //todo consider moving to a new file
    inline int calcStitchRecPos(int i, int j) {
        //i<j
        return (j * (j - 1)) / 2 + j - i - 1;
    }

    inline int calcRawTranscriptPos(int i, int dist, int num) {

        return i * num + dist;
    }

    void RawTranscript::init(const Window &win, int i, int matchScore) {
        previousTranscriptId = -1;
        score = win.aligns[i].length * matchScore;
        newAlignId = i;
        matches = win.aligns[i].length;
        lastExonLength = win.aligns[i].length;
        mismatches = 0;
        numExon = 1;
        StartAlignId = i;
        sjStrand = -1;
    }

    Transcript
    Stitching::convertRawTranscriptToTranscript(const RawTranscript &rt, const Window &win, int startAlignId) {
        Transcript t;
        t.strand = win.direction;
        t.score = rt.score;
        t.matched = rt.matches;
        t.mismatched = rt.mismatches;


        t.exons.resize(rt.numExon);
        t.sj.resize(rt.numExon - 1);

        int numIns = 0;
        int numDel = 0;
        // collect alignments
        WindowAlign &lastAlign = win.aligns[rt.newAlignId];
        t.exons[rt.numExon - 1] = Exon{
                .genomeStart = lastAlign.genomeStart,
                .length = lastAlign.length + rt.extendedLengthBackward,
                .readStart = lastAlign.readStart
        };

        int nowRawTranscriptId = rt.previousTranscriptId;
        int rightAlignId = rt.newAlignId;
        int alignPos = rt.numExon - 1;

        while (nowRawTranscriptId != -1) {
            const RawTranscript &nowRawT = rawTranscripts_[nowRawTranscriptId];
            const WindowAlign &nowAlign = win.aligns[nowRawT.newAlignId];
            const StitchRecord &nowStitch = stitchRecords_[calcStitchRecPos(nowRawT.newAlignId, rightAlignId)];
            if (nowStitch.type == StitchRecord::PERFECT_MATCH || nowStitch.type == StitchRecord::SAME_GAP) {
                // no new exon
                t.exons[alignPos].length += nowStitch.latterExonLengthShift;
                t.exons[alignPos].readStart -= nowStitch.latterExonLengthShift;
                t.exons[alignPos].genomeStart -= nowStitch.latterExonLengthShift;
            } else {
                --alignPos;
                t.exons[alignPos].genomeStart = nowAlign.genomeStart;
                t.exons[alignPos].length = nowAlign.length + nowStitch.formerExonLengthShift;
                t.exons[alignPos].readStart = nowAlign.readStart;
                t.exons[alignPos + 1].genomeStart -= nowStitch.latterExonLengthShift;
                t.exons[alignPos + 1].length += nowStitch.latterExonLengthShift;
                t.exons[alignPos + 1].readStart -= nowStitch.latterExonLengthShift;

                //calc sj
                int genomeGap = (t.exons[alignPos + 1].genomeStart) -
                                (t.exons[alignPos].genomeStart + t.exons[alignPos].length);
                int readGap =
                        (t.exons[alignPos + 1].readStart) - (t.exons[alignPos].readStart + t.exons[alignPos].length);
                int sjGap = genomeGap - readGap; // -ins, +del or sj

                //todo handle sj collecting later
                t.sj[alignPos] = SpliceJunction{
                        .type = nowStitch.spliceJunctionType,
                        .length = sjGap,
                };
            }

            nowRawTranscriptId = nowRawT.previousTranscriptId;
            rightAlignId = nowRawT.newAlignId;
        }

        t.exons[0].genomeStart -= rt.extendedLengthForward;
        t.exons[0].length += rt.extendedLengthForward;
        t.exons[0].readStart -= rt.extendedLengthForward;
        t.genomeStart = t.exons[0].genomeStart;
        t.readStart = t.exons[0].readStart;
        t.genomeStart = t.exons[0].genomeStart;
        return t;
    }


    void Stitching::stitchSingle(const Window &window) {
        int numAligns = window.numAligns;
        if (numAligns == 0) return;

        // try to extend each alignment in the window
        for (int i = 0; i < numAligns; ++i) {
            extendWindowAlign(window, window.aligns[i], extendRecords_[0][i], 0); // forward
            extendWindowAlign(window, window.aligns[i], extendRecords_[1][i], 1); // backward
        }

        // max score estimate
        {
            int64_t nowReadStart = 0;
            int64_t nowScore = 0;
            int64_t maxPossibleScore = 0;
            for (int i = 0; i < numAligns; ++i) {
                if (nowReadStart < window.aligns[i].readStart) {
                    maxPossibleScore += nowScore;
                    nowScore = 0;
                    nowReadStart = window.aligns[i].readStart;
                }
                nowScore = std::max(window.aligns[i].length + extendRecords_[0][i].maxExtensionScore +
                                    extendRecords_[1][i].maxExtensionScore, nowScore);
            }
            maxPossibleScore += nowScore;
            if (maxPossibleScore < outFilterScoreMin_) return; // cannot produce a valid transcript, skip this window
        }

        for (int i = 0; i < numAligns; ++i) {
            if (window.aligns[i].readStart + window.aligns[i].length >= read_->length - 5) {
                //todo replace the hardcoded 5 with a parameter
                for (int j = i + 1; j < numAligns; ++j) {
                    StitchRecord &record = stitchRecords_[calcStitchRecPos(i, j)];
                    record.type = StitchRecord::CANNOT_STITCH; // cannot stitch if one alignment reaches the end of the read
                }
            }
            for (int j = i + 1; j < numAligns; ++j) {
                StitchRecord &record = stitchRecords_[calcStitchRecPos(i, j)];
                stitchWindowAlign(window, window.aligns[i], window.aligns[j], record);
            }
        }

        // DP process
        // calculate the best stitching beginning with i and ending with j, for all alignments i,j
        for (int dist = 0; dist < numAligns; ++dist) {
            // calculate [i,i+dist]
            for (int i = 0; i < numAligns - dist; ++i) {
                RawTranscript &t = rawTranscripts_[calcRawTranscriptPos(i, dist, numAligns)];
                if (dist == 0) {
                    t.init(window, i, MATCH_SCORE_);
                } else {
                    t.score = transcriptMinScore;
                    t.StartAlignId = i;
                    for (int j = 0; j < dist; ++j) {
                        const RawTranscript &formerT = rawTranscripts_[calcRawTranscriptPos(i, j, numAligns)];
                        if (formerT.score < 0) continue;
                        const StitchRecord &sr = stitchRecords_[calcStitchRecPos(i + j, i + dist)];
                        if (sr.type == StitchRecord::CANNOT_STITCH) continue;
                        // filter short exon
                        // handle sjstrand consistency
                        if (sr.type == StitchRecord::SPLICE_JUNCTION) {
                            if (sr.isj == -1) {
                               //todo replace the hardcoded 5 with a parameter
                                if ( formerT.lastExonLength + sr.formerExonLengthShift - sr.shiftLeft < 5) continue; 
                                if ( window.aligns[i + dist].length + sr.latterExonLengthShift - sr.shiftRight < 5) continue;
                            }else {
                                //sjdb
                                if (formerT.lastExonLength +sr.formerExonLengthShift < 3) continue;
                                if ( window.aligns[i + dist].length + sr.latterExonLengthShift < 3) continue;
                            }

                            if (formerT.sjStrand != -1){
                                if (formerT.sjStrand != 2 - sr.spliceJunctionType % 2) continue;
                            }
                            
                        }

                        if (formerT.mismatches + sr.mismatches > maxMismatch_) continue;
                        if (formerT.numExon >= maxExons_) continue;

                        int nowScore = formerT.score + sr.score;
                        if (nowScore > t.score || (nowScore == t.score && formerT.numExon > t.numExon)) {
                            t.score = nowScore;
                            t.mismatches = formerT.mismatches + sr.mismatches;
                            t.matches = formerT.matches + sr.matches;
                            t.previousTranscriptId = calcRawTranscriptPos(i, j, numAligns);
                            t.newAlignId = i + dist;
                            if (sr.type == StitchRecord::PERFECT_MATCH || sr.type == StitchRecord::SAME_GAP)
                                t.numExon = formerT.numExon;
                            else t.numExon = formerT.numExon + 1;
                            if (sr.type == StitchRecord::SPLICE_JUNCTION) {
                                t.sjStrand = 2 - sr.spliceJunctionType % 2;
                            }else {
                                t.sjStrand = formerT.sjStrand;
                            }
                            t.lastExonLength = window.aligns[i + dist].length + sr.latterExonLengthShift;
                        }
                    }
                }
            }
        }

        //extend all best [i,j] transcript
        int firstDirToExtend = window.direction;
        int localBestScore = transcriptMinScore;

        for (int i = 0; i < numAligns; ++i) {
            for (int j = 0; j < numAligns - i; ++j) {
                RawTranscript &t = rawTranscripts_[calcRawTranscriptPos(i, j, numAligns)];
                if (t.score < 0) continue;
                int maxMismatchRemain = maxMismatch_ - t.mismatches;
                if (maxMismatchRemain < 0) {
                    t.score = transcriptMinScore;
                    continue;
                }

                int dir = firstDirToExtend;
                int endAlignId = dir == 0 ? i : i + j;
                const ExtendRecord &e1 = extendRecords_[dir][endAlignId];
                int extendLength = 0;

                const ExtendRecord::singleExtendRecord &extendInfo1 = e1.maxExtensionLengthWithMismatch[std::min(
                        e1.maxMismatch, maxMismatchRemain)];
                int score;

                maxMismatchRemain -= extendInfo1.mismatches;
                score = extendInfo1.matched * MATCH_SCORE_ + extendInfo1.mismatches * MISMATCH_PENALTY_;
                t.score += score;
                t.mismatches += extendInfo1.mismatches;
                t.matches += extendInfo1.matched;
                extendLength += extendInfo1.length;

                if (dir == 0) {
                    t.extendedLengthForward = extendLength;
                } else {
                    t.extendedLengthBackward = extendLength;
                }

                dir = 1 - dir;
                endAlignId = dir == 0 ? i : i + j;
                const ExtendRecord &e2 = extendRecords_[dir][endAlignId];
                const ExtendRecord::singleExtendRecord &extendInfo2 = e2.maxExtensionLengthWithMismatch[std::min(
                        e2.maxMismatch, maxMismatchRemain)];
                extendLength = 0;
                score = extendInfo2.matched * MATCH_SCORE_ + extendInfo2.mismatches * MISMATCH_PENALTY_;
                t.score += score;
                t.mismatches += extendInfo2.mismatches;
                t.matches += extendInfo2.matched;
                extendLength = extendInfo2.length;

                if (dir == 0) {
                    t.extendedLengthForward = extendLength;
                } else {
                    t.extendedLengthBackward = extendLength;
                }

                //finalize the transcript
                //add length penalty
                t.score += int(ceil(log2(double(window.aligns[i + j].genomeStart + window.aligns[i + j].length -
                                                window.aligns[i].genomeStart + t.extendedLengthForward + t.extendedLengthBackward)) * -0.25 - 0.5));

                if (t.score > localBestScore)
                    localBestScore = t.score;


            }
        }

        if (localBestScore > maxTranscriptScore_) maxTranscriptScore_ = localBestScore;
        else if (localBestScore < std::max(maxTranscriptScore_ - multimapScoreRange_, outFilterScoreMin_)) return;
        int scoreThreshold = std::max(localBestScore - multimapScoreRange_, outFilterScoreMin_);

        // add the good transcript to the result transcript set

        std::string chrName = genomeIndex_.genome_.chromosomes_[window.chrIndex].name;
        int64_t chrStart = genomeIndex_.genome_.chromosomes_[window.chrIndex].start;

        for (int i = 0; i < numAligns; ++i) {
            for (int j = 0; j < numAligns - i; ++j) {
                RawTranscript &t = rawTranscripts_[calcRawTranscriptPos(i, j, numAligns)];
                if (t.score < scoreThreshold) continue;
                if (t.score < outFilterScoreMin_) continue;
                if (t.matches < outFilterMatchMin_) continue;

                Transcript transcript = convertRawTranscriptToTranscript(t, window, i);
                transcript.chr = chrName;
                transcript.posInChr = transcript.exons[0].genomeStart - chrStart;
                transcript.readLength = read_->length;
                transcript.CIGAR = transcript.getCIGAR();

                //insert transcript
                bool alreadyExists = false;
                for (int iT = 0; iT < numGoodTranscripts_; ++iT) {
                    if (transcripts_[iT].genomeStart == transcript.genomeStart &&
                        transcripts_[iT].exons.back().genomeStart + transcripts_[iT].exons.back().length ==
                        transcript.exons.back().genomeStart + transcript.exons.back().length &&
                        transcript.CIGAR == transcripts_[iT].CIGAR) {

                        // same transcript already exists
                        alreadyExists = true;
                        break;
                    }
                }
                if (alreadyExists) continue;

                // find position to insert
                int posToInsert = numGoodTranscripts_;
                if (numGoodTranscripts_ < transcriptStoredMax_) {
                    ++numGoodTranscripts_;
                    for (posToInsert = numGoodTranscripts_ - 1; posToInsert > 0; --posToInsert) {
                        if (transcripts_[posToInsert - 1].score < t.score ||
                            (transcripts_[posToInsert - 1].score == t.score &&
                             transcripts_[posToInsert - 1].matched < t.matches)) {
                            transcripts_[posToInsert] = std::move(transcripts_[posToInsert - 1]);
                        } else {
                            break;
                        }
                    }
                } else {
                    if (transcripts_[posToInsert - 1].score < t.score ||
                        (transcripts_[posToInsert - 1].score == t.score &&
                         transcripts_[posToInsert - 1].matched < t.matches)) {
                        posToInsert = -1; // full, cannot insert
                    } else {
                        for (posToInsert = numGoodTranscripts_ - 1; posToInsert > 0; --posToInsert) {
                            if (transcripts_[posToInsert - 1].score < t.score ||
                                (transcripts_[posToInsert - 1].score == t.score &&
                                 transcripts_[posToInsert - 1].matched < t.matches)) {
                                transcripts_[posToInsert] = std::move(transcripts_[posToInsert - 1]);
                            } else {
                                break;
                            }
                        }
                    }
                }

                if (posToInsert != -1) {
                    Transcript &nowT = transcripts_[posToInsert];
                    nowT = transcript;
                }

            }
        }

    }

    inline int calcPairedRawTranscriptPos(int i, int dist, int numFrag0, int numFrag1, int iFrag) {
        if (iFrag == 0) {
            return i * numFrag0 + dist;
        } else {
            return numFrag0 * numFrag0 + i * numFrag1 + dist;
        }

    }

    void RawTranscriptPaired::init(const Window &win, int i, int matchScore) {
        //todo fix
    }

    void Stitching::stitchPaired(const Window &window) {
        // notice: only handle the case without overlapping
        int numAligns = window.numAligns;
        int numFrag0 = window.numFirstFragAligns;
        int numFrag1 = window.numSecondFragAligns;
        if (numAligns == 0) return;

        // try to extend each alignment in the window
        for (int i = 0; i < numAligns; ++i) {
            extendWindowAlign(window, window.aligns[i], extendRecords_[0][i], 0); // forward
            extendWindowAlign(window, window.aligns[i], extendRecords_[1][i], 1); // backward
        }

        // todo may implement some algorithm for paired-end
        // max score estimate
        {
            int64_t nowReadStart = 0;
            int64_t nowScore = 0;
            int64_t maxPossibleScore = 0;
            for (int i = 0; i < numAligns; ++i) {
                if (nowReadStart < window.aligns[i].readStart) {
                    maxPossibleScore += nowScore;
                    nowScore = 0;
                    nowReadStart = window.aligns[i].readStart;
                }
                nowScore = std::max(window.aligns[i].length + extendRecords_[0][i].maxExtensionScore +
                                    extendRecords_[1][i].maxExtensionScore, nowScore);
            }
            maxPossibleScore += nowScore;
            if (maxPossibleScore < outFilterScoreMin_) return; // cannot produce a valid transcript, skip this window
        }

        //stitch separately
        //todo better memory usage
        for (int i = 0; i < numFrag0; ++i) {
            for (int j = i + 1; j < numFrag0; ++j) {
                StitchRecord &record = stitchRecords_[calcStitchRecPos(i, j)];
                stitchWindowAlign(window, window.aligns[j], window.aligns[i], record);
            }
        }

        for (int i = numFrag0; i < numAligns; ++i) {
            for (int j = i + 1; j < numAligns; ++j) {
                StitchRecord &record = stitchRecords_[calcStitchRecPos(i, j)];
                stitchWindowAlign(window, window.aligns[j], window.aligns[i], record);
            }
        }

        //local DP of each fragment
        for (int dist = 0; dist < numFrag0; ++dist) {
            for (int i = 0; i < numFrag0; ++i) {
                RawTranscriptPaired &t = rawTranscriptsPaired_[calcPairedRawTranscriptPos(i, dist, numFrag0, numFrag1,
                                                                                          0)];
                if (dist == 0) {
                    t.init(window, i, MATCH_SCORE_);
                    t.iFragment = 0;
                } else {
                    t.score = transcriptMinScore;
                    t.StartAlignId = i;
                    for (int j = 0; j < dist; ++j) {
                        const RawTranscriptPaired &formerT = rawTranscriptsPaired_[calcPairedRawTranscriptPos(i, j,
                                                                                                              numFrag0,
                                                                                                              numFrag1,
                                                                                                              0)];
                        if (formerT.score < 0) continue;
                        const StitchRecord &sr = stitchRecords_[calcStitchRecPos(i + j, i + dist)];
                        if (sr.type == StitchRecord::CANNOT_STITCH) continue;

                        if (formerT.mismatches + sr.mismatches > maxMismatch_) continue;
                        if (formerT.numExon >= maxExons_) continue;

                        int nowScore = formerT.score + sr.score;
                        if (nowScore > t.score || (nowScore == t.score && formerT.numExon > t.numExon)) {
                            t.score = nowScore;
                            t.mismatches = formerT.mismatches + sr.mismatches;
                            t.matches = formerT.matches + sr.matches;
                            t.previousTranscriptId = calcRawTranscriptPos(i, j, numAligns);
                            t.newAlignId = i + dist;
                            if (sr.type == StitchRecord::PERFECT_MATCH || sr.type == StitchRecord::SAME_GAP)
                                t.numExon = formerT.numExon;
                            else t.numExon = formerT.numExon + 1;
                        }
                    }
                }
            }
        }

        for (int dist = 0; dist < numFrag1; ++dist) {
            for (int i = numFrag0; i < numAligns; ++i) {
                RawTranscriptPaired &t = rawTranscriptsPaired_[calcPairedRawTranscriptPos(i, dist, numFrag0, numFrag1,
                                                                                          0)];
                if (dist == 0) {
                    t.init(window, i, MATCH_SCORE_);
                    t.iFragment = 0;
                } else {
                    t.score = transcriptMinScore;
                    t.StartAlignId = i;
                    for (int j = 0; j < dist; ++j) {
                        const RawTranscriptPaired &formerT = rawTranscriptsPaired_[calcPairedRawTranscriptPos(i, j,
                                                                                                              numFrag0,
                                                                                                              numFrag1,
                                                                                                              0)];
                        if (formerT.score < 0) continue;
                        const StitchRecord &sr = stitchRecords_[calcStitchRecPos(i + j, i + dist)];
                        if (sr.type == StitchRecord::CANNOT_STITCH) continue;

                        if (formerT.mismatches + sr.mismatches > maxMismatch_) continue;
                        if (formerT.numExon >= maxExons_) continue;

                        int nowScore = formerT.score + sr.score;
                        if (nowScore > t.score || (nowScore == t.score && formerT.numExon > t.numExon)) {
                            t.score = nowScore;
                            t.mismatches = formerT.mismatches + sr.mismatches;
                            t.matches = formerT.matches + sr.matches;
                            t.previousTranscriptId = calcRawTranscriptPos(i, j, numAligns);
                            t.newAlignId = i + dist;
                            if (sr.type == StitchRecord::PERFECT_MATCH || sr.type == StitchRecord::SAME_GAP)
                                t.numExon = formerT.numExon;
                            else t.numExon = formerT.numExon + 1;
                        }
                    }
                }
            }
        }

        //try to combine fragments
        //simple enumeration
        /*for (int i0 = 0; i0 < numFrag0; ++i0){
            for (int j1 = numFrag1; j1 <numAligns; ++j1){
                for (int j0 = i0; j0 < numFrag0; ++j0){
                    for (int i1 = numFrag1; i1 < numAligns; ++i1){
                        //try stitch [i0,j0] and [i1,j1]
                        RawTranscriptPaired &t0 = rawTranscriptsPaired_[calcPairedRawTranscriptPos(i0,j0 - i0,numFrag0,numFrag1,0)];
                        RawTranscriptPaired &t1 = rawTranscriptsPaired_[calcPairedRawTranscriptPos(i1,j1 - i1,numFrag0,numFrag1,1)];
                        int maxMismatchRemain = maxMismatch_ - t0.mismatches - t1.mismatches;
                        if (maxMismatchRemain < 0) continue;
                        if (!examineSpliceJunction(window,t0,t1)) continue;
                        //extend and finalize
                        //todo notice the extension direction!!!
                        const ExtendRecord &e00 = extendRecords_[0][i0];
                        const ExtendRecord &e01 = extendRecords_[1][j0];
                        const ExtendRecord &e10 = extendRecords_[0][i1];
                        const ExtendRecord &e11 = extendRecords_[1][j1];



                    }
                }
            }
        }*/


        /*for (int j0 = 0; j0 < numFrag0; ++j0){
            for (int i1 = numFrag0; i1 <numAligns; ++i1){
                // detect overlap
                const WindowAlign &a0 = window.aligns[j0]; // the last align of fragment 0
                const WindowAlign &a1 = window.aligns[i1]; // the first align of fragment 1
                int endAlign0 = a0.genomeStart + a0.length;
                if (endAlign0 >= a1.genomeStart){
                    // overlap, handle splice junction compatibility
                    // find the best compatible pair covering the overlap

                    // first, find a window align in fragment0 that overlaps with a1
                    // only need to check those before j0
                    for (int k = j0; k >= 0; --k){

                    }
                }

            }
        }*/



        /*for (int i0 = 0; i0 < numFrag0; ++i0){
            for (int j0 = i0; j0 < numFrag0; ++j0){
                //Fragment 0 aligns begin with i0 and end with j0
                const RawTranscriptPaired &t0 = rawTranscriptsPaired_[calcPairedRawTranscriptPos(i0,j0 - i0,numFrag0,numFrag1,0)];
                if (t0.score < 0) continue; // no valid transcript
                int endAlignI0 = window.aligns[i0].genomeStart + window.aligns[j0].length;
                int endAlignJ0 = window.aligns[j0].genomeStart + window.aligns[j0].length;
                // now find a compatible transcript in fragment 1
                for (int i1 = numFrag1; i1 < numAligns; ++i1){
                    int startAlignI1 = window.aligns[i1].genomeStart;
                    int endAlignI1 = window.aligns[i1].genomeStart + window.aligns[i1].length;
                    if (endAlignI1 < window.aligns[i0].genomeStart) continue; // too left

                }
            }
        }*/

        /*for (int j0 = 0;j0 < numFrag0; ++j0){
            for (int i1 = numFrag1;i1 <numAligns; ++i1){
                // find the best compatible pair if overlap

                const WindowAlign &a0 = window.aligns[j0]; // the last align of fragment 0
                const WindowAlign &a1 = window.aligns[i1]; // the first align of fragment 1
                int startGenomePos = a1.genomeStart;
                int endGenomePos = a0.genomeStart + a0.length;
                if (endGenomePos <= startGenomePos) continue; // no overlap, skip
            }
        }*/






    }


}