#include "StitchingManagement.h"
#include "../utils/exceptions.h"
#include "../utils/types.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include <chrono>


namespace rna {
    inline int64_t log2(int64_t n) {
        if (n <= 0) return -1; // invalid input
        int64_t ans = 0;
        while (n > 0) {
            n >>= 1;
            ans++;
        }
        return ans - 1; // log2(n) = ans - 1
    }

    StitchingManagement::StitchingManagement(const stitchingConfig &config,
                                             const GenomeIndexPrefix &genomeIndex)
            : config_(config), genomeIndex_(genomeIndex) {
        winBinMap_ = new int32_t[(genomeIndex_.suffixArray.fullLength_ >> (config.winBinSizeLog - 1)) + 2];
        for (size_t i = 0; i < (genomeIndex_.suffixArray.fullLength_ >> (config.winBinSizeLog - 1)) + 2; ++i) {
            winBinMap_[i] = -1; // initialize all bins to -1
        }
    }

    StitchingManagement::~StitchingManagement() {
        delete[] winBinMap_;
    }

    void StitchingManagement::processAlignments(std::vector<Align> &alignments, ReadPtr read) {

        if (alignments.empty()) {
            throw RNAException(
                    ExceptionType::INVALID_INPUT,
                    "No alignments provided for stitching"
            );
        }


        bestTranscript_ = std::make_shared<Transcript>();
        read_ = std::move(read);
        windows_.reserve(config_.maxWindows);
        transcripts_.reserve(100);

        identifyAnchors(alignments);

        createWindows(alignments);

        assignAlignmentsToWindows(alignments);

        generateTranscripts();

    }


    void StitchingManagement::clear() {
        //clear previous status
        windows_.clear();
        for (size_t i = 0; i < (genomeIndex_.suffixArray.fullLength_ >> (config_.winBinSizeLog - 1)) + 2; ++i) {
            winBinMap_[i] = -1; // initialize all bins to -1
        }
        transcripts_.clear();
        sjdb.clear();
        bestTranscript_.reset();
        totalAnchors_ = 0;
    }


    void StitchingManagement::identifyAnchors(std::vector<Align> &alignments) {
        for (auto &alignment: alignments) {
            // Check whether it is an anchor
            if (alignment.rep <= config_.maxAnchorRep) {
                alignment.isAnchor = true;
                totalAnchors_++;
            } else {
                alignment.isAnchor = false;
            }
        }
    }

    inline StitchingManagement::PositiveStrandAlign
    StitchingManagement::convertAlignToPositiveStrand(const Align &align, size_t SAi) {
        int64_t location = genomeIndex_.suffixArray[SAi];
        if (location > genomeIndex_.genomeLength) {
            // reverse strand alignment
            int64_t genomeStart = genomeIndex_.genomeLength * 2 - location - align.length;
            int64_t readStart = read_->length - align.readStart - align.length;
            return {genomeStart, readStart,align.length, 1 - align.direction,align.isAnchor};
        } else {
            // positive strand alignment
            return {location, align.readStart,align.length, align.direction,align.isAnchor};
        }
    }

    void StitchingManagement::createWindows(const std::vector<Align> &alignments) {
        windows_.reserve(config_.maxWindows);
        for (const auto &align: alignments) {
            if (!align.isAnchor) continue;

            for (size_t i = align.leftSAIndex; i <= align.rightSAIndex; ++i) {
                // Create Windows
                auto positiveAlign = convertAlignToPositiveStrand(align, i);
                auto location = positiveAlign.genomeStart;
                auto chrId = genomeIndex_.genome->getPosChrIndex(location);
                Genome::Chromosome chr = genomeIndex_.genome->chromosomes_[chrId];
                int64_t chrStartBin = chr.start >> config_.winBinSizeLog;
                int64_t chrEndBin = (chr.start + chr.length - 1) >> config_.winBinSizeLog;

                Window window;
                window.chrIndex = chrId;
                window.direction = positiveAlign.direction;
                int64_t baseBin = location >> config_.winBinSizeLog;
                int64_t baseBinKey = (baseBin << 1) | window.direction;

                if (winBinMap_[baseBinKey] != -1) {
                    // already exists, skip
                    continue;
                }


                // try merging with left existing windows
                int64_t binKey = baseBinKey;
                int64_t leftBound = (std::max(chrStartBin,baseBin - config_.winAnchorDistBins) << 1)| window.direction;
                int32_t leftOverlap = -1;
                for (binKey = binKey - 2; binKey >= leftBound; binKey -= 2){
                    if (winBinMap_[binKey] != -1) {
                        leftOverlap = winBinMap_[binKey];
                        break;
                    }
                }
                int32_t nowWindowIndex;
                if (leftOverlap != -1) {
                    // merge with existing window
                    nowWindowIndex = leftOverlap;
                    for (binKey = binKey + 2;binKey <= baseBinKey; binKey += 2) {
                        winBinMap_[binKey] = nowWindowIndex;
                    }

                } else {
                    //create a new window
                    nowWindowIndex = windows_.size();
                }
                winBinMap_[baseBinKey] = nowWindowIndex;

                // try merging with right existing windows
                int64_t rightBound = (std::min(chrEndBin, baseBin + config_.winAnchorDistBins) << 1) | window.direction;
                int32_t rightOverlap = -1;
                for (binKey = baseBinKey + 2; binKey <= rightBound; binKey += 2) {
                    if (winBinMap_[binKey] != -1) {
                        rightOverlap = winBinMap_[binKey];
                    }
                }
                if (rightOverlap != -1) {
                    binKey = baseBinKey + 2;
                    //extend to reach the overlapping window
                    while (winBinMap_[binKey] != rightOverlap){
                        winBinMap_[binKey] = nowWindowIndex;
                        binKey += 2;
                    }
                    // merge with right window
                    while (winBinMap_[binKey] == rightOverlap){
                        winBinMap_[binKey] = nowWindowIndex;
                        binKey += 2;
                    }

                    window.endBin = binKey >> 1;
                    // kill right window
                    windows_[rightOverlap].startBin = 1;
                    windows_[rightOverlap].endBin = 0;
                }else{
                    window.endBin = baseBin;
                }
                window.startBin = baseBin;
                // update window information
                if (leftOverlap == -1){
                    //new window
                    windows_.emplace_back(window);
                }else{
                    windows_[leftOverlap].endBin = window.endBin;
                }

            }
        }


        // flank existing windows
        for (int32_t i = 0; i < windows_.size();++i){
            auto &win = windows_[i];
            if(win.startBin > win.endBin) continue;
            // flank the window
            auto &chr = genomeIndex_.genome->chromosomes_[win.chrIndex];
            int64_t chrStartBin = chr.start >> config_.winBinSizeLog;
            int64_t chrEndBin = (chr.start + chr.length - 1) >> config_.winBinSizeLog;
            int64_t leftBin = std::max(chrStartBin, win.startBin - config_.flankSize);
            int64_t rightBin = std::min(chrEndBin, win.endBin + config_.flankSize);
            int64_t leftKeyBound = (win.startBin << 1)|win.direction;
            for(int64_t binKey = (leftBin << 1)| win.direction; binKey < leftKeyBound;binKey += 2){
                winBinMap_[binKey] = i;
            }
            int64_t rightKeyBound = (rightBin << 1)|win.direction;
            for(int64_t binKey = (win.endBin << 1)| win.direction; binKey < rightKeyBound;binKey += 2){
                winBinMap_[binKey] = i;
            }
            win.startBin = leftBin;
            win.endBin = rightBin;
            // reserve space for alignments
            win.aligns.reserve(config_.maxSeedPerWindows);
        }


    }

    void StitchingManagement::createChrWindows(const std::vector<Align> &alignments) {
        //create windows for each chromosome
        for (int i = 0; i < genomeIndex_.genome->chromosomes_.size(); ++i) {
            const auto &Chr = genomeIndex_.genome->chromosomes_[i];
            // +
            Window window;
            window.chrIndex = i;
            window.direction = 0; // default direction

            window.aligns.reserve(50);
            windows_.emplace_back(window);
            // -
            Window windowReverse;
            windowReverse.chrIndex = i;
            windowReverse.direction = 1; // reverse direction
            windowReverse.aligns.reserve(50);
            windows_.emplace_back(windowReverse);
        }
    }

    void StitchingManagement::assignSingleAlignment(Window& win,const PositiveStrandAlign &a){

        // todo all windows aligns can be stored in a huge array in StitchingManagement, only store pointer and alignNum in Windows.
        // when window is full, check if it needs to be replaced
        if(a.length < win.minLengthWhenFull && !a.isAnchor) return; // ignore too short no anchor alignment

        //detect overlap
        std::vector<WindowAlign> &aligns = win.aligns;
        int alignNum = aligns.size();
        for(int i = 0; i < alignNum; ++i) {
            if (a.genomeStart + aligns[i].readStart == aligns[i].genomeStart + a.readStart \
                && ((a.readStart >= aligns[i].readStart ) && a.readStart < aligns[i].readStart + aligns[i].length)\
                ||(a.readStart + a.length > aligns[i].readStart && a.readStart + a.length <= aligns[i].readStart + aligns[i].length)) {
                // overlap

                // same
                if (a.genomeStart == aligns[i].genomeStart && a.length == aligns[i].length) return;


                if (a.length > aligns[i].length){
                    //delete the old alignment and insert the new one
                    //keep the array sorted by readStart

                    //find the position to insert
                    //only need to search one side
                    if (a.readStart > aligns[i].readStart){
                        for (int j = i+1;j<alignNum;++j){
                            if(aligns[j].readStart > a.readStart){
                                aligns[j-1] = WindowAlign{
                                    a.readStart,
                                    a.genomeStart,
                                    a.length,
                                    a.length + log2(a.length),
                                    a.isAnchor
                                };
                                break;
                            }else aligns[j-1] = aligns[j];
                        }
                    }else {
                        for (int j = i-1; j>=0;--j){
                            if(aligns[j].readStart <= a.readStart){
                                aligns[j+1] = WindowAlign{
                                    a.readStart,
                                    a.genomeStart,
                                    a.length,
                                    a.length + log2(a.length),
                                    a.isAnchor
                                };
                                break;
                            }else aligns[j+1] = aligns[j];
                        }
                        // reaching here means that a.readStart is the smallest, insert at the beginning
                        aligns[0] = WindowAlign{
                            a.readStart,
                            a.genomeStart,
                            a.length,
                            a.length + log2(a.length),
                            a.isAnchor
                        };
                    }
                }
                // if not, do nothing
                return;
            }
        }

        if(a.isAnchor) ++ win.numAnchors; // anchor must be added to the window

        // handle the case that there are too many seeds in the window
        if (alignNum == config_.maxSeedPerWindows){
            // calculate minLengthWhenFull and the alignment to remove
            win.minLengthWhenFull = std::numeric_limits<int>::max();
            int removePos = -1;
            for (int i = 0; i < alignNum; ++i){
                if (!aligns[i].isAnchor && aligns[i].length < win.minLengthWhenFull){
                    win.minLengthWhenFull = aligns[i].length;
                    removePos = i;
                }
            }
            if (removePos == -1) return; // all are anchors, cannot add new alignment
            if (win.minLengthWhenFull >= a.length && !a.isAnchor) return; // new alignment is too short, ignore

            //remove the removePos alignment and add the new one

            for (int i = 0; i <= removePos; ++i){
                if (aligns[i].readStart > a.readStart){
                    // add here
                    // move [i,removePos-1] to [i+1,removePos]
                    for (int j = i+1;j<= removePos;++j) aligns[j] = aligns[j-1];
                    aligns[i] = WindowAlign{
                        a.readStart,
                        a.genomeStart,
                        a.length,
                        a.length + log2(a.length),
                        a.isAnchor
                    };
                    return;
                }
            }
            // reaching here means that a.readStart > alignToBeRemoved.readStart
            for (int i = removePos;i<alignNum-1;++i) {
                if (aligns[i+1].readStart > a.readStart) {
                    // add to i
                    aligns[i] = WindowAlign{
                        a.readStart,
                        a.genomeStart,
                        a.length,
                        a.length + log2(a.length),
                        a.isAnchor
                    };
                    return;
                }else aligns[i] = aligns[i+1];
            }
            // reaching here means that a.readStart is the largest
            aligns[alignNum-1] = WindowAlign{
                a.readStart,
                a.genomeStart,
                a.length,
                a.length + log2(a.length),
                a.isAnchor
            };
            return;
        }


        // simple case
        // find the position to insert
        aligns.emplace_back(WindowAlign()); // add a temp element to increase vector size, to be removed in array implementation
        for (int i = alignNum; i>0; --i){
            if (aligns[i-1].readStart <= a.readStart){
                aligns[i] = WindowAlign{
                    a.readStart,
                    a.genomeStart,
                    a.length,
                    a.length + log2(a.length),
                    a.isAnchor
                };
                return;
            }else aligns[i] = aligns[i-1];
        }
        // reaching here means that a.readStart is the smallest, insert at the beginning
        aligns[0] = WindowAlign{
            a.readStart,
            a.genomeStart,
            a.length,
            a.length + log2(a.length),
            a.isAnchor
        };
    }

    void StitchingManagement::assignAlignmentsToWindows(const std::vector<Align> &alignments) {
        if (genomeIndex_.config.twoDirections) {
            for (const auto &align: alignments) {
                if (align.rep > config_.maxRep) continue; // skip alignments with too many reps

                for (size_t i = align.leftSAIndex; i <= align.rightSAIndex; ++i) {

                    auto a = convertAlignToPositiveStrand(align, i);
                    size_t binKey = ((a.genomeStart >> config_.winBinSizeLog) << 1) |a.direction;
                    auto winId = winBinMap_[binKey];
                    if (winId == -1) continue;
                    // find the corresponding window
                    Window &win = windows_[winId];
                    assignSingleAlignment(win,a);
                }
            }
        } else {
            std::cout<<"Error: not implemented"<<std::endl;
            exit(1);
        }

        // no need to sort alignments in each window, as they are inserted in order

    }


    void StitchingManagement::assignAlignmentsToChrWindows(const std::vector<Align> &alignments) {
        if (genomeIndex_.config.twoDirections) {
            for (const auto &align: alignments) {
                if (align.rep > 10000) continue; // skip too many reps
                for (size_t i = align.leftSAIndex; i <= align.rightSAIndex; ++i) {
                    int64_t readStart = align.readStart;
                    int direction = align.direction;
                    GenomePos location = genomeIndex_.suffixArray[i];
                    if (location > genomeIndex_.genomeLength) {
                        // reverse strand alignment
                        direction = 1 - align.direction;
                        location = genomeIndex_.genomeLength * 2 - location - align.length;
                        readStart = read_->length - align.readStart - align.length;

                    }

                    int64_t chrIndex = genomeIndex_.genome->getPosChrIndex(location);
                    Window &win = windows_[(chrIndex << 1) | direction];
                    win.aligns.emplace_back(WindowAlign{
                            readStart,
                            location,
                            align.length,
                            MATCH_SCORE * align.length + log2(align.length),
                            align.isAnchor
                    });
                    if (align.isAnchor) {
                        win.numAnchors++;
                    }
                }
            }
        } else {
            for (const auto &align: alignments) {
                if (align.length < config_.minAlignLength) continue; // skip too short alignments
                for (size_t i = align.leftSAIndex; i <= align.rightSAIndex; ++i) {
                    GenomePos location = genomeIndex_.suffixArray[i];
                    //GenomePos location = 0;//for test
                    int64_t chrIndex = genomeIndex_.genome->getPosChrIndex(location);
                    Window &win = windows_[(chrIndex << 1) | align.direction];
                    win.aligns.emplace_back(WindowAlign{
                            align.readStart,
                            location,
                            align.length,
                            MATCH_SCORE * align.length + log2(align.length),
                            align.isAnchor
                    });
                    if (align.isAnchor) {
                        win.numAnchors++;
                    }
                }
            }
        }

        // sort alignments in each window
        for (auto &window: windows_) {
            if (!window.aligns.empty()) std::sort(window.aligns.begin(), window.aligns.end());
        }
    }

    // Generate transcripts for each window and find the best transcript
    void StitchingManagement::generateTranscripts() {

        for (auto &window: windows_) {

            stitchWindowAligns(window);

        }

        for (auto t: transcripts_) {
            bestTranscript_ = bestTranscript_->score < t.score ? std::make_shared<Transcript>(t) : bestTranscript_;
        }
    }


    void StitchingManagement::stitchWindowAligns(Window &window) {
        //simple DP
        size_t nAligns = window.aligns.size();
        if (nAligns == 0) return;  // no alignment
        if (window.numAnchors == 0) return; // no anchor, cannot stitch

        Transcript bestTranscriptOfCurrentWindow = Transcript();
        std::vector<Transcript> formerTranscripts;
        formerTranscripts.resize(nAligns);
        for (size_t i = 0; i < nAligns; ++i) {
            formerTranscripts[i].strand = window.direction;
        }
        // can be optimized by storing and reusing the 5' and 3' extension results
        for (size_t i = 0; i < nAligns; ++i) {
            Transcript bestTranscriptWithI = Transcript();
            bestTranscriptWithI.strand = window.direction;
            int64_t score = 0;// stitch score
            score = stitchAlignmentToTranscript(window.aligns[i], bestTranscriptWithI);
            for (size_t j = 0; j < i; ++j) {
                Transcript t = formerTranscripts[j];//copy
                score = stitchAlignmentToTranscript(window.aligns[i], t);
                if (score < -1000000) continue; // invalid stitch
                bestTranscriptWithI = bestTranscriptWithI.score < t.score ? t : bestTranscriptWithI;
            }
            formerTranscripts[i] = bestTranscriptWithI;
            finalizeTranscript(bestTranscriptWithI, window.chrIndex);
            bestTranscriptOfCurrentWindow = bestTranscriptWithI.score < bestTranscriptOfCurrentWindow.score ?
                                            bestTranscriptOfCurrentWindow : bestTranscriptWithI;
        }
        if (bestTranscriptOfCurrentWindow.exons.empty()) return;
        transcripts_.emplace_back(bestTranscriptOfCurrentWindow);
    }

    //return: {matched,unmatched}
    inline std::pair<int64_t, int64_t> compareGenomeRead(const std::string &genomeSeq,
                                                         const std::string &readSeq,
                                                         int64_t genomeStart,
                                                         int64_t readStart,
                                                         int64_t length) {
        if (genomeStart < 0 || readStart < 0 || length <= 0) {
            return {0, 0}; // invalid parameters, return zero matches
        }
        int64_t matched = 0;
        int64_t unmatched = 0;

        //If there are Ns in the same pos of both read and genome, they will not be counted as matches.

        int64_t i;
        for (i = 0; i < length; i++) {
            if (genomeSeq[genomeStart + i] == readSeq[readStart + i]) {
                matched++;
            } else {
                if (genomeSeq[genomeStart + i] != 'N' && readSeq[readStart + i] != 'N')
                    unmatched++;
            }

        }

        return {matched, unmatched};
    }

    //return {junctionType, penaltyScore}
    inline std::pair<int, int64_t> StitchingManagement::checkJunctionMotif(const std::string &genomeSeq,
                                                                           int64_t leftPos,
                                                                           int64_t rightPos) {
        // may optimize this function later, use hash and array
        if (genomeSeq[leftPos + 1] == 'G' && genomeSeq[leftPos + 2] == 'T') {
            if (genomeSeq[rightPos - 1] == 'A' && genomeSeq[rightPos] == 'G')
                return {1, 0};
            else if (genomeSeq[rightPos - 1] == 'A' && genomeSeq[rightPos] == 'T')
                return {6, SCORE_GAP_ATAC};
        } else if (genomeSeq[leftPos + 1] == 'C' && genomeSeq[leftPos + 2] == 'T') {
            if (genomeSeq[rightPos - 1] == 'A' && genomeSeq[rightPos] == 'C')
                return {2, 0};
            else if (genomeSeq[rightPos - 1] == 'G' && genomeSeq[rightPos] == 'C')
                return {4, SCORE_GAP_GCAG};
        } else if (genomeSeq[leftPos + 1] == 'A' && genomeSeq[leftPos + 2] == 'T') {
            if (genomeSeq[rightPos - 1] == 'A' && genomeSeq[rightPos] == 'C')
                return {5, SCORE_GAP_ATAC};
        } else if (genomeSeq[leftPos + 1] == 'G' && genomeSeq[leftPos + 2] == 'C') {
            if (genomeSeq[rightPos - 1] == 'A' && genomeSeq[rightPos] == 'G')
                return {3, SCORE_GAP_GCAG};
        }
        return {0, SCORE_GAP_NON_CANONICAL};
    }


    // stitch a new piece of alignment to the transcript and return the STITCH score. The alignment's match score is not included
    // however, when maintaining the transcript, the alignment's match score is added to the transcript's total score
    int64_t StitchingManagement::stitchAlignmentToTranscript(
            WindowAlign &alignment,
            Transcript &transcript) {
        std::string &rSeq = read_->sequence[transcript.strand];
        if (transcript.exons.empty()) {
            // add the first alignment
            transcript.readStart = alignment.readStart;
            transcript.genomeStart = alignment.genomeStart;
            transcript.exons.push_back(Exon{
                    alignment.genomeStart,
                    alignment.length,
                    alignment.readStart,
                    alignment.score
            });
            transcript.score = alignment.score;
            transcript.matched = alignment.length;
            transcript.aligns.push_back(alignment);

            return 0; // no stitching happened
        }

        int64_t lastReadPos = transcript.exons.back().readStart + transcript.exons.back().length - 1;
        int64_t lastGenomePos = transcript.exons.back().start + transcript.exons.back().length - 1;
        int64_t alignmentReadEnd = alignment.readStart + alignment.length - 1;
        int64_t alignmentReadStart = alignment.readStart;
        int64_t alignmentGenomeEnd = alignment.genomeStart + alignment.length - 1;
        int64_t alignmentGenomeStart = alignment.genomeStart;
        int64_t length = alignment.length;
        int64_t junctionType = 0;
        int64_t junctionPos;
        int64_t ins = 0;
        int64_t del = 0;
        // check if the alignment can be stitched
        if (alignmentReadEnd <= lastReadPos) {
            // fully overlap, cannot stitch
            return -1000001;
        }
        if (alignmentGenomeEnd <= lastGenomePos) {
            // fully overlap, cannot stitch
            return -1000002;
        }

        if (alignmentReadStart <= lastReadPos) {
            alignmentGenomeStart += (lastReadPos - alignmentReadStart + 1);
            alignmentReadStart = lastReadPos + 1;
            length = alignmentReadEnd - alignmentReadStart + 1;
        }

        // calculate gap
        int64_t readGap = alignmentReadStart - (lastReadPos + 1);
        int64_t genomeGap = alignmentGenomeStart - (lastGenomePos + 1);
        int64_t lastIntronBase = alignmentGenomeStart - readGap - 1; // if all gap belongs to acceptor
        int64_t score = 0;
        int64_t nMatch = 0, nMismatch = 0;


        if (readGap == genomeGap) {
            // perfect match
            if (readGap == 0) score = alignment.score;
            else {
                // same length gap, maybe there are some mismatches
                auto matchResult = compareGenomeRead(
                        genomeIndex_.genome->sequence_,
                        rSeq,
                        lastGenomePos + 1,
                        lastReadPos + 1,
                        readGap
                );
                score = matchResult.first * MATCH_SCORE + matchResult.second * MISMATCH_PENALTY;
                transcript.matched += matchResult.first;
                transcript.unmatched += matchResult.second;
            }
        } else if (genomeGap > readGap) {
            // genome gap larger than read gap
            // maybe there is a deletion
            // an intron may be present
            del = genomeGap - readGap;
            if (del > MAX_INTRON_LENGTH) {
                // gaps too large, cannot stitch
                return -1000003;
            }


            int64_t junctionScore = 0;
            junctionPos = 1;
            int junctionCan = 0;

            //find the best left pos
            do {
                junctionPos--;
                if (rSeq[lastReadPos + junctionPos] !=
                    genomeIndex_.genome->sequence_[lastIntronBase + junctionPos] &&
                    rSeq[lastReadPos + junctionPos] ==
                    genomeIndex_.genome->sequence_[lastGenomePos + junctionPos]) {
                    junctionScore += MISMATCH_PENALTY;
                }
            } while (junctionScore + SCORE_STITCH_SJ_SHIFT >= 0 &&
                     transcript.exons.back().length + junctionPos > 1);

            junctionScore = 0;
            int64_t maxScore = -1000000;
            int64_t maxScoreJunctionReadPos = 0;
            int64_t maxScoreJunctionCan = 0;
            int64_t maxScorePenalty = 0;
            do {
                if (rSeq[lastReadPos + junctionPos] !=
                    genomeIndex_.genome->sequence_[lastIntronBase + junctionPos] &&
                    rSeq[lastReadPos + junctionPos] ==
                    genomeIndex_.genome->sequence_[lastGenomePos + junctionPos]) {
                    junctionScore += MATCH_SCORE;
                }
                if (rSeq[lastReadPos + junctionPos] ==
                    genomeIndex_.genome->sequence_[lastIntronBase + junctionPos] &&
                    rSeq[lastReadPos + junctionPos] !=
                    genomeIndex_.genome->sequence_[lastGenomePos + junctionPos]) {
                    junctionScore += MISMATCH_PENALTY;
                }
                auto [jType, junctionPenalty] = checkJunctionMotif(
                        genomeIndex_.genome->sequence_,
                        lastGenomePos + junctionPos,
                        lastIntronBase + junctionPos
                );
                if (junctionScore + junctionPenalty > maxScore) {
                    maxScore = junctionScore + junctionPenalty;
                    maxScoreJunctionReadPos = junctionPos;
                    maxScoreJunctionCan = jType;
                    maxScorePenalty = junctionPenalty;
                }
                junctionPos++;
            } while (junctionPos < alignmentReadEnd - lastReadPos);
            junctionPos = maxScoreJunctionReadPos;
            // score donor acceptor
            for (int64_t i = std::min((int64_t) 1, maxScoreJunctionReadPos + 1);
                 i <= std::max(readGap, maxScoreJunctionReadPos); ++i) {
                uint64_t g = (i < maxScoreJunctionReadPos) ? (lastGenomePos + i) : (lastIntronBase + i);
                if (rSeq[lastReadPos + i] != genomeIndex_.genome->sequence_[g]) {
                    nMismatch++;
                    score += MISMATCH_PENALTY;
                    if (i < 1 || i > readGap) {
                        score -= MATCH_SCORE;
                        nMatch--;
                    }
                } else {
                    if (i >= 1 && i <= readGap) {
                        nMatch++;
                        score += MATCH_SCORE;
                    }
                }
            }

            // score the gap
            if (del > MIN_INTRON_LENGTH) {
                //intron
                score += GAP_OPEN_PENALTY + maxScorePenalty;

            } else {
                maxScoreJunctionCan = -1;//mark deletion
                score += DEL_OPEN_PENALTY + del * DEL_EXTEND_PENALTY;
            }
            junctionType = maxScoreJunctionCan;
            transcript.matched += nMatch;
            transcript.unmatched += nMismatch;

        } else if (readGap > genomeGap) {
            // insertion
            ins = readGap - genomeGap;
            junctionPos = 0;
            // when genomeGap == 0, no need to stitch
            if (genomeGap < 0) {
                //overlapping
                score += genomeGap * MATCH_SCORE;//reduce score
            } else if (genomeGap > 0) {
                int64_t tmpScore = 0, maxScore = 0;
                for (int64_t jPos = 1; jPos <= genomeGap; ++jPos) {
                    tmpScore += (rSeq[lastReadPos + jPos] ==
                                 genomeIndex_.genome->sequence_[lastGenomePos + jPos]) ?
                                MATCH_SCORE : -MATCH_SCORE;
                    tmpScore += (rSeq[lastReadPos + ins + jPos] ==
                                 genomeIndex_.genome->sequence_[lastGenomePos + jPos]) ?
                                -MATCH_SCORE : MATCH_SCORE;
                    if (tmpScore > maxScore) {
                        maxScore = tmpScore;
                        junctionPos = jPos;
                    }
                }

                //score donor and acceptor
                for (int64_t i = 1; i <= genomeGap; i++) {
                    int64_t r = lastReadPos + i + (i <= junctionPos ? 0 : ins);
                    if (genomeIndex_.genome->sequence_[lastGenomePos + i] !=
                        rSeq[r]) {
                        nMismatch++;
                        score += MISMATCH_PENALTY;
                    } else {
                        nMatch++;
                        score += MATCH_SCORE;
                    }
                }
            }
            score += ins * INS_EXTEND_PENALTY + INS_OPEN_PENALTY;
            transcript.matched += nMatch;
            transcript.unmatched += nMismatch;
            junctionType = -2; // mark insertion
        }

        // add to transcript, if there are no insertions or deletions, just extend the last exon
        if (genomeGap == readGap) {
            transcript.aligns.back().length += length + readGap;
            transcript.aligns.back().score += alignment.score;
            transcript.exons.back().length += length + readGap;
            transcript.exons.back().score += alignment.score;

        } else {
            transcript.aligns.push_back(alignment);
            transcript.exons.back().length += junctionPos;
            if (del > 0) {
                transcript.exons.push_back(Exon{
                        lastIntronBase + 1 + junctionPos,
                        alignmentReadEnd - lastReadPos - junctionPos,
                        lastReadPos + 1 + junctionPos,
                        alignment.score
                });
                score += int64_t(-0.25 * log2(del));

            } else if (ins > 0) {
                transcript.exons.push_back(Exon{
                        lastGenomePos + junctionPos + 1,
                        alignmentGenomeEnd - alignmentGenomeStart - junctionPos - ins,
                        lastReadPos + junctionPos + ins + 1,
                        alignment.score
                });
            }
            transcript.sj.push_back(
                    SpliceJunction{junctionType, genomeGap - readGap, false});// cannot handle annotated junctions yet

        }

        transcript.score += alignment.score;//add alignment score
        transcript.score += score;// add stitch score
        transcript.matched += transcript.exons.back().length; // update matched bases

        return score;
    }


    //{lengthExtended, mismatch}
    // dir: true for 5' extension, false for 3' extension
    std::pair<int64_t, int64_t> StitchingManagement::extendAlign(WindowAlign &a, bool dir) {
        // extend the alignment in the given direction
        if (dir) {
            //5' extension

        }
    }

    void StitchingManagement::finalizeTranscript(Transcript &t, const int chrId) {
        //try to extend the transcript forward(5') and backward(3')
        std::string chrName = genomeIndex_.genome->chromosomes_[chrId].name;
        int64_t posInChr = t.genomeStart - genomeIndex_.genome->chromosomes_[chrId].start;
        if (t.readStart > 0) {
            // extend 5' end
            int64_t extendLength = std::min(t.readStart, posInChr); // to ensure not to extend beyond current chromosome
            t.readStart -= extendLength;
            t.genomeStart -= extendLength;
            t.exons[0].start -= extendLength;
            t.exons[0].length += extendLength;
            t.exons[0].readStart -= extendLength;
            auto [matched, unmatched] = compareGenomeRead(
                    genomeIndex_.genome->sequence_,
                    read_->sequence[t.strand],
                    t.genomeStart,
                    t.readStart,
                    extendLength
            );
            t.score += matched * MATCH_SCORE + unmatched * MISMATCH_PENALTY;
            t.matched += matched;
            t.unmatched += unmatched;
            posInChr -= extendLength;
        }
        int64_t readEnd = t.exons.back().readStart + t.exons.back().length - 1;
        int64_t genomeEnd = t.exons.back().start + t.exons.back().length - 1;
        if (readEnd < read_->length - 1) {
            // extend 3' end
            int64_t extendLength = read_->length - readEnd - 1;
            t.exons.back().length += extendLength;
            t.aligns.back().length += extendLength;
            auto [matched, unmatched] = compareGenomeRead(
                    genomeIndex_.genome->sequence_,
                    read_->sequence[t.strand],
                    genomeEnd + 1,
                    readEnd + 1,
                    extendLength
            );
            t.score += matched * MATCH_SCORE + unmatched * MISMATCH_PENALTY;
            t.matched += matched;
            t.unmatched += unmatched;
        }

        t.chr = chrName;
        t.posInChr = posInChr;
    }


    TranscriptPtr StitchingManagement::getBestTranscript() const {
        return bestTranscript_;
    }


    std::vector<SJDBOutput> StitchingManagement::getSJDB() const {
        return sjdb;
    }

} // namespace rna

