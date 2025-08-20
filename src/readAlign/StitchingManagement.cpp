#include "StitchingManagement.h"
#include "../utils/exceptions.h"
#include "../utils/types.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include <chrono>

#define SIMPLE_DELETION -1;
#define SIMPLE_INSERTION -2;


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
        auto winBinMapRaw = new int32_t[(genomeIndex_.suffixArray.fullLength_ >> (config.winBinSizeLog - 1)) + 2];
        for (size_t i = 0; i < (genomeIndex_.suffixArray.fullLength_ >> (config.winBinSizeLog - 1)) + 2; ++i) {
            winBinMapRaw[i] = -1; // initialize all bins to -1
        }
        winBinMap_[0] = winBinMapRaw; //dirction 0
        winBinMap_[1] = winBinMapRaw + (genomeIndex_.suffixArray.fullLength_ >> config.winBinSizeLog) + 1;
        allWindowAligns_ = new WindowAlign[config_.maxWindows * config_.maxSeedPerWindows];
        nowStitchingRecord_ = new StitchingRecord[config_.maxSeedPerWindows * (config_.maxSeedPerWindows - 1) / 2];
        nowExtensionRecord_[0] = new ExtensionRecord[config_.maxSeedPerWindows];
        nowExtensionRecord_[1] = new ExtensionRecord[config_.maxSeedPerWindows];
        nowRawTranscript_ = new RawTranscript[config_.maxSeedPerWindows *
                                              config_.maxSeedPerWindows]; // actually we only need half of this, but we allocate the full size for simplicity
        goodTranscripts_ = new Transcript[config_.transcriptStoredMax];
    }

    StitchingManagement::~StitchingManagement() {
        delete[] winBinMap_[0];
        delete[] allWindowAligns_;
        delete[] nowStitchingRecord_;
        delete[] nowExtensionRecord_[0];
        delete[] nowExtensionRecord_[1];
        delete[] nowRawTranscript_;
        delete[] goodTranscripts_;
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
        maxTranscriptScore_ = 0;
        numGoodTranscripts_ = 0;
        outFilterMatchMin_ = int(double(read_->length) * config_.outFilterMatchMinOverLRead);
        outFilterScoreMin_ = int64_t(double(read_->length) * config_.outFilterScoreMinOverLRead);

        identifyAnchors(alignments);

        createWindows(alignments);

        assignAlignmentsToWindows(alignments);

        generateTranscripts();

    }


    void StitchingManagement::refreshWinBinMap() {
        for (const auto &win: windows_) {
            if (win.startBin > win.endBin) continue;
            int32_t *winBinArray = winBinMap_[win.direction];
            for (int64_t winBin = win.startBin; winBin <= win.endBin; ++winBin) {
                winBinArray[winBin] = -1; // reset bin map
            }

        }
    }

    void StitchingManagement::clear() {
        //clear previous status
        refreshWinBinMap();
        windows_.clear();
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
            return {genomeStart, readStart, align.length, 1 - align.direction, align.isAnchor};
        } else {
            // positive strand alignment
            return {location, align.readStart, align.length, align.direction, align.isAnchor};
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
                int64_t baseBinKey = baseBin;
                int64_t dir = positiveAlign.direction;

                if (winBinMap_[dir][baseBinKey] != -1) {
                    // already exists, skip
                    continue;
                }


                // try merging with left existing windows
                int64_t binKey = baseBinKey;
                int64_t leftBound = std::max(chrStartBin, baseBin - config_.winAnchorDistBins);
                int32_t leftOverlap = -1;
                for (binKey = binKey - 1; binKey >= leftBound; binKey--) {
                    if (winBinMap_[dir][binKey] != -1) {
                        leftOverlap = winBinMap_[dir][binKey];
                        break;
                    }
                }
                int32_t nowWindowIndex;
                if (leftOverlap != -1) {
                    // merge with existing window
                    nowWindowIndex = leftOverlap;
                    for (binKey = binKey + 1; binKey <= baseBinKey; binKey++) {
                        winBinMap_[dir][binKey] = nowWindowIndex;
                    }

                } else {
                    //create a new window
                    nowWindowIndex = windows_.size();
                }
                winBinMap_[dir][baseBinKey] = nowWindowIndex;

                // try merging with right existing windows
                int64_t rightBound = std::min(chrEndBin, baseBin + config_.winAnchorDistBins);
                int32_t rightOverlap = -1;
                for (binKey = baseBinKey + 1; binKey <= rightBound; binKey++) {
                    if (winBinMap_[dir][binKey] != -1) {
                        rightOverlap = winBinMap_[dir][binKey];
                    }
                }
                if (rightOverlap != -1) {
                    binKey = baseBinKey + 1;
                    //extend to reach the overlapping window
                    while (winBinMap_[dir][binKey] != rightOverlap) {
                        winBinMap_[dir][binKey] = nowWindowIndex;
                        binKey++;
                    }
                    // merge with right window
                    while (winBinMap_[dir][binKey] == rightOverlap) {
                        winBinMap_[dir][binKey] = nowWindowIndex;
                        binKey++;
                    }

                    window.endBin = binKey;
                    // kill right window
                    windows_[rightOverlap].startBin = 1;
                    windows_[rightOverlap].endBin = 0;
                } else {
                    window.endBin = baseBin;
                }
                window.startBin = baseBin;
                // update window information
                if (leftOverlap == -1) {
                    //new window
                    windows_.emplace_back(window);
                } else {
                    windows_[leftOverlap].endBin = window.endBin;
                }

            }
        }



        // flank existing windows
        for (int32_t i = 0; i < windows_.size(); ++i) {
            auto &win = windows_[i];
            if (win.startBin > win.endBin) continue;
            // flank the window
            auto &chr = genomeIndex_.genome->chromosomes_[win.chrIndex];
            int64_t chrStartBin = chr.start >> config_.winBinSizeLog;
            int64_t chrEndBin = (chr.start + chr.length - 1) >> config_.winBinSizeLog;
            int64_t leftBin = std::max(chrStartBin, win.startBin - config_.flankSize);
            int64_t rightBin = std::min(chrEndBin, win.endBin + config_.flankSize);
            int64_t leftKeyBound = win.startBin;
            int dir = win.direction;
            for (int64_t binKey = leftBin; binKey < leftKeyBound; binKey++) {
                winBinMap_[dir][binKey] = i;
            }
            int64_t rightKeyBound = rightBin;
            for (int64_t binKey = win.endBin; binKey < rightKeyBound; binKey++) {
                winBinMap_[dir][binKey] = i;
            }
            win.startBin = leftBin;
            win.endBin = rightBin;
            // reserve space for alignments
            win.aligns = allWindowAligns_ + i * config_.maxSeedPerWindows;
            memset(win.aligns, 0, sizeof(WindowAlign) * config_.maxSeedPerWindows);
            win.alignNum = 0;
        }


    }


    void StitchingManagement::assignSingleAlignment(Window &win, const PositiveStrandAlign &a) const {

        // when window is full, check if it needs to be replaced
        if (a.length < win.minLengthWhenFull && !a.isAnchor) return; // ignore too short no anchor alignment

        //detect overlap
        WindowAlign *aligns = win.aligns;
        int alignNum = win.alignNum;
        for (int i = 0; i < alignNum; ++i) {
            if (a.genomeStart + aligns[i].readStart == aligns[i].genomeStart + a.readStart \
 && ((a.readStart >= aligns[i].readStart) && a.readStart < aligns[i].readStart + aligns[i].length)\
 || (a.readStart + a.length > aligns[i].readStart &&
     a.readStart + a.length <= aligns[i].readStart + aligns[i].length)) {
                // overlap

                // same
                if (a.genomeStart == aligns[i].genomeStart && a.length == aligns[i].length) return;


                if (a.length > aligns[i].length) {
                    //delete the old alignment and insert the new one
                    //keep the array sorted by readStart

                    //find the position to insert
                    //only need to search one side
                    if (a.readStart > aligns[i].readStart) {
                        for (int j = i + 1; j < alignNum; ++j) {
                            if (aligns[j].readStart > a.readStart) {
                                aligns[j - 1] = WindowAlign{
                                        a.readStart,
                                        a.genomeStart,
                                        a.length,
                                        a.length,
                                        a.isAnchor
                                };
                                break;
                            } else aligns[j - 1] = aligns[j];
                        }
                    } else {
                        for (int j = i - 1; j >= 0; --j) {
                            if (aligns[j].readStart <= a.readStart) {
                                aligns[j + 1] = WindowAlign{
                                        a.readStart,
                                        a.genomeStart,
                                        a.length,
                                        a.length,
                                        a.isAnchor
                                };
                                break;
                            } else aligns[j + 1] = aligns[j];
                        }
                        // reaching here means that a.readStart is the smallest, insert at the beginning
                        aligns[0] = WindowAlign{
                                a.readStart,
                                a.genomeStart,
                                a.length,
                                a.length,
                                a.isAnchor
                        };
                    }
                }
                // if not, do nothing
                return;
            }
        }

        if (a.isAnchor) ++win.numAnchors; // anchor must be added to the window

        // handle the case that there are too many seeds in the window
        if (alignNum == config_.maxSeedPerWindows) {
            // calculate minLengthWhenFull and the alignment to remove
            win.minLengthWhenFull = std::numeric_limits<int>::max();
            int removePos = -1;
            for (int i = 0; i < alignNum; ++i) {
                if (!aligns[i].isAnchor && aligns[i].length < win.minLengthWhenFull) {
                    win.minLengthWhenFull = aligns[i].length;
                    removePos = i;
                }
            }
            if (removePos == -1) return; // all are anchors, cannot add new alignment
            if (win.minLengthWhenFull >= a.length && !a.isAnchor) return; // new alignment is too short, ignore

            //remove the removePos alignment and add the new one

            for (int i = 0; i <= removePos; ++i) {
                if (aligns[i].readStart > a.readStart) {
                    // add here
                    // move [i,removePos-1] to [i+1,removePos]
                    for (int j = i + 1; j <= removePos; ++j) aligns[j] = aligns[j - 1];
                    aligns[i] = WindowAlign{
                            a.readStart,
                            a.genomeStart,
                            a.length,
                            a.length,
                            a.isAnchor
                    };
                    return;
                }
            }
            // reaching here means that a.readStart > alignToBeRemoved.readStart
            for (int i = removePos; i < alignNum - 1; ++i) {
                if (aligns[i + 1].readStart > a.readStart) {
                    // add to i
                    aligns[i] = WindowAlign{
                            a.readStart,
                            a.genomeStart,
                            a.length,
                            a.length,
                            a.isAnchor
                    };
                    return;
                } else aligns[i] = aligns[i + 1];
            }
            // reaching here means that a.readStart is the largest
            aligns[alignNum - 1] = WindowAlign{
                    a.readStart,
                    a.genomeStart,
                    a.length,
                    a.length,
                    a.isAnchor
            };
            return;
        }


        // simple case
        // find the position to insert
        ++win.alignNum;
        for (int i = alignNum; i > 0; --i) {
            if (aligns[i - 1].readStart <= a.readStart) {
                aligns[i] = WindowAlign{
                        a.readStart,
                        a.genomeStart,
                        a.length,
                        a.length,
                        a.isAnchor
                };
                return;
            } else aligns[i] = aligns[i - 1];
        }
        // reaching here means that a.readStart is the smallest, insert at the beginning
        aligns[0] = WindowAlign{
                a.readStart,
                a.genomeStart,
                a.length,
                a.length,
                a.isAnchor
        };

    }

    void StitchingManagement::assignAlignmentsToWindows(const std::vector<Align> &alignments) {
        if (genomeIndex_.config.twoDirections) {
            for (const auto &align: alignments) {
                if (align.rep > config_.maxRep) continue; // skip alignments with too many reps

                for (size_t i = align.leftSAIndex; i <= align.rightSAIndex; ++i) {

                    auto a = convertAlignToPositiveStrand(align, i);
                    size_t binKey = (a.genomeStart >> config_.winBinSizeLog);
                    auto winId = winBinMap_[a.direction][binKey];
                    if (winId == -1) continue;
                    // find the corresponding window
                    Window &win = windows_[winId];
                    assignSingleAlignment(win, a);
                }
            }
        } else {
            std::cout << "Error: not implemented" << std::endl;
            exit(1);
        }

        // no need to sort alignments in each window, as they are inserted in order

    }


    // Generate transcripts for each window and find the best transcript
    void StitchingManagement::generateTranscripts() {

        for (auto &window: windows_) {

            stitchWindowAligns(window);

        }

        for (const auto &t: transcripts_) {
            bestTranscript_ = bestTranscript_->score < t.score ? std::make_shared<Transcript>(t) : bestTranscript_;
        }
    }


    void StitchingManagement::stitchWindowAligns(Window &window) {
        //simple DP
        size_t nAligns = window.alignNum;
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
                int jType = -1;
                int64_t junctionPenalty = 0;
                if (del >= MIN_INTRON_LENGTH) {
                    auto res = checkJunctionMotif(
                            genomeIndex_.genome->sequence_,
                            lastGenomePos + junctionPos,
                            lastIntronBase + junctionPos
                    );
                    jType = res.first;
                    junctionPenalty = res.second;
                }


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


    StitchingRecord
    StitchingManagement::stitchingBetweenWindowAligns(const rna::WindowAlign &a1, const rna::WindowAlign &a2,
                                                      int windowDir) {
        //calculate the stitching score between two alignments
        //stitch a2 to a1
        StitchingRecord record;
        const char *genomeSeq = genomeIndex_.genome->sequence_.data();
        const char *readSeq = read_->sequence[windowDir].data();

        int64_t a2ReadEnd = a2.readStart + a2.length - 1;
        int64_t a2GenomeEnd = a2.genomeStart + a2.length - 1;
        int64_t a1ReadEnd = a1.readStart + a1.length - 1;
        int64_t a1GenomeEnd = a1.genomeStart + a1.length - 1;
        int64_t a2ReadStart = a2.readStart;
        int64_t a2GenomeStart = a2.genomeStart;
        int64_t a2Length = a2.length;

        if (a2ReadEnd < a1ReadEnd || a2GenomeEnd < a1GenomeEnd) {
            // r or g fully overlap, cannot stitch
            record.type = StitchingRecord::CANNOT_STITCH;
            return record;
        }

        if (a2.readStart < a1ReadEnd) {
            // a2 5' overlaps a1 3'
            a2GenomeStart += (a1ReadEnd - a2.readStart + 1);
            a2ReadStart = a1ReadEnd + 1;
            a2Length = a2ReadEnd - a2ReadStart + 1;
        }

        int64_t readGap = a2ReadStart - (a1ReadEnd + 1);
        int64_t genomeGap = a2GenomeStart - (a1GenomeEnd + 1);

        int64_t nMatch = 0, nMismatch = 0;
        bool delFlag = false, insFlag = false;
        int64_t Del = 0, Ins = 0; // deletion and insertion length
        int64_t junctionReadPos = 0;
        int64_t junctionType = 0;
        int64_t lastIntronBase = a2GenomeStart - readGap - 1; // if all gap belongs to acceptor

        if (genomeGap == 0 && readGap == 0) {
            // perfect match, no gap
            record.type = StitchingRecord::PERFECT_MATCH;
            record.score = (a2.length - a2Length) * MATCH_SCORE;
            record.formerExonLengthShift = a2Length;
            record.latterExonLengthShift = -a2.length;
            record.matches = a2Length - a2.length;
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

            record.type = StitchingRecord::SAME_GAP;
            record.score = (nMatch + a2Length) * MATCH_SCORE + nMismatch * MISMATCH_PENALTY;
            record.formerExonLengthShift = 0; // if readGap > 0, there's no overlap, so no length shift
            record.latterExonLengthShift = 0;
            record.matches = nMatch + a2Length;
            record.mismatches = nMismatch;
        } else if (genomeGap > readGap) {
            // means there is a deletion or an intron (long deletion)
            delFlag = true;
            Del = genomeGap - readGap;
            record.score = MATCH_SCORE * a2Length;
            if (Del > MAX_INTRON_LENGTH) {
                record.type = StitchingRecord::CANNOT_STITCH;// too large deletion
            }

            // Now, find the best junction pos
            int64_t junctionReadPosTemp = 1;
            int64_t scoreTemp = 0;
            do {
                junctionReadPosTemp--;
                if (readSeq[a1ReadEnd + junctionReadPosTemp] !=
                    genomeSeq[lastIntronBase + junctionReadPosTemp] &&
                    readSeq[a1ReadEnd + junctionReadPosTemp] ==
                    genomeSeq[a2GenomeEnd + junctionReadPosTemp]) {
                    scoreTemp -= MATCH_SCORE;
                }
            } while (scoreTemp + SCORE_STITCH_SJ_SHIFT >= 0 && a1.length + junctionReadPosTemp > 1);

            int64_t maxDelScore = -1000000;
            int64_t junctionPenalty = 0;
            scoreTemp = 0;

            do {
                if (readSeq[a1ReadEnd + junctionReadPosTemp] !=
                    genomeSeq[lastIntronBase + junctionReadPosTemp] &&
                    readSeq[a1ReadEnd + junctionReadPosTemp] ==
                    genomeSeq[a2GenomeEnd + junctionReadPosTemp]) {
                    scoreTemp += MATCH_SCORE;
                }
                if (readSeq[a1ReadEnd + junctionReadPosTemp] ==
                    genomeSeq[lastIntronBase + junctionReadPosTemp] &&
                    readSeq[a1ReadEnd + junctionReadPosTemp] !=
                    genomeSeq[a2GenomeEnd + junctionReadPosTemp]) {
                    scoreTemp -= MATCH_SCORE;
                }

                int junctionTypeTemp = SIMPLE_DELETION;
                int64_t junctionPenaltyTemp = 0;
                int64_t maxDelScoreTemp = scoreTemp;

                if (Del >= MIN_INTRON_LENGTH) {
                    auto res = checkJunctionMotif(
                            genomeSeq,
                            a2GenomeEnd + junctionReadPosTemp,
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

            // todo maintain repeat length for junctions


            // score donor and acceptor
            for (int64_t i = std::min((int64_t) 1, junctionReadPos + 1); i <= std::max(readGap, junctionReadPos); ++i) {
                size_t pos = (i <= junctionReadPos) ? (a1GenomeEnd + i) : (lastIntronBase + i);
                if (genomeSeq[pos] != 'N' && readSeq[a1ReadEnd + i] != 'N') {
                    if (genomeSeq[pos] != readSeq[a1ReadEnd + i]) {
                        nMismatch++;
                        record.score -= MATCH_SCORE;
                        if (i < 1 || i > readGap) {
                            record.score -= MATCH_SCORE;
                            nMatch--;
                        }
                    } else {
                        if (i >= 1 && i <= readGap) {
                            nMatch++;
                            record.score += MATCH_SCORE;
                        }
                    }
                }
            }

            // score the gap
            // todo check if the junction is annotated

            if (Del > MIN_INTRON_LENGTH) {
                record.score += GAP_OPEN_PENALTY + junctionPenalty;
                record.type = StitchingRecord::SPLICE_JUNCTION;

            } else {
                junctionType = SIMPLE_DELETION; // mark deletion
                record.score += DEL_OPEN_PENALTY + Del * DEL_EXTEND_PENALTY;
                record.type = StitchingRecord::DELETION;
            }
            record.formerExonLengthShift = junctionReadPos;
            record.latterExonLengthShift = (a2ReadEnd - a1ReadEnd - junctionReadPos) - a2Length;
            record.matches = nMatch + a2Length;
            record.mismatches = nMismatch;
            record.stitchingType = {junctionType, Del, false};

        } else if (readGap > genomeGap) {
            record.score = MATCH_SCORE * a2Length;
            insFlag = true;
            Ins = readGap - genomeGap;
            if (genomeGap == 0) {
                junctionReadPos = 0;
            } else if (genomeGap < 0) {
                // overlapping
                junctionReadPos = 0;
                record.score -= (-genomeGap) * MATCH_SCORE; // reduce score
            } else {

                int64_t maxInsScore = 0;
                int64_t scoreTemp = 0;
                int64_t junctionReadPosTemp = 1;
                // find the best junction pos
                for (junctionReadPosTemp = 1; junctionReadPosTemp <= genomeGap; junctionReadPosTemp++) {
                    if (genomeSeq[a1GenomeEnd + junctionReadPosTemp] != 'N') {
                        scoreTemp += (readSeq[a1ReadEnd + junctionReadPosTemp] ==
                                      genomeSeq[a1GenomeEnd + junctionReadPosTemp]) ? MATCH_SCORE : -MATCH_SCORE;
                        scoreTemp += (readSeq[a1ReadEnd + Ins + junctionReadPosTemp] ==
                                      genomeSeq[a1GenomeEnd + junctionReadPosTemp]) ? -MATCH_SCORE : MATCH_SCORE;
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
                            record.score -= MATCH_SCORE;
                        } else {
                            nMatch++;
                            record.score += MATCH_SCORE;
                        }
                    }
                }

                // todo implement the alignInsertionFlush parameter

                record.score += Ins * INS_EXTEND_PENALTY + INS_OPEN_PENALTY;
                junctionType = SIMPLE_INSERTION; // mark insertion
                record.type = StitchingRecord::INSERTION;
                record.formerExonLengthShift = junctionReadPos;
                record.latterExonLengthShift = (a2ReadEnd - a1ReadEnd - junctionReadPos - Ins) - a2Length;
                record.matches = nMatch + a2Length;
                record.mismatches = nMismatch;
            }
        }


        return record;
    }

    // extend the alignment
    ExtensionRecord StitchingManagement::extendWindowAlign(const WindowAlign &a, int windowDir, int extendDir) {
        const char *genomeSeq = genomeIndex_.genome->sequence_.data();
        const char *readSeq = read_->sequence[windowDir].data();
        ExtensionRecord res = ExtensionRecord();
        res.maxExtensionLengthWithMismatch = new ExtensionRecord::singleExtensionRecord[config_.maxMismatch + 1];
        int mismatchCount = 0;
        int matchCount = 0;
        if (extendDir == 0) {
            // extend forward
            int64_t extendLength = a.readStart;
            for (int i = 1; i <= extendLength; ++i) {
                int64_t genomePos = a.genomeStart - i;
                int64_t readPos = a.readStart - i;
                if (genomeSeq[genomePos] == 'N' || readSeq[readPos] == 'N') continue;
                if (genomeSeq[genomePos] == '#') break;
                if (genomeSeq[genomePos] != readSeq[readPos]) {
                    res.maxExtensionLengthWithMismatch[mismatchCount] = {
                            .length = i - 1,
                            .matched = matchCount,
                    };
                    ++mismatchCount;
                    if (mismatchCount > config_.maxMismatch) break; // stop extending if too many mismatches
                } else {
                    ++matchCount;
                }
            }

            if (mismatchCount <= config_.maxMismatch) {
                res.maxMismatch = mismatchCount;
                res.maxExtensionLengthWithMismatch[mismatchCount] = {
                        .length = (int) extendLength,
                        .matched = matchCount
                };
            } else {
                res.maxMismatch = config_.maxMismatch;
            }

        } else {
            // extend backward
            int64_t extendLength = read_->length - a.readStart - a.length;
            int64_t pieceEndReadPos = a.readStart + a.length - 1;
            int64_t pieceEndGenomePos = a.genomeStart + a.length - 1;
            for (int i = 1; i <= extendLength; ++i) {
                int64_t readPos = pieceEndReadPos + i;
                int64_t genomePos = pieceEndGenomePos + i;
                if (genomeSeq[genomePos] == 'N' || readSeq[readPos] == 'N') continue;
                if (genomeSeq[genomePos] == '#') break;
                if (genomeSeq[genomePos] != readSeq[readPos]) {
                    res.maxExtensionLengthWithMismatch[mismatchCount] = {
                            .length = i - 1,
                            .matched = matchCount,
                    };
                    ++mismatchCount;
                    if (mismatchCount > config_.maxMismatch) break; // stop extending if too many mismatches
                } else {
                    ++matchCount;
                }
            }

            if (mismatchCount <= config_.maxMismatch) {
                res.maxMismatch = mismatchCount;
                res.maxExtensionLengthWithMismatch[mismatchCount] = {
                        .length = (int) extendLength,
                        .matched = matchCount
                };
            } else {
                res.maxMismatch = config_.maxMismatch;
            }

        }


        return res;
    }


    inline int getStitchingRecordIndex(int a, int b, int num) {
        // calculate the index of the stitching record between two alignments a,b (0-based)
        if (a >= b) return -1; // invalid
        int t = (num - a);
        return (num * (num - 1) / 2) - (t * (t - 1) / 2) + (b - a - 1);
    }

    inline int getRawTranscriptIndex(int i, int j, int num) {
        return i * num + j;
    }

    void StitchingManagement::stitchWindowsAlignNew(Window &window) {
        //calculate the stitching score between all alignments in the window
        int nAligns = window.alignNum;
        if (nAligns == 0) return;  // no alignment
        if (window.numAnchors == 0) return; // no anchor, cannot stitch

        int nowIndex = 0;
        for (int i = 0; i < nAligns; ++i) {
            for (int j = i + 1; j < nAligns; ++j) {
                nowStitchingRecord_[nowIndex] = stitchingBetweenWindowAligns(window.aligns[i], window.aligns[j],
                                                                             window.direction);
                ++nowIndex;
            }
        }

        for (int i = 0; i < nAligns; ++i) {
            nowExtensionRecord_[0][i] = extendWindowAlign(window.aligns[i], window.direction, 0);
            nowExtensionRecord_[1][i] = extendWindowAlign(window.aligns[i], window.direction, 1);
        }

        // DP, calculate the best stitching beginning with i and ending with j, for all alignments i,j
        // this may require n^3 time complexity, but the Constant is very low. We have stored all important information in nowStitchingRecord_ and nowExtensionRecord_
        // filter by max_exons
        for (int k = 0; k < nAligns; ++k) {
            // calculate [i,i+k] for all i
            for (int i = 0; i < nAligns - k; ++i) {
                if (k == 0) {
                    RawTranscript &t = nowRawTranscript_[getRawTranscriptIndex(i, 0, nAligns)];
                    t.previousTranscriptId = -1;
                    t.newAlignId = i;
                    t.score = window.aligns[i].score;
                    t.mismatches = 0;
                    t.matches = window.aligns[i].length;
                    t.exonCount = 1;
                } else {
                    RawTranscript &t = nowRawTranscript_[getRawTranscriptIndex(i, k, nAligns)];
                    t.score = -1000000;
                    for (int j = 0; j < k; ++j) {
                        //find [i,i+j], try stitching i+j and i+k
                        int formerTranscriptId = getRawTranscriptIndex(i, j, nAligns);
                        RawTranscript &formerT = nowRawTranscript_[formerTranscriptId];
                        if (formerT.score < 0) continue; // invalid transcript
                        StitchingRecord &stitchingR = nowStitchingRecord_[getStitchingRecordIndex(i + j, i + k,
                                                                                                  nAligns)];
                        if (stitchingR.type == StitchingRecord::CANNOT_STITCH) continue; // cannot stitch
                        if (formerT.mismatches + stitchingR.mismatches > config_.maxMismatch)
                            continue; // too many mismatches
                        if (formerT.exonCount == config_.maxExons) continue; // too many exons

                        int64_t nowScore = formerT.score + stitchingR.score;
                        if (nowScore > t.score || (nowScore == t.score && formerT.exonCount < t.exonCount)) {
                            t.score = formerT.score + stitchingR.score;
                            t.mismatches = formerT.mismatches + stitchingR.mismatches;
                            t.matches = formerT.matches + stitchingR.matches;
                            t.previousTranscriptId = formerTranscriptId;
                            t.newAlignId = i + k; // the new alignment is i+k
                            if (stitchingR.type != StitchingRecord::PERFECT_MATCH &&
                                stitchingR.type != StitchingRecord::SAME_GAP) {
                                t.exonCount = formerT.exonCount + 1; // increase exon count
                            } else {
                                t.exonCount = formerT.exonCount; // just extend the last exon
                            }
                        }
                    }
                }
            }
        }

        // now, we get all best [i,j] raw transcript, try to extend them and finalize them
        // todo filter them by score, matches and mismatches
        int firstDirToExtend = window.direction;// 5' first
        int64_t localBestScore = -1000000;


        for (int i = 0; i < nAligns; ++i) {
            for (int j = 0; j < nAligns - i; ++j) {
                RawTranscript &t = nowRawTranscript_[getRawTranscriptIndex(i, j, nAligns)];
                int maxMismatchRemaining = config_.maxMismatch - t.mismatches;
                // this cannot happen, but just in case of bugs
                if (maxMismatchRemaining < 0) {
                    t.score = -1000000; // too many mismatches, invalid transcript
                    continue;
                }
                int dir = firstDirToExtend;
                int endAlignId = dir == 0 ? i : i + j;
                ExtensionRecord &e = nowExtensionRecord_[firstDirToExtend][endAlignId];
                int extendLength = 0;
                if (e.maxMismatch <= maxMismatchRemaining) {
                    maxMismatchRemaining -= e.maxMismatch;
                    int score = e.maxExtensionLengthWithMismatch[e.maxMismatch].matched * MATCH_SCORE -
                                e.maxMismatch * MATCH_SCORE;
                    t.score += score;
                    t.mismatches += e.maxMismatch;
                    extendLength = e.maxExtensionLengthWithMismatch[e.maxMismatch].length;
                } else {
                    int score = e.maxExtensionLengthWithMismatch[maxMismatchRemaining].matched * MATCH_SCORE -
                                maxMismatchRemaining * MATCH_SCORE;
                    t.score += score;
                    t.mismatches += maxMismatchRemaining;
                    extendLength = e.maxExtensionLengthWithMismatch[maxMismatchRemaining].length;
                }

                if (dir == 0) {
                    t.extendedLengthForward = extendLength;
                } else {
                    t.extendedLengthBackward = extendLength;
                }

                dir = 1 - dir; // switch direction
                endAlignId = dir == 0 ? i : i + j;
                e = nowExtensionRecord_[dir][endAlignId];
                if (e.maxMismatch <= maxMismatchRemaining) {
                    maxMismatchRemaining -= e.maxMismatch;
                    int score = e.maxExtensionLengthWithMismatch[e.maxMismatch].matched * MATCH_SCORE -
                                e.maxMismatch * MATCH_SCORE;
                    t.score += score;
                    t.mismatches += e.maxMismatch;
                    extendLength = e.maxExtensionLengthWithMismatch[e.maxMismatch].length;
                } else {
                    int score = e.maxExtensionLengthWithMismatch[maxMismatchRemaining].matched * MATCH_SCORE -
                                maxMismatchRemaining * MATCH_SCORE;
                    t.score += score;
                    t.mismatches += maxMismatchRemaining;
                    extendLength = e.maxExtensionLengthWithMismatch[maxMismatchRemaining].length;
                }

                if (dir == 0) {
                    t.extendedLengthForward = extendLength;
                } else {
                    t.extendedLengthBackward = extendLength;
                }

                //extension finished

                // finalize the transcript by adding the length penalty
                t.score -= int64_t(log2(window.aligns[i + j].genomeStart + window.aligns[i + j].length -
                                        window.aligns[i].genomeStart) / 4);


                if (t.score > localBestScore) localBestScore = t.score;

            }
        }

        if (localBestScore > maxTranscriptScore_) maxTranscriptScore_ = localBestScore;
        else if (localBestScore < std::max(maxTranscriptScore_ - config_.multimapScoreRange, (int64_t) 0))
            return; // no good transcript
        int64_t scoreThreshold = std::max(maxTranscriptScore_ - config_.multimapScoreRange, (int64_t) 0);

        // add the good transcript to the result transcript set

        std::string chrName = genomeIndex_.genome->chromosomes_[window.chrIndex].name;
        int64_t chrStart = genomeIndex_.genome->chromosomes_[window.chrIndex].start;

        for (int i = 0; i < nAligns; ++i) {
            for (int j = 0; j < nAligns - i; ++j) {
                RawTranscript &t = nowRawTranscript_[getRawTranscriptIndex(i, j, nAligns)];
                if (t.score < scoreThreshold) continue; // not a good transcript
                if (t.score < outFilterScoreMin_) continue;
                if (t.matches < outFilterMatchMin_) continue;

                // convert the RawTranscript to a Transcript object
                // find the pos to insert
                // todo will it be better to store in ascending order?

                ++numGoodTranscripts_;
                for (int iT = numGoodTranscripts_ - 1; iT >= 0; --iT) {
                    Transcript &curT = goodTranscripts_[iT];
                    if (curT.score < t.score || (curT.score == t.score && curT.matched < t.matches)) {
                        if (iT + 1 < config_.transcriptStoredMax){
                            goodTranscripts_[iT + 1] = std::move(curT);
                        }
                    }else {
                        // found the pos to insert
                        curT.chr = chrName;
                        curT.strand = window.direction;
                        curT.matched = t.matches;
                        curT.unmatched = t.mismatches;
                        curT.score = t.score;
                        curT.readStart = window.aligns[t.newAlignId].readStart - t.extendedLengthForward;
                        curT.genomeStart = window.aligns[t.newAlignId].genomeStart - t.extendedLengthForward;
                        curT.posInChr = curT.genomeStart - chrStart;
                        curT.exons.resize(t.exonCount);
                        curT.sj.resize(t.exonCount - 1);
                        // recover the exons and splice junctions, aligns no longer needed

                        // add the last exon
                        WindowAlign& lastWA = window.aligns[t.newAlignId];
                        curT.exons[t.exonCount - 1] = Exon{
                            lastWA.genomeStart,
                            lastWA.length+ t.extendedLengthBackward,
                            lastWA.readStart
                        };

                        int nowRawTranscriptIndex = t.previousTranscriptId;
                        int rightAlignId = t.newAlignId;
                        int alignPos = t.exonCount - 1; // for convenience
                        while (nowRawTranscriptIndex != -1){

                            RawTranscript& nowRawT = nowRawTranscript_[nowRawTranscriptIndex];
                            WindowAlign& nowWindowAlign = window.aligns[nowRawT.newAlignId];
                            StitchingRecord& stitchingRecord = nowStitchingRecord_[getStitchingRecordIndex(
                                    nowRawT.newAlignId, rightAlignId, nAligns)];
                            if (stitchingRecord.type == StitchingRecord::PERFECT_MATCH ||
                                stitchingRecord.type == StitchingRecord::SAME_GAP) {
                                // just extend the exon
                                curT.exons[alignPos].length += nowRawT.extendedLengthBackward;
                            } else {
                                //add a new exon and a new splice junction
                                --alignPos;
                                curT.exons
                            }


                        }





                    }
                }


            }
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
        t.score += int64_t(-0.25 * log2(t.exons.back().start + t.exons.back().length - t.exons[0].start));

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

