#include "Stitching.h"
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


    Stitching::Stitching(const StitchingConfig &config,
                         const GenomeIndex &genomeIndex)
            : config_(config), genomeIndex_(genomeIndex) {
        auto winBinMapRaw = new int32_t[(genomeIndex_.suffixArray.fullLength_ >> (config.winBinSizeLog - 1)) +
                                        2];
        for (size_t i = 0; i < (genomeIndex_.suffixArray.fullLength_ >> (config.winBinSizeLog - 1)) + 2; ++i) {
            winBinMapRaw[i] = -1; // initialize all bins to -1
        }
        winBinMap_[0] = winBinMapRaw; //dirction 0
        winBinMap_[1] = winBinMapRaw + (genomeIndex_.suffixArray.fullLength_ >> config.winBinSizeLog) + 1;


        // data structures for stitching
        allWindowAligns_ = new WindowAlign[config_.maxWindows * config_.maxSeedPerWindows];
        nowStitchingRecord_ = new StitchingRecord[config_.maxSeedPerWindows * (config_.maxSeedPerWindows - 1) / 2];
        nowExtensionRecord_[0] = new ExtensionRecord[config_.maxSeedPerWindows];
        nowExtensionRecord_[1] = new ExtensionRecord[config_.maxSeedPerWindows];
        nowRawTranscript_ = new RawTranscript[config_.maxSeedPerWindows *
                                              config_.maxSeedPerWindows]; // actually we only need half of this, but we allocate the full size for simplicity
        goodTranscripts_ = new Transcript[config_.transcriptStoredMax];
        allSingleExtensionRecord_ = new ExtensionRecord::singleExtensionRecord[config_.maxSeedPerWindows * 2 *
                                                                               (1 + config_.maxMismatch)];
        allFragmentStitchingRecord_ = new FragmentStitchingRecord[(config_.maxSeedPerWindows * config_.maxSeedPerWindows) *
                                                                  (config_.maxSeedPerWindows * config_.maxSeedPerWindows)/16 + 1];

        int singleExtensionRecordPtr = 0;
        for (int i = 0; i < config_.maxSeedPerWindows; ++i) {
            nowExtensionRecord_[0][i].maxExtensionLengthWithMismatch =
                    allSingleExtensionRecord_ + singleExtensionRecordPtr;
            singleExtensionRecordPtr += (1 + config_.maxMismatch);
            nowExtensionRecord_[1][i].maxExtensionLengthWithMismatch =
                    allSingleExtensionRecord_ + singleExtensionRecordPtr;
            singleExtensionRecordPtr += (1 + config_.maxMismatch);
        }

        //initialize the outFilter parameters

        outFilterScoreMin_ = 0;
        outFilterMatchMin_ = 0;
    }

    Stitching::~Stitching() {
        delete[] winBinMap_[0];
        delete[] allWindowAligns_;
        delete[] nowStitchingRecord_;
        delete[] nowExtensionRecord_[0];
        delete[] nowExtensionRecord_[1];
        delete[] nowRawTranscript_;
        delete[] goodTranscripts_;
        delete[] allSingleExtensionRecord_;
        delete[] allFragmentStitchingRecord_;
    }

    void Stitching::processAlignments(
            Align *alignments, int alignNum, ReadPtr read) {

        status = SUCCESS;
        if (alignNum == 0) {
            status = FAILED_NO_ALIGNMENTS;
            return;
        }



        read_ = std::move(read);
        readLength = read_->length;
        windows_.reserve(config_.maxWindows);
        transcripts_.reserve(100);
        maxTranscriptScore_ = 0;
        numGoodTranscripts_ = 0;
        outFilterMatchMin_ = int(double(read_->length) * config_.outFilterMatchMinOverLRead);
        outFilterScoreMin_ = int64_t(double(read_->length) * config_.outFilterScoreMinOverLRead);

        identifyAnchors(alignments, alignNum);

        createWindows(alignments, alignNum);

        assignAlignmentsToWindows(alignments, alignNum);

        //generateTranscripts();

        generateTranscriptsNew();

        if (numGoodTranscripts_ == 0) {
            status = FAILED_NO_GOOD_TRANSCRIPT;
            return;
        }

        int trueNumGoodTranscripts = numGoodTranscripts_;

        for (int i = 0; i < numGoodTranscripts_; ++i) {
            if (goodTranscripts_[i].score < maxTranscriptScore_ - config_.multimapScoreRange) {
                trueNumGoodTranscripts = i;
                break;
            }
        }

        if (trueNumGoodTranscripts > config_.outFilterMultimapMax) {
            status = FAILED_TOO_MANY_TRANSCRIPTS;
            return;
        }
        numGoodTranscripts_ = trueNumGoodTranscripts;
    }


    void Stitching::refreshWinBinMap() {
        for (const auto &win: windows_) {
            if (win.startBin > win.endBin) continue;
            int32_t *winBinArray = winBinMap_[win.direction];
            for (int64_t winBin = win.startBin; winBin <= win.endBin; ++winBin) {
                winBinArray[winBin] = -1; // reset bin map
            }

        }
    }

    void Stitching::clear() {
        //clear previous status
        refreshWinBinMap();
        windows_.clear();
        transcripts_.clear();
        totalAnchors_ = 0;
        numGoodTranscripts_ = 0;
    }


    void Stitching::identifyAnchors(Align *alignments, int alignNum) {
        for (int i = 0; i < alignNum; ++i) {
            // Check whether it is an anchor
            Align &alignment = alignments[i];
            if (alignment.rep <= config_.maxAnchorRep) {
                alignment.isAnchor = true;
                totalAnchors_++;
            } else {
                alignment.isAnchor = false;
            }
        }
    }

    inline Stitching::PositiveStrandAlign
    Stitching::convertAlignToPositiveStrand(const Align &align, size_t SAi) {
        int64_t location = genomeIndex_.suffixArray[SAi];
        if (location > genomeIndex_.genomeLength) {
            // reverse strand alignment
            int64_t genomeStart = genomeIndex_.genomeLength * 2 - location - align.length;
            int64_t readStart = read_->length - align.readStart - align.length;
            return {genomeStart, readStart, align.length, 1 - align.direction, align.isAnchor,align.iFragment};
        } else {
            // positive strand alignment
            return {location, align.readStart, align.length, align.direction, align.isAnchor,align.iFragment};
        }
    }

    std::pair<bool, Stitching::sjSplitRecord>
    Stitching::trySplitSJ(const Align &align, int64_t location) {
        if (genomeIndex_.gtf == nullptr) return {false, {}};
        // try to split the alignment
        const int64_t genomeLength = genomeIndex_.genome->originalGenomeLength;
        const int64_t sjSeqLength = genomeIndex_.gtf->sjdbSeqLengthSingleStrand_;
        location -= genomeLength;
        sjSplitRecord ans{};
        if (location < sjSeqLength) {
            // positive strand sj
            ans.direction = 0;
        } else {
            // negative strand sj
            ans.direction = 1;
            location = 2 * sjSeqLength - location - align.length;
        }


        int64_t startInSj = location % genomeIndex_.gtf->config_.sjdbLength;
        if (startInSj < genomeIndex_.gtf->config_.sjdbOverhang &&
            startInSj + align.length > genomeIndex_.gtf->config_.sjdbOverhang) {
            // can be split
            ans.isj = location / genomeIndex_.gtf->config_.sjdbLength;
            ans.donorStart = genomeIndex_.genome->sjDonorStart_[ans.isj] + startInSj;
            ans.donorLength = genomeIndex_.gtf->config_.sjdbOverhang - startInSj;
            ans.acceptorStart = genomeIndex_.genome->sjAcceptorStart_[ans.isj];
            ans.acceptorLength = align.length - ans.donorLength;
            return {true, ans};
        } else {
            return {false, {}};
        }
    }

    void Stitching::createRawWindowsSingle(const PositiveStrandAlign &positiveAlign) {
        //create unflanked window for a single alignment
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
            return;
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

    void Stitching::createWindows(const Align *alignments, int alignNum) {
        for (int alignInd = 0; alignInd < alignNum; ++alignInd) {
            const Align &align = alignments[alignInd];
            if (!align.isAnchor) continue;

            for (size_t i = align.leftSAIndex; i <= align.rightSAIndex; ++i) {
                // handle read crossing sj
                if (genomeIndex_.suffixArray[i] > genomeIndex_.genome->originalGenomeLength) {
                    // alignment on sj record
                    auto [crossSJ, sjRecord] = trySplitSJ(align, genomeIndex_.suffixArray[i]);
                    if (!crossSJ) continue;
                    //create windows for both sides of the sj
                    PositiveStrandAlign donor{
                            sjRecord.donorStart,
                            align.readStart,
                            sjRecord.donorLength,
                            align.direction,
                            align.isAnchor,
                            align.iFragment
                    };
                    PositiveStrandAlign acceptor{
                            sjRecord.acceptorStart,
                            align.readStart + sjRecord.donorLength,
                            sjRecord.acceptorLength,
                            align.direction,
                            align.isAnchor,
                            align.iFragment
                    };
                    createRawWindowsSingle(donor);
                    createRawWindowsSingle(acceptor);
                    continue;
                }
                // Create Windows
                auto positiveAlign = convertAlignToPositiveStrand(align, i);
                createRawWindowsSingle(positiveAlign);
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
            for (int64_t binKey = win.endBin; binKey <= rightKeyBound; binKey++) {
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


    bool Stitching::assignSingleAlignment(Window &win, const PositiveStrandAlign &a, int64_t isj) const {

        // when window is full, check if it needs to be replaced
        if (a.length <= win.minLengthWhenFull && !a.isAnchor) return false; // ignore too short no anchor alignment

        //detect overlap
        WindowAlign *aligns = win.aligns;
        int alignNum = win.alignNum;
        for (int i = 0; i < alignNum; ++i) {
            if (aligns[i].isj == isj && a.genomeStart + aligns[i].readStart == aligns[i].genomeStart + a.readStart \
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
                        for (int j = i + 1; j < alignNum; ++j) {
                            if (aligns[j].readStart > a.readStart) {
                                aligns[j - 1] = WindowAlign{
                                        a.readStart,
                                        a.genomeStart,
                                        a.length,
                                        a.length,
                                        a.isAnchor,
                                        isj,
                                        a.iFragment
                                };
                                return true;
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
                                        a.isAnchor,
                                        isj,
                                        a.iFragment
                                };
                                return true;
                            } else aligns[j + 1] = aligns[j];
                        }

                        // reaching here means that a.readStart is the smallest, insert at the beginning
                        aligns[0] = WindowAlign{
                                a.readStart,
                                a.genomeStart,
                                a.length,
                                a.length,
                                a.isAnchor,
                                isj,
                                a.iFragment
                        };
                        return true;


                    }
                }
                // if not, do nothing
                return false;
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
            if (removePos == -1) return false; // all are anchors, cannot add new alignment


            //remove all shortest non-anchor alignments and insert the new one

            int newAlignNum = 0;
            int posToInsert = -1;
            for (int i = 0; i < alignNum; ++i) {
                if (aligns[i].length > win.minLengthWhenFull || aligns[i].isAnchor) {
                    if (aligns[i].readStart >= a.readStart && posToInsert == -1) {
                        posToInsert = newAlignNum++;
                        aligns[newAlignNum] = aligns[i];
                    } else {
                        aligns[newAlignNum++] = aligns[i];
                    }
                }
            }
            win.alignNum = newAlignNum;
            if (win.minLengthWhenFull >= a.length && !a.isAnchor) return false; // new alignment is too short, ignore
            if (posToInsert == -1) posToInsert = newAlignNum++;
            // reaching here means that a.readStart is the largest



            for (int i = newAlignNum - 1; i > posToInsert; --i) {
                aligns[i] = aligns[i - 1];
            }
            aligns[posToInsert] = WindowAlign{
                    a.readStart,
                    a.genomeStart,
                    a.length,
                    a.length,
                    a.isAnchor,
                    isj,
                    a.iFragment
            };
            win.alignNum = newAlignNum;

            return true;
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
                        a.isAnchor,
                        isj,
                        a.iFragment
                };
                return true;
            } else aligns[i] = aligns[i - 1];
        }
        // reaching here means that a.readStart is the smallest, insert at the beginning
        aligns[0] = WindowAlign{
                a.readStart,
                a.genomeStart,
                a.length,
                a.length,
                a.isAnchor,
                isj,
                a.iFragment
        };
        return true;

    }


    void Stitching::assignAlignmentsToWindows(const Align *alignments, int alignNum) {
        if (genomeIndex_.config.twoDirections) {
            for (int i = 0; i < alignNum; ++i) {
                const auto &align = alignments[i];
                if (align.rep > config_.maxRep)
                    continue; // skip alignments with too many reps

                for (size_t j = align.leftSAIndex; j <= align.rightSAIndex; ++j) {
                    int64_t location = genomeIndex_.suffixArray[j];
                    if (location > genomeIndex_.genome->originalGenomeLength) {
                        auto [crossSJ, sjRecord] = trySplitSJ(align, location);
                        if (!crossSJ) continue;
                        // split the alignment into two parts
                        PositiveStrandAlign donor{
                                sjRecord.donorStart,
                                align.readStart,
                                sjRecord.donorLength,
                                align.direction,
                                align.isAnchor,
                                align.iFragment
                        };
                        PositiveStrandAlign acceptor{
                                sjRecord.acceptorStart,
                                align.readStart + sjRecord.donorLength,
                                sjRecord.acceptorLength,
                                align.direction,
                                align.isAnchor,
                                align.iFragment
                        };
                        size_t binKey = (donor.genomeStart >> config_.winBinSizeLog);
                        auto winId = winBinMap_[donor.direction][binKey];
                        if (winId == -1) continue;
                        Window &winDonor = windows_[winId];
                        assignSingleAlignment(winDonor, donor, sjRecord.isj);
                        binKey = (acceptor.genomeStart >> config_.winBinSizeLog);
                        winId = winBinMap_[acceptor.direction][binKey];
                        if (winId == -1) continue;
                        Window &winAcceptor = windows_[winId];
                        assignSingleAlignment(winAcceptor, acceptor, sjRecord.isj);
                    }

                    auto a = convertAlignToPositiveStrand(align, j);
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


    void Stitching::generateTranscriptsNew() {
        for (auto &window: windows_) {
            stitchWindowsAlignNew(window);
        }

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
                return {6, SCORE_GAP_ATAC};
        } else if (left == 8) { //CT
            if (right == 1) //AC
                return {2, 0};
            else if (right == 11) //GC
                return {4, SCORE_GAP_GCAG};
        } else if (left == 3) { //AT
            if (right == 1) // AC
                return {5, SCORE_GAP_ATAC};
        } else if (left == 11) { //GC
            if (right == 2)
                return {3, SCORE_GAP_GCAG};
        }
        return {0, SCORE_GAP_NON_CANONICAL};
    }


    // stitch a new piece of alignment to the transcript and return the STITCH score. The alignment's match score is not included
    // however, when maintaining the transcript, the alignment's match score is added to the transcript's total score



    void
    Stitching::stitchingBetweenWindowAligns(const rna::WindowAlign &a1, const rna::WindowAlign &a2,
                                            int windowDir, StitchingRecord &record) {
        //calculate the stitching score between two alignments
        //stitch a2 to a1

        const int64_t a2ReadEnd = a2.readStart + a2.length - 1;
        const int64_t a2GenomeEnd = a2.genomeStart + a2.length - 1;
        const int64_t a1ReadEnd = a1.readStart + a1.length - 1;
        const int64_t a1GenomeEnd = a1.genomeStart + a1.length - 1;

        if (a1.iFragment != a2.iFragment){
            // stitching between paired reads
            // handle mates extension while DP
            record.type = StitchingRecord::CROSS_FRAGMENTS;
            record.score = a2.length * MATCH_SCORE;
            record.matches = a2.length;
            record.mismatches = 0;
            record.stitchingType.type = -3;
            record.formerExonLengthShift = 0;
            record.latterExonLengthShift = 0;
            return;
        }



        if(a1.isj != -1 && a1.isj == a2.isj && a2.readStart == a1ReadEnd + 1 && a1GenomeEnd + 1 < a2.genomeStart ){
            // annotated splice junction
            if (genomeIndex_.genome->sjdb[a2.isj].motif == 0
                &&(a2.length <= genomeIndex_.genome->sjdb[a2.isj].shiftRight||
                a1.length <= genomeIndex_.genome->sjdb[a2.isj].shiftLeft)) {
                // too large repeats around non-canonical
                record.type = StitchingRecord::CANNOT_STITCH;
                return;
            }
            record.type = StitchingRecord::SPLICE_JUNCTION;
            record.score = SCORE_ANNOTATED_SJ + a2.length * MATCH_SCORE;
            record.formerExonLengthShift = 0;
            record.latterExonLengthShift = 0;
            record.matches = a2.length;
            record.mismatches = 0;
            record.stitchingType = {genomeIndex_.genome->sjdb[a2.isj].motif, a2.genomeStart - (a1GenomeEnd + 1), true};
            return;
        }


        if (a2ReadEnd <= a1ReadEnd || a2GenomeEnd <= a1GenomeEnd) {
            // r or g fully overlap, cannot stitch
            record.type = StitchingRecord::CANNOT_STITCH;
            return;
        }

        int64_t a2ReadStart = a2.readStart;
        int64_t a2GenomeStart = a2.genomeStart;
        int64_t a2Length = a2.length;
        std::string &genomeSeq = genomeIndex_.genome->sequence_;
        std::string &readSeq = read_->sequence[windowDir];

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
            record.type = StitchingRecord::PERFECT_MATCH;
            record.score = (a2.length - a2Length) * MATCH_SCORE;
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

            record.type = StitchingRecord::SAME_GAP;
            record.score = (nMatch + a2Length) * MATCH_SCORE + nMismatch * MISMATCH_PENALTY;
            record.formerExonLengthShift = -a1.length;
            record.latterExonLengthShift = a1.length + readGap;
            record.matches = nMatch + a2Length;
            record.mismatches = nMismatch;
        } else if (genomeGap > readGap) {
            // means there is a deletion or an intron (long deletion)
            Del = genomeGap - readGap;
            record.score = MATCH_SCORE * a2Length;
            if (Del > MAX_INTRON_LENGTH) {
                record.type = StitchingRecord::CANNOT_STITCH;// too large deletion
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
                    genomeSeq[a1GenomeEnd + junctionReadPosTemp]) {
                    scoreTemp += MATCH_SCORE;
                }
                if (readSeq[a1ReadEnd + junctionReadPosTemp] ==
                    genomeSeq[lastIntronBase + junctionReadPosTemp] &&
                    readSeq[a1ReadEnd + junctionReadPosTemp] !=
                    genomeSeq[a1GenomeEnd + junctionReadPosTemp]) {
                    scoreTemp -= MATCH_SCORE;
                }

                int junctionTypeTemp = SIMPLE_DELETION;
                int64_t junctionPenaltyTemp = 0;
                int64_t maxDelScoreTemp = scoreTemp;

                if (Del >= MIN_INTRON_LENGTH) {
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


            if (junctionType < 0 || (MAX_MISMATCH_FOR_SJ[(junctionType + 1) / 2] == -1 ||
                                     nMismatch <= MAX_MISMATCH_FOR_SJ[(junctionType + 1) / 2])) {
                // valid junction
            } else {
                // not a valid junction
                record.type = StitchingRecord::CANNOT_STITCH;
                return;
            }

            // score the gap
            // todo check if the junction is annotated
            if (genomeIndex_.genome->sjdbNum > 0) {

            }

            if (Del > MIN_INTRON_LENGTH) {
                record.score += GAP_OPEN_PENALTY + junctionPenalty;
                record.type = StitchingRecord::SPLICE_JUNCTION;

            } else {
                record.score += DEL_OPEN_PENALTY + Del * DEL_EXTEND_PENALTY;
                record.type = StitchingRecord::DELETION;
            }
            record.formerExonLengthShift = junctionReadPos;
            record.latterExonLengthShift = (a2ReadEnd - a1ReadEnd - junctionReadPos) - a2.length;
            record.matches = nMatch + a2Length;
            record.mismatches = nMismatch;
            int sjStrand = 0;
            if (junctionType > 0) sjStrand = 2- (junctionType % 2);
            record.stitchingType = {junctionType, Del, false,0,0,sjStrand};

        } else if (readGap > genomeGap) {
            //insert
            record.score = MATCH_SCORE * a2Length;
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

            }

            record.score += Ins * INS_EXTEND_PENALTY + INS_OPEN_PENALTY;
            junctionType = SIMPLE_INSERTION; // mark insertion
            record.type = StitchingRecord::INSERTION;
            record.formerExonLengthShift = junctionReadPos;
            record.latterExonLengthShift = (a2ReadEnd - a1ReadEnd - junctionReadPos - Ins) - a2.length;
            record.matches = nMatch + a2Length;
            record.mismatches = nMismatch;
            record.stitchingType = {junctionType, Ins, false};
        }


    }

    // extend the alignment
    void
    Stitching::extendWindowAlign(const WindowAlign &a, int windowDir, int extendDir, ExtensionRecord &res) {
        std::string &genomeSeq = genomeIndex_.genome->sequence_;
        std::string &readSeq = read_->sequence[windowDir];


        int mismatchCount = 0;
        int matchCount = 0;
        int maxExtendScore = 0;
        int maxScoreMismatch = 0;
        int maxScoreLength = 0;
        int maxScoreMatch = 0;
        int extendScore = 0;
        if (extendDir == 0) {
            // extend forward
            int64_t extendLength = a.readStart;
            int l = 1;
            for (l = 1; l <= extendLength; ++l) {
                int64_t genomePos = a.genomeStart - l;
                int64_t readPos = a.readStart - l;
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
                    if (mismatchCount > config_.maxMismatch) break; // stop extending if too many mismatches
                } else {
                    ++matchCount;
                    ++extendScore;
                }
            }

            if (mismatchCount <= config_.maxMismatch) {
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
                res.maxMismatch = config_.maxMismatch;
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
                    if (mismatchCount > config_.maxMismatch) break; // stop extending if too many mismatches
                } else {
                    ++matchCount;
                    ++extendScore;
                }
            }

            if (mismatchCount <= config_.maxMismatch) {
                if (extendScore > maxExtendScore) {
                    maxExtendScore = extendScore;
                    maxScoreLength = (int) l-1;
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
                res.maxMismatch = config_.maxMismatch;
            }

            res.maxExtensionScore = maxExtendScore;

        }

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


    void Stitching::stitchWindowsAlignSingle(rna::Window &window) {
        int nAligns = window.alignNum;
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
                    t.startAlignId = i;
                    t.crossFragments = false;
                    t.firstFragmentMatchEnd = 0;
                } else {
                    RawTranscript &t = nowRawTranscript_[getRawTranscriptIndex(i, k, nAligns)];
                    t.score = -1000000;
                    t.startAlignId = i;
                    for (int j = 0; j < k; ++j) {
                        //find [i,i+j], try stitching i+j and i+k
                        int formerTranscriptId = getRawTranscriptIndex(i, j, nAligns);
                        RawTranscript &formerT = nowRawTranscript_[formerTranscriptId];
                        if (formerT.score < 0) continue; // invalid transcript
                        StitchingRecord &stitchingR = nowStitchingRecord_[getStitchingRecordIndex(i + j, i + k,
                                                                                                  nAligns)];
                        if (stitchingR.type == StitchingRecord::CANNOT_STITCH) continue; // cannot stitch

                        //handle fragment crossing
                        if (stitchingR.type == StitchingRecord::CROSS_FRAGMENTS) {
                            // handle alignEndsProtrude
                            if (window.aligns[i].readStart + window.aligns[i+k].genomeStart + ALIGN_ENDS_PROTRUDE < window.aligns[i].genomeStart ) continue;

                            //todo need to fix
                            if (window.aligns[i+k].genomeStart + readLength - window.aligns[i+k].readStart +readLength < window.aligns[i+j].genomeStart +window.aligns[i+j].length) continue;

                            // todo check junction consistency
                            // todo add parameter alignMateGapMax
                            const ExtensionRecord &e1 = nowExtensionRecord_[1][i + j]; // extend backward of former alignment
                            int nowMismatchToleranceRemaining = config_.maxMismatch - formerT.mismatches;
                            auto extendedInfo1 = e1.maxExtensionLengthWithMismatch[std::min(e1.maxMismatch, nowMismatchToleranceRemaining)];
                            nowMismatchToleranceRemaining -= extendedInfo1.mismatches;
                            const ExtensionRecord &e2 = nowExtensionRecord_[0][i + k]; // extend forward of latter alignment
                            int64_t maxLatterExtensionLength = window.aligns[i].readStart + window.aligns[i+k].genomeStart - window.aligns[i].genomeStart;

                            int64_t maxLatterFragmentMismatch = std::min(e2.maxMismatch, nowMismatchToleranceRemaining);
                            auto extendedInfo2 = e2.maxExtensionLengthWithMismatch[std::min(e2.maxMismatch, nowMismatchToleranceRemaining)];

                            //todo binary search
                            if (extendedInfo2.length > maxLatterExtensionLength){
                                int nowMismatch = extendedInfo2.mismatches;
                                while (extendedInfo2.length > maxLatterExtensionLength){
                                    if (nowMismatch == 0){
                                        extendedInfo2.length = maxLatterExtensionLength;
                                        extendedInfo2.mismatches = 0;
                                        extendedInfo2.matched = maxLatterExtensionLength;
                                        break;
                                    }else {
                                        nowMismatch--;
                                        extendedInfo2 = e2.maxExtensionLengthWithMismatch[nowMismatch];
                                        if (extendedInfo2.length <= maxLatterExtensionLength){
                                            extendedInfo2.length = maxLatterExtensionLength;
                                            extendedInfo2.mismatches = nowMismatch + 1;
                                            extendedInfo2.matched = maxLatterExtensionLength - extendedInfo2.mismatches;
                                            break;
                                        }
                                    }
                                }
                            }


                            int64_t nowScore = formerT.score + stitchingR.score + extendedInfo1.matched * MATCH_SCORE + extendedInfo1.mismatches * MISMATCH_PENALTY
                                               + extendedInfo2.matched * MATCH_SCORE + extendedInfo2.mismatches * MISMATCH_PENALTY;

                            if (nowScore > t.score || (nowScore == t.score && formerT.exonCount > t.exonCount)){
                                t.score = nowScore;
                                t.mismatches = formerT.mismatches + stitchingR.mismatches + extendedInfo1.mismatches + extendedInfo2.mismatches;
                                t.matches = formerT.matches + stitchingR.matches + extendedInfo1.matched + extendedInfo2.matched;
                                t.previousTranscriptId = formerTranscriptId;
                                t.newAlignId = i + k;
                                t.exonCount = formerT.exonCount + 1;
                                t.crossFragments = true;
                                t.extendedLengthFragmentFormer = extendedInfo1.length;
                                t.extendedLengthFragmentLatter = extendedInfo2.length;
                                t.firstFragmentMatchEnd = window.aligns[i+j].genomeStart + window.aligns[i+j].length;

                            }


                            continue;

                        }

                        if (formerT.mismatches + stitchingR.mismatches > config_.maxMismatch)
                            continue; // too many mismatches
                        if (formerT.exonCount == config_.maxExons) continue; // too many exons

                        int64_t nowScore = formerT.score + stitchingR.score;
                        if (nowScore > t.score || (nowScore == t.score && formerT.exonCount > t.exonCount)) {
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
                            t.crossFragments = formerT.crossFragments;
                            t.extendedLengthFragmentFormer = formerT.extendedLengthFragmentFormer;
                            t.extendedLengthFragmentLatter = formerT.extendedLengthFragmentLatter;
                            t.firstFragmentMatchEnd = formerT.firstFragmentMatchEnd;
                        }
                    }
                }
            }
        }
    }

    void Stitching::stitchBetweenFragments(const rna::RawTranscript &t1, const rna::RawTranscript &t2,
                                           rna::RawTranscript &record) {

    }

    inline int getPairedStitchingRecordIndex(const int s1,const int e1,const int s2,const int e2,const int nFirstFragmentAligns,const int nSecondFragmentAligns) {
        // calculate the index of the stitching record between two alignments a,b (0-based)
        return (s1 * nFirstFragmentAligns + e1) * (nSecondFragmentAligns * (nSecondFragmentAligns -1)/2) + (s2 * nSecondFragmentAligns + e2);
    }

    void Stitching::stitchWindowsAlignPaired(rna::Window &window) {
        int nAligns = window.alignNum;
        int nFirstFragmentAligns = window.firstFragmentAlignNum;
        int nSecondFragmentAligns = window.secondFragmentAlignNum;

        int64_t maxFirstFragmentScore = -1000000;
        int64_t maxSecondFragmentScore = -1000000;

        //Stitch first fragment and second fragment separately
        for (int k = 0; k < nFirstFragmentAligns; ++k){
            for (int i = 0; i < nFirstFragmentAligns - k; ++i){
                if (k == 0){
                    RawTranscript &t = nowRawTranscript_[getRawTranscriptIndex(i, 0, nFirstFragmentAligns)];
                    t.previousTranscriptId = -1;
                    t.newAlignId = i;
                    t.score = window.aligns[i].score;
                    t.mismatches = 0;
                    t.matches = window.aligns[i].length;
                    t.exonCount = 1;
                    t.startAlignId = i;
                    t.crossFragments = false;
                    t.firstFragmentMatchEnd = 0;
                } else {
                    RawTranscript &t = nowRawTranscript_[getRawTranscriptIndex(i, k, nFirstFragmentAligns)];
                    t.score = -1000000;
                    t.startAlignId = i;
                    for (int j = 0; j < k; ++j) {
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
                        if (nowScore > t.score || (nowScore == t.score && formerT.exonCount > t.exonCount)) {
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
                            t.crossFragments = formerT.crossFragments;
                            t.firstFragmentMatchEnd = formerT.firstFragmentMatchEnd;
                            maxFirstFragmentScore = std::max(maxFirstFragmentScore, t.score);
                        }
                    }
                }
            }
        }


        for (int k = 0; k < nSecondFragmentAligns; ++k){
            for (int i = nFirstFragmentAligns; i < nAligns - k; ++i){
                if (k == 0){
                    RawTranscript &t = nowRawTranscript_[getRawTranscriptIndex(i, 0, nFirstFragmentAligns)];
                    t.previousTranscriptId = -1;
                    t.newAlignId = i;
                    t.score = window.aligns[i].score;
                    t.mismatches = 0;
                    t.matches = window.aligns[i].length;
                    t.exonCount = 1;
                    t.startAlignId = i;
                    t.crossFragments = false;
                    t.firstFragmentMatchEnd = 0;
                } else {
                    RawTranscript &t = nowRawTranscript_[getRawTranscriptIndex(i, k, nFirstFragmentAligns)];
                    t.score = -1000000;
                    t.startAlignId = i;
                    for (int j = 0; j < k; ++j) {
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
                        if (nowScore > t.score || (nowScore == t.score && formerT.exonCount > t.exonCount)) {
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
                            t.crossFragments = formerT.crossFragments;
                            t.firstFragmentMatchEnd = formerT.firstFragmentMatchEnd;
                            maxSecondFragmentScore = std::max(maxSecondFragmentScore, t.score);
                        }
                    }
                }
            }
        }

        // Now, stitch first fragment transcripts and second fragment transcripts

        // end genome pos of the first fragment: end of i+k
        // start genome pos of the second fragment: start of j

        if (maxFirstFragmentScore + maxSecondFragmentScore < outFilterScoreMin_) {
            return;
        }

        int64_t maxFragmentStitchingScore = -1000000;

        // todo may be optimized
        for (int i = 0; i < nFirstFragmentAligns; ++i){
            for (int ki = 0; ki < nFirstFragmentAligns - i; ++ki){
                //iterate all first fragment transcripts [i,i+ki]
                RawTranscript &firstFragmentT = nowRawTranscript_[getRawTranscriptIndex(i, ki, nFirstFragmentAligns)];
                if (firstFragmentT.score < 0) continue; // invalid transcript
                for (int j = nFirstFragmentAligns; j < nAligns; ++j){
                    for (int kj = 0; kj < nAligns - j; ++kj){
                        RawTranscript &secondFragmentT = nowRawTranscript_[getRawTranscriptIndex(j, kj, nFirstFragmentAligns)];
                        if (secondFragmentT.score < 0) continue; // invalid transcript
                        if (firstFragmentT.score +secondFragmentT.score < maxFragmentStitchingScore) continue; // cannot reach best score
                        if (firstFragmentT.score + secondFragmentT.score < outFilterScoreMin_) continue; // too low score
                        if (firstFragmentT.mismatches + secondFragmentT.mismatches > config_.maxMismatch) continue; // too many mismatches

                        // check fragment consistency
                        int64_t firstFragmentEnd = window.aligns[i + ki].genomeStart + window.aligns[i + ki].length;
                        int64_t secondFragmentStart = window.aligns[j].genomeStart;
                        int64_t secondFragmentEnd = window.aligns[j + kj].genomeStart + window.aligns[j + kj].length;

                        // todo handle alignEndsProtrude
                        if (secondFragmentEnd - firstFragmentEnd < 0) continue; // invalid

                        // check overlap
                        int64_t fragmentGap = secondFragmentStart - firstFragmentEnd;
                        int idx = getPairedStitchingRecordIndex(i, i + ki, j, j + kj, nFirstFragmentAligns, nSecondFragmentAligns);
                        if (fragmentGap < 0){
                            // todo handle overlap, need to check whether the overlapped region is consistent

                            // check consistency
                            int64_t overlapLength = -fragmentGap;
                            int64_t end1 = firstFragmentEnd;
                            int64_t start2 = secondFragmentStart;

                            bool isConsistent = true;
                            while(true){
                                int64_t wa1GenomeStart = window.aligns[firstFragmentT.newAlignId].genomeStart;
                                int prevTranscriptId1 = firstFragmentT.previousTranscriptId;
                                if (wa1GenomeStart > start2){
                                    if (prevTranscriptId1 == -1){
                                        // alignment result of fragment 2 starts before fragment 1, means that a sj is in 2 but not in 1.
                                        isConsistent = false;
                                        break;
                                    }
                                    //todo go to previous alignment of fragment 1
                                }

                            }
                            // then check whether the overlapped region has the same sj




                        }else{
                            //no overlap, stitch directly

                            allFragmentStitchingRecord_[idx].score = firstFragmentT.score + secondFragmentT.score;
                            allFragmentStitchingRecord_[idx].mismatches = firstFragmentT.mismatches + secondFragmentT.mismatches;
                            allFragmentStitchingRecord_[idx].matches = firstFragmentT.matches + secondFragmentT.matches;
                            allFragmentStitchingRecord_[idx].firstFragmentRecord = &firstFragmentT;
                            allFragmentStitchingRecord_[idx].secondFragmentRecord = &secondFragmentT;
                        }



                    }
                }
            }
        }


    }





    void Stitching::stitchWindowsAlignNew(Window &window) {
        //calculate the stitching score between all alignments in the window
        int nAligns = window.alignNum;
        if (nAligns == 0) return;  // no alignment
        if (window.numAnchors == 0) return; // no anchor, cannot stitch

        //todo filter alignSJOverhangMin

        //todo when stitching between paired reads, would it be better to extend the alignments after stitching frag0 and frag1?

        for (int i = 0; i < nAligns; ++i) {
            extendWindowAlign(window.aligns[i], window.direction, 0, nowExtensionRecord_[0][i]);
            extendWindowAlign(window.aligns[i], window.direction, 1, nowExtensionRecord_[1][i]);
        }

        // judge whether it is hopeful to produce a valid transcript
        //todo there might be a better way to estimate the max possible score
        {
            int64_t nowReadStart = 0;
            int64_t nowScore = 0;
            int64_t maxPossibleScore = 0;
            window.firstFragmentAlignNum = 0;
            window.secondFragmentAlignNum = 0;
            for (int i = 0; i < nAligns; ++i) {
                if (nowReadStart < window.aligns[i].readStart) {
                    maxPossibleScore += nowScore;
                    nowScore = 0;
                    nowReadStart = window.aligns[i].readStart;
                }
                nowScore = std::max(window.aligns[i].length + nowExtensionRecord_[0][i].maxExtensionScore +
                                    nowExtensionRecord_[1][i].maxExtensionScore, nowScore);
                if (window.aligns[i].iFragment == 0) ++window.firstFragmentAlignNum;
                else ++window.secondFragmentAlignNum;
            }
            maxPossibleScore += nowScore;
            if (maxPossibleScore < outFilterScoreMin_) return; // cannot produce a valid transcript, skip this window
        }


        int nowIndex = 0;
        for (int i = 0; i < nAligns; ++i) {
            for (int j = i + 1; j < nAligns; ++j) {
                stitchingBetweenWindowAligns(window.aligns[i], window.aligns[j],
                                             window.direction, nowStitchingRecord_[nowIndex]);
                ++nowIndex;
            }
        }



        // DP, calculate the best stitching beginning with i and ending with j, for all alignments i,j
        // this may require n^3 time complexity, but the Constant is very low. We have stored all important information in nowStitchingRecord_ and nowExtensionRecord_
        // filter by max_exons
        if (config_.isPaired) {
            stitchWindowsAlignPaired(window);
        } else {
            stitchWindowsAlignSingle(window);
        }

        // now, we get all best [i,j] raw transcript, try to extend them and finalize them
        // filter them by score, matches and mismatches
        int firstDirToExtend = window.direction;// 5' first
        int64_t localBestScore = -1000000;


        for (int i = 0; i < nAligns; ++i) {
            for (int j = 0; j < nAligns - i; ++j) {
                RawTranscript &t = nowRawTranscript_[getRawTranscriptIndex(i, j, nAligns)];
                if (t.score < 0) continue;
                int maxMismatchRemaining = config_.maxMismatch - t.mismatches;
                // this cannot happen, but just in case of bugs
                if (maxMismatchRemaining < 0) {
                    t.score = -1000000; // too many mismatches, invalid transcript
                    continue;
                }



                int dir = firstDirToExtend;
                int endAlignId = dir == 0 ? i : i + j;
                const ExtensionRecord &e1 = nowExtensionRecord_[firstDirToExtend][endAlignId];
                int extendLength = 0;


                ExtensionRecord::singleExtensionRecord extendInfo;
                int score;

                extendInfo = e1.maxExtensionLengthWithMismatch[std::min(e1.maxMismatch, maxMismatchRemaining)];
                maxMismatchRemaining -= extendInfo.mismatches;
                score = extendInfo.matched * MATCH_SCORE - extendInfo.mismatches * MATCH_SCORE;
                t.score += score;
                t.mismatches += extendInfo.mismatches;
                t.matches += extendInfo.matched;
                extendLength = extendInfo.length;

                if (dir == 0) {
                    t.extendedLengthForward = extendLength;
                } else {
                    t.extendedLengthBackward = extendLength;
                }

                dir = 1 - dir; // switch direction
                endAlignId = dir == 0 ? i : i + j;
                const ExtensionRecord &e2 = nowExtensionRecord_[dir][endAlignId];

                extendInfo = e2.maxExtensionLengthWithMismatch[std::min(e2.maxMismatch, maxMismatchRemaining)];
                maxMismatchRemaining -= extendInfo.mismatches;
                score = extendInfo.matched * MATCH_SCORE - extendInfo.mismatches * MATCH_SCORE;
                t.score += score;
                t.mismatches += extendInfo.mismatches;
                t.matches += extendInfo.matched;
                extendLength = extendInfo.length;

                if (dir == 0) {
                    t.extendedLengthForward = extendLength;
                } else {
                    t.extendedLengthBackward = extendLength;
                }

                //extension finished

                // finalize the transcript by adding the length penalty
                t.score += int64_t(ceil(log2(double(window.aligns[i + j].genomeStart + window.aligns[i + j].length -
                                                    window.aligns[i].genomeStart)) * -0.25 - 0.5));


                if (t.score > localBestScore)
                    localBestScore = t.score;

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
                // check multiplicity

                bool alreadyExists = false;

                Transcript curT;
                curT.chr = chrName;
                curT.strand = window.direction;
                curT.matched = t.matches;
                curT.unmatched = t.mismatches;
                curT.score = t.score;
                curT.exons.resize(t.exonCount);
                curT.sj.resize(t.exonCount - 1);
                // recover the exons and splice junctions, aligns no longer needed

                // add the last exon
                WindowAlign &lastWA = window.aligns[t.newAlignId];
                curT.exons[t.exonCount - 1] = Exon{
                        lastWA.genomeStart,
                        lastWA.length + t.extendedLengthBackward,
                        lastWA.readStart
                };

                int nowRawTranscriptIndex = t.previousTranscriptId;
                int rightAlignId = t.newAlignId;
                int alignPos = t.exonCount - 1; // for convenience
                while (nowRawTranscriptIndex != -1) {

                    const RawTranscript &nowRawT = nowRawTranscript_[nowRawTranscriptIndex];
                    const WindowAlign &nowWindowAlign = window.aligns[nowRawT.newAlignId];
                    const StitchingRecord &stitchingRecord = nowStitchingRecord_[getStitchingRecordIndex(
                            nowRawT.newAlignId, rightAlignId, nAligns)];
                    if (stitchingRecord.type == StitchingRecord::CROSS_FRAGMENTS){
                        // handle fragment crossing
                        if (t.extendedLengthFragmentLatter > 0) {
                            // extend the latter exon
                            curT.exons[alignPos].length += t.extendedLengthFragmentLatter;
                            curT.exons[alignPos].readStart -= t.extendedLengthFragmentLatter;
                            curT.exons[alignPos].start -= t.extendedLengthFragmentLatter;
                        }
                        // add a new exon for the former fragment
                        --alignPos;
                        curT.exons[alignPos].start = nowWindowAlign.genomeStart;
                        curT.exons[alignPos].length = nowWindowAlign.length + t.extendedLengthFragmentFormer;
                        curT.exons[alignPos].readStart = nowWindowAlign.readStart;
                        curT.sj[alignPos] = stitchingRecord.stitchingType;

                        nowRawTranscriptIndex = nowRawT.previousTranscriptId;
                        rightAlignId = nowRawT.newAlignId;
                        continue;

                    }else if (stitchingRecord.type == StitchingRecord::PERFECT_MATCH ||
                        stitchingRecord.type == StitchingRecord::SAME_GAP) {
                        // just extend the exon
                        curT.exons[alignPos].length += stitchingRecord.latterExonLengthShift;
                        curT.exons[alignPos].readStart -= stitchingRecord.latterExonLengthShift;
                        curT.exons[alignPos].start -= stitchingRecord.latterExonLengthShift;
                    } else {
                        //add a new exon and a new splice junction
                        --alignPos;
                        curT.exons[alignPos].start = nowWindowAlign.genomeStart;
                        curT.exons[alignPos].length = nowWindowAlign.length + stitchingRecord.formerExonLengthShift;
                        curT.exons[alignPos].readStart = nowWindowAlign.readStart;
                        curT.exons[alignPos + 1].start -= stitchingRecord.latterExonLengthShift;
                        curT.exons[alignPos + 1].length += stitchingRecord.latterExonLengthShift;
                        curT.exons[alignPos + 1].readStart -= stitchingRecord.latterExonLengthShift;
                        curT.sj[alignPos] = stitchingRecord.stitchingType;
                    }

                    nowRawTranscriptIndex = nowRawT.previousTranscriptId;
                    rightAlignId = nowRawT.newAlignId;
                }

                curT.exons[0].start -= t.extendedLengthForward;
                curT.exons[0].length += t.extendedLengthForward;
                curT.exons[0].readStart -= t.extendedLengthForward;
                curT.readStart = curT.exons[0].readStart;
                curT.genomeStart = curT.exons[0].start;
                curT.posInChr = curT.genomeStart - chrStart;
                curT.chrStartPos = chrStart;
                curT.readLength = readLength;
                curT.CIGAR = curT.getCIGAR();

                for (int iT = 0; iT < numGoodTranscripts_; ++iT) {
                    if (goodTranscripts_[iT].genomeStart == curT.genomeStart &&
                        goodTranscripts_[iT].exons.back().start + goodTranscripts_[iT].exons.back().length ==
                        curT.exons.back().start + curT.exons.back().length &&
                        curT.CIGAR == goodTranscripts_[iT].CIGAR) {

                        // same transcript already exists
                        alreadyExists = true;
                        break;
                    }
                }
                if (alreadyExists) continue; // skip this transcript, already exists

                int posToInsert = numGoodTranscripts_;
                if (numGoodTranscripts_ < config_.transcriptStoredMax) {
                    ++numGoodTranscripts_;
                    for (posToInsert = numGoodTranscripts_ - 1; posToInsert > 0; --posToInsert) {
                        if (goodTranscripts_[posToInsert - 1].score < t.score ||
                            (goodTranscripts_[posToInsert - 1].score == t.score &&
                             goodTranscripts_[posToInsert - 1].matched < t.matches)) {
                            goodTranscripts_[posToInsert] = std::move(goodTranscripts_[posToInsert - 1]);
                        } else {
                            break;
                        }
                    }
                } else {
                    if (goodTranscripts_[posToInsert - 1].score < t.score ||
                        (goodTranscripts_[posToInsert - 1].score == t.score &&
                         goodTranscripts_[posToInsert - 1].matched < t.matches)) {
                        posToInsert = -1; // full, cannot insert
                    } else {
                        for (posToInsert = numGoodTranscripts_ - 1; posToInsert > 0; --posToInsert) {
                            if (goodTranscripts_[posToInsert - 1].score < t.score ||
                                (goodTranscripts_[posToInsert - 1].score == t.score &&
                                 goodTranscripts_[posToInsert - 1].matched < t.matches)) {
                                goodTranscripts_[posToInsert] = std::move(goodTranscripts_[posToInsert - 1]);
                            } else {
                                break;
                            }
                        }
                    }
                }


                if (posToInsert != -1) {
                    Transcript &nowT = goodTranscripts_[posToInsert];
                    nowT = curT;
                }

            }
        }


    }




} // namespace rna

