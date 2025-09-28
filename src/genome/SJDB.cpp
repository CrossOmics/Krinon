#include "SJDB.h"
#include "../utils/types.h"
#include "GenomeIndex.h"
#include "Genome.h"
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include "omp.h"

#define GTAG 1
#define CTAC 2
#define GCAG 3
#define CTGC 4
#define ATAC 5
#define GTAT 6
#define NON_CANONICAL 0

#define MARK_END -999
const int SJDB_PADDING_LENGTH = 20; // todo move into define

namespace rna {
    // load exons from GTF file, store in exonLoci_
    void GTF::loadGTF(const std::string &gtfFile, Genome &genome) {
        std::ifstream sjdbStreamIn(gtfFile);
        char sjdbReadBuffer[65536];
        sjdbStreamIn.rdbuf()->pubsetbuf(sjdbReadBuffer, sizeof sjdbReadBuffer);
        std::map<std::string, int64_t> transcriptIDNumber, geneIDNumber;
        if (sjdbStreamIn.fail()) {
            std::cerr << "Error: cannot open GTF file " << gtfFile << "\n";
            exit(1);
        }

        if (genome.chromosomeNameToIndex_.empty()) {
            for (int i = 0; i < genome.chromosomes_.size(); ++i) {
                genome.chromosomeNameToIndex_[genome.chromosomes_[i].name] = i;
            }
        }

        int exonN = 0;
        while (sjdbStreamIn.good()) {//count the number of exons

            std::string oneLine, chr1, ddd2, featureType;
            getline(sjdbStreamIn, oneLine);
            std::istringstream oneLineStream(oneLine);

            oneLineStream >> chr1 >> ddd2 >> featureType;
            if (chr1.substr(0, 1) != "#" && featureType == config_.sjdbGTFfeatureExon) {
                exonN++;
            };
        };
        if (exonN == 0) return;

        sjdbStreamIn.clear();
        sjdbStreamIn.seekg(0, std::ios::beg);

        exonLoci_.resize(exonN);
        exonN = 0;
        while (sjdbStreamIn.good()) {//read the exons
            std::string oneLine, chr1, ddd2, featureType;
            getline(sjdbStreamIn, oneLine);
            std::istringstream oneLineStream(oneLine);

            oneLineStream >> chr1 >> ddd2 >> featureType;
            if (chr1.substr(0, 1) != "#" && featureType == config_.sjdbGTFfeatureExon) {
                //exon line
                if (!config_.sjdbGTFChrPrefix.empty()) chr1 = config_.sjdbGTFChrPrefix + chr1;

                if (genome.chromosomeNameToIndex_.find(chr1) == genome.chromosomeNameToIndex_.end()) {
                    std::cerr << "Warning: chromosome " << chr1 << " not found in genome, skipping this exon\n";
                    continue;
                }

                size_t exonStart, exonEnd;
                char str1;
                oneLineStream >> exonStart >> exonEnd >> ddd2 >> str1 >> ddd2;
                if (exonEnd > genome.chromosomes_[genome.chromosomeNameToIndex_[chr1]].length) {
                    std::cerr << "Warning: exon end " << exonEnd << " exceeds chromosome " << chr1
                              << " length, skipping this exon\n";
                    continue;
                }
                std::string nowLine;
                getline(oneLineStream, nowLine);
                std::replace(nowLine.begin(), nowLine.end(), ';', ' ');
                std::replace(nowLine.begin(), nowLine.end(), '=', ' ');
                std::replace(nowLine.begin(), nowLine.end(), '\t', ' ');
                std::replace(nowLine.begin(), nowLine.end(), '\"', ' ');

                //trID, gID, gName, gBiotype
                std::vector<std::vector<std::string>> exAttrNames(
                        {{config_.sjdbGTFTagExonParentTranscriptId}, {config_.sjdbGTFTagExonParentGene},
                         config_.sjdbGTFTagExonParentGeneName, config_.sjdbGTFTagExonParentGeneType});

                std::vector<std::string> exAttr; //trID, gID, gName, gBiotype
                exAttr.resize(exAttrNames.size());

                for (int ii = 0; ii < exAttrNames.size(); ii++) {
                    for (auto &attr1: exAttrNames[ii]) {//scan through possible names
                        size_t pos1 = nowLine.find(" " + attr1 + " "); //attribute name is separated by spaces
                        if (pos1 != std::string::npos)
                            pos1 = nowLine.find_first_not_of(' ', pos1 + attr1.size() + 1);
                        if (pos1 != std::string::npos) {
                            exAttr[ii] = nowLine.substr(pos1, nowLine.find_first_of(' ', pos1) - pos1);
                        }
                    }
                }

                if (exAttr[0].empty()) {//no transcript ID
                    logFile_ << "WARNING: while processing pGe.sjdbGTFfile=" << gtfFile
                             << ": no transcript_id for line:\n";
                    logFile_ << oneLine << "\n";
                    logFile_.flush();
                    exAttr[0] = "tr_" + chr1 + "_" + std::to_string(exonStart) + "_" + std::to_string(exonEnd) + "_" +
                                std::to_string(exonN); //unique name for the transcript
                }

                if (exAttr[1].empty()) {//no gene ID
                    logFile_ << "WARNING: while processing pGe.sjdbGTFfile=" << gtfFile << ": no gene_id for line:\n";
                    logFile_ << oneLine << "\n";
                    logFile_.flush();
                    exAttr[1] = "MissingGeneID";
                }

                if (exAttr[2].empty()) {//no gene name
                    exAttr[2] = exAttr[1];
                }

                if (exAttr[3].empty()) {//no gene type
                    exAttr[3] = "MissingGeneType";
                }

                transcriptIDNumber.insert(std::pair<std::string, size_t>(exAttr[0], transcriptIDNumber.size()));
                if (transcriptID_.size() < transcriptIDNumber.size()) {
                    transcriptID_.push_back(exAttr[0]);
                    if (str1 == '+') transcriptStrand_.push_back(1);
                    else if (str1 == '-') transcriptStrand_.push_back(2);
                    else transcriptStrand_.push_back(0);
                }

                geneIDNumber.insert(std::pair<std::string, size_t>(exAttr[1], geneIDNumber.size()));
                if (geneID_.size() < geneIDNumber.size()) {//new gene is added
                    geneID_.push_back(exAttr[1]);
                    geneAttr_.emplace_back(exAttr[2], exAttr[3]);
                }

                exonLoci_[exonN] = {transcriptIDNumber[exAttr[0]], static_cast<int64_t>(exonStart +
                                                                                        genome.chromosomes_[genome.chromosomeNameToIndex_[chr1]].start -
                                                                                        1),
                                    static_cast<int64_t>(exonEnd +
                                                         genome.chromosomes_[genome.chromosomeNameToIndex_[chr1]].start -
                                                         1), geneIDNumber[exAttr[1]]};
                exonN++;

            }
        }

        if (exonN == 0) std::cerr << "Warning: no exons found in GTF file " << gtfFile << "\n";
        exonLoci_.resize(exonN);
    }


    int GTF::fillSjdbLoci(const std::string &dirOut, Genome &genome) {
        // from exonLoci_ to sjdbLoci_
        // return the number of junctions
        std::sort(exonLoci_.begin(), exonLoci_.end());
        size_t exonNum = exonLoci_.size();
        if (exonNum == 0) return 0;
        //todo output exon-gene info

        std::vector<exonTrLoci> extrLoci;
        extrLoci.resize(exonNum);

        size_t trExon = 0;
        for (size_t i = 0; i <= exonNum; ++i) {
            if (i == exonNum || exonLoci_[i].trID != exonLoci_[trExon].trID) {
                //process exons from trExon to i-1
                for (size_t j = trExon; j < i; ++j) {
                    extrLoci[j].trEnd = exonLoci_[i - 1].end;
                }
                if (i == exonNum) break;
                trExon = i;
            }
            extrLoci[i].trStart = exonLoci_[trExon].start;
            extrLoci[i].trID = exonLoci_[i].trID;
            extrLoci[i].exonStart = exonLoci_[i].start;
            extrLoci[i].exonEnd = exonLoci_[i].end;
            extrLoci[i].geneId = exonLoci_[i].geneId;
        }
        std::sort(extrLoci.begin(), extrLoci.end());

        //todo output exon transcript info

        std::vector<sjStride> sjLoci;
        sjLoci.reserve(exonNum);
        size_t trID = exonLoci_[0].trID;

        for (size_t i = 1; i < exonNum; ++i) {
            if (trID == exonLoci_[i].trID) {
                size_t chr = genome.getPosChrIndex(exonLoci_[i].start);
                if (exonLoci_[i].start <= exonLoci_[i - 1].end + 1) {
                    //touching or overlapping exons
                    continue;
                }
                sjLoci.emplace_back(exonLoci_[i - 1].end + 1, exonLoci_[i].start - 1, transcriptStrand_[trID],
                                    exonLoci_[i].geneId);
            } else {
                trID = exonLoci_[i].trID;
            }
        }

        std::sort(sjLoci.begin(), sjLoci.end());
        char strandChar[3] = {'.', '+', '-'};
        size_t sjNum = sjLoci.size();
        size_t prevSjNum = sjdbLoci_.size();
        sjdbLoci_.reserve(sjNum + prevSjNum);
        for (size_t i = 0; i < sjNum; ++i) {
            if (i == 0 || sjLoci[i] != sjLoci[i - 1]) {
                size_t chr = genome.getPosChrIndex(sjLoci[i].sjStart);
                sjdbLoci_.emplace_back(genome.chromosomes_[chr].name,
                                       sjLoci[i].sjStart - genome.chromosomes_[chr].start + 1,
                                       sjLoci[i].sjEnd - genome.chromosomes_[chr].start + 1,
                                       strandChar[sjLoci[i].sjStrand], 20);
                sjdbLoci_.back().gene.insert(sjLoci[i].geneId);
            } else {
                sjdbLoci_.back().gene.insert(sjLoci[i].geneId);
            }
        }

        // todo output sjdbList

        return (int) (sjdbLoci_.size() - prevSjNum);

    }

    void GTF::insertJunctions(Genome &genome, GenomeIndex &genomeIndex) {
        if (sjdbLoci_.empty()) return;
        size_t sjNum = sjdbLoci_.size();
        sjdbSeq_.resize(2 * sjdbLoci_.size() * config_.sjdbLength);
        auto *sjdbStart = new size_t[sjNum];
        auto *sjdbEnd = new size_t[sjNum];
        auto *sjdbMotif = new uint8_t[sjNum];
        auto *sjdbShiftLeft = new uint8_t[sjNum];
        auto *sjdbShiftRight = new uint8_t[sjNum];

        std::string prevChr;
        size_t iChr = 0;
        for (size_t i = 0; i < sjNum; ++i) {
            if (prevChr != sjdbLoci_[i].chr) {
                for (iChr = 0; iChr < genome.chromosomes_.size(); ++iChr) {
                    if (genome.chromosomes_[iChr].name == sjdbLoci_[i].chr) break;
                }
                if (iChr >= genome.chromosomes_.size()) {
                    // todo report error
                    std::cerr << "Warning: chromosome " << sjdbLoci_[i].chr
                              << " not found in genome, skipping junction\n";
                    continue;
                }
                prevChr = sjdbLoci_[i].chr;
            }

            sjdbStart[i] = sjdbLoci_[i].start + genome.chromosomes_[iChr].start - 1;
            sjdbEnd[i] = sjdbLoci_[i].end + genome.chromosomes_[iChr].start - 1;

            //todo can be reusable in stitching
            //judge the junction motif
            const char left1 = genome.sequence_[sjdbStart[i]];
            const char left2 = genome.sequence_[sjdbStart[i] + 1];
            const char right1 = genome.sequence_[sjdbEnd[i] - 1];
            const char right2 = genome.sequence_[sjdbEnd[i]];
            if (left1 == 'G' && left2 == 'T' && right1 == 'A' && right2 == 'G') sjdbMotif[i] = GTAG;
            else if (left1 == 'C' && left2 == 'T' && right1 == 'A' && right2 == 'C') sjdbMotif[i] = CTAC;
            else if (left1 == 'G' && left2 == 'C' && right1 == 'A' && right2 == 'G') sjdbMotif[i] = GCAG;
            else if (left1 == 'C' && left2 == 'T' && right1 == 'G' && right2 == 'C') sjdbMotif[i] = CTGC;
            else if (left1 == 'A' && left2 == 'T' && right1 == 'A' && right2 == 'C') sjdbMotif[i] = ATAC;
            else if (left1 == 'G' && left2 == 'T' && right1 == 'A' && right2 == 'T') sjdbMotif[i] = GTAT;
            else sjdbMotif[i] = NON_CANONICAL;

            //calc repeat length
            size_t jjL = 0;
            size_t jjR = 0;

            while (jjL < sjdbStart[i] - 1
                   && genome.sequence_[sjdbStart[i] - 1 - jjL] == genome.sequence_[sjdbEnd[i] - jjL]
                   && genome.sequence_[sjdbStart[i] - 1 - jjL] != '#'
                   && genome.sequence_[sjdbStart[i] - 1 - jjL] != 'N'
                   && jjL < 255) {
                jjL++;
            }

            while (sjdbEnd[i] + jjR < genome.sequence_.size()
                   && genome.sequence_[sjdbStart[i] + jjR] == genome.sequence_[sjdbEnd[i] + 1 + jjR]
                   && genome.sequence_[sjdbStart[i] + jjR] != '#'
                   && genome.sequence_[sjdbStart[i] + jjR] != 'N'
                   && jjR < 255) {
                jjR++;
            }

            sjdbShiftLeft[i] = (uint8_t) jjL;
            sjdbShiftRight[i] = (uint8_t) jjR;

            //todo report too long repeats

            sjdbStart[i] -= jjL;
            sjdbEnd[i] -= jjL;
        }

        //todo can be improved
        std::vector<sjdbSortRecord> sjdbSortRecords;
        sjdbSortRecords.resize(sjNum);
        for (size_t i = 0; i < sjNum; ++i) {
            size_t shift = 0;
            switch (sjdbLoci_[i].strand) {
                case '+':
                    shift = 0;
                    break;
                case '-':
                    shift = genome.sequence_.length();
                    break;
                default:
                    shift = genome.sequence_.length() * 2;
            }
            sjdbSortRecords[i] = {sjdbStart[i] + shift, sjdbEnd[i] + shift, i};
        }
        std::sort(sjdbSortRecords.begin(), sjdbSortRecords.end());

        std::vector<size_t> sjdbNewIndex;
        sjdbNewIndex.reserve(sjNum);

        for (size_t i = 0; i < sjNum; ++i) {
            size_t prevIndex;
            size_t nowIndex = sjdbSortRecords[i].originalIndex;
            size_t nowSjNum = sjdbNewIndex.size();
            if (nowSjNum > 0) {
                prevIndex = sjdbNewIndex.back();
            }

            if (nowSjNum == 0 || sjdbStart[nowIndex] != sjdbStart[prevIndex]
                || sjdbEnd[nowIndex] != sjdbEnd[prevIndex]) {
                sjdbNewIndex.push_back(nowIndex);
            } else if (sjdbLoci_[nowIndex].priority < sjdbLoci_[prevIndex].priority) {
                continue;
            } else if (sjdbLoci_[nowIndex].priority > sjdbLoci_[prevIndex].priority) {
                sjdbNewIndex.back() = nowIndex;
            } else if ((sjdbMotif[nowIndex] != NON_CANONICAL && sjdbMotif[prevIndex] == NON_CANONICAL) ||
                       ((sjdbMotif[nowIndex] != NON_CANONICAL) == (sjdbMotif[prevIndex] != NON_CANONICAL) &&
                        sjdbShiftLeft[nowIndex] < sjdbShiftLeft[prevIndex])) {
                sjdbNewIndex.back() = nowIndex;
            }
        }
        sjdbSortRecords.reserve(sjdbNewIndex.size());
        for (size_t i = 0; i < sjdbNewIndex.size(); ++i) {
            size_t shift = sjdbMotif[sjdbNewIndex[i]] == NON_CANONICAL ? 0 : sjdbShiftLeft[sjdbNewIndex[i]];
            sjdbSortRecords[i] = {
                    sjdbStart[sjdbNewIndex[i]] + shift,
                    sjdbEnd[sjdbNewIndex[i]] + shift,
                    sjdbNewIndex[i]
            };
        }

        std::sort(sjdbSortRecords.begin(), sjdbSortRecords.end());

        genome.sjdb.reserve(sjdbNewIndex.size());

        for (size_t i = 0; i < sjdbNewIndex.size(); ++i) {
            size_t nowIndex = sjdbSortRecords[i].originalIndex;
            bool replace = false;

            size_t nowSjNum = genome.sjdb.size();
            if (nowSjNum > 0 && genome.sjdb.back().start == sjdbStart[nowIndex] &&
                genome.sjdb.back().end == sjdbEnd[nowIndex]) {
                size_t prevIndex = sjdbSortRecords[i - 1].originalIndex;
                if (sjdbLoci_[nowIndex].priority < sjdbLoci_[prevIndex].priority) {
                    continue;
                } else if (sjdbLoci_[nowIndex].priority > sjdbLoci_[prevIndex].priority) {
                    replace = true;
                } else if (genome.sjdb.back().strand > 0 && sjdbLoci_[nowIndex].strand == '.') {
                    continue;
                } else if (genome.sjdb.back().strand == 0 && sjdbLoci_[nowIndex].strand != '.') {
                    replace = true;
                } else if (genome.sjdb.back().strand > 0 && sjdbMotif[nowIndex] == NON_CANONICAL) {
                    genome.sjdb.back().strand = 0;
                    continue;
                } else if ((genome.sjdb.back().motif != NON_CANONICAL && sjdbMotif[nowIndex] == NON_CANONICAL) ||
                           (genome.sjdb.back().motif % 2 == (2 - genome.sjdb.back().strand))) {
                    continue;
                } else {
                    replace = true;
                }
            }
            uint8_t strand = 0;
            switch (sjdbLoci_[i].strand) {
                case '+':
                    strand = 1;
                    break;
                case '-':
                    strand = 2;
                    break;
                default:
                    strand = 0;
            }
            if (replace) {
                if (strand == 0) {
                    if (genome.sjdb.back().motif != NON_CANONICAL) {
                        strand = 2 - (genome.sjdb.back().motif % 2);
                    }
                }

                genome.sjdb.back() = {
                        .start = sjdbSortRecords[i].start,
                        .end = sjdbSortRecords[i].end,
                        .motif = sjdbMotif[nowIndex],
                        .shiftLeft = sjdbShiftLeft[nowIndex],
                        .shiftRight = sjdbShiftRight[nowIndex],
                        .strand = strand
                };
            } else {
                genome.sjdb.push_back({
                                              .start = sjdbSortRecords[i].start,
                                              .end = sjdbSortRecords[i].end,
                                              .motif = sjdbMotif[nowIndex],
                                              .shiftLeft = sjdbShiftLeft[nowIndex],
                                              .shiftRight = sjdbShiftRight[nowIndex],
                                              .strand = strand
                                      });
            }
        }
        genome.sjdbNum = genome.sjdb.size();
        genome.sjDonorStart_.reserve(genome.sjdbNum);
        genome.sjAcceptorStart_.reserve(genome.sjdbNum);

        //todo output final sjdb info

        size_t sjGStart = 0;
        for (size_t i = 0; i < genome.sjdbNum; ++i) {
            genome.sjDonorStart_[i] = genome.sjdb[i].start - config_.sjdbOverhang;
            genome.sjAcceptorStart_[i] = genome.sjdb[i].end + 1;

            if (genome.sjdb[i].motif == NON_CANONICAL) {
                genome.sjDonorStart_[i] += genome.sjdb[i].shiftLeft;
                genome.sjAcceptorStart_[i] += genome.sjdb[i].shiftLeft;
            }

            memcpy((void *) (sjdbSeq_.c_str() + sjGStart), genome.sequence_.c_str() + genome.sjDonorStart_[i],
                   config_.sjdbOverhang);
            memcpy((void *) (sjdbSeq_.c_str() + sjGStart + config_.sjdbOverhang),
                   genome.sequence_.c_str() + genome.sjAcceptorStart_[i], config_.sjdbOverhang);
            memset((void *) (sjdbSeq_.c_str() + sjGStart + config_.sjdbOverhang * 2),'#',SJDB_PADDING_LENGTH);
            sjGStart += config_.sjdbLength;
        }

        delete[] sjdbStart;
        delete[] sjdbEnd;
        delete[] sjdbMotif;
        delete[] sjdbShiftLeft;
        delete[] sjdbShiftRight;

        modifyIndex(genome, genomeIndex);

    }

    inline char complement(char base) {
        switch (base) {
            case 'A':
                return 'T';
            case 'T':
                return 'A';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            default:
                return base; // For 'N' or any other character
        }
    }


    struct insertRecord {
        size_t pos;
        size_t sjPos;
    };

    class insertRecordComparator {
    private:
        const char *sjdbSeq_;
    public:
        explicit insertRecordComparator(const char *sjdbSeq) : sjdbSeq_(sjdbSeq) {}

        bool operator()(const insertRecord &a, const insertRecord &b) const {
            if (a.pos != b.pos) return a.pos < b.pos;
            const char *seqA = sjdbSeq_ + a.sjPos;
            const char *seqB = sjdbSeq_ + b.sjPos;
            int i = 0;
            while (true) {
                if (seqA[i] != seqB[i]) {
                    if (seqA[i] == 'N' || seqA[i] == '#') return false;
                    if (seqB[i] == 'N' || seqB[i] == '#') return true;
                    return seqA[i] < seqB[i];
                }
                if (seqA[i] == '#') return a.sjPos < b.sjPos;
                i++;
            }
        }
    };


    void GTF::modifyIndex(Genome &genome, GenomeIndex &genomeIndex) {
        //change SA and index
        int64_t sjSeqLength = config_.sjdbLength * genome.sjdbNum;
        sjdbSeqLengthSingleStrand_ = sjSeqLength;
        for (int i = 0; i < sjSeqLength; ++i) {
            sjdbSeq_[2 * sjSeqLength - 1 - i] = complement(sjdbSeq_[i]);
        }
        sjdbSeq_ += std::string('#',SJDB_PADDING_LENGTH); //padding to avoid overflow

        // find insert positions in SA
        std::vector<insertRecord> insertPos;
        insertPos.resize(2 * genome.sjdbNum * config_.sjdbLength);
#pragma omp parallel num_threads(8)
#pragma omp for schedule (dynamic, 1000) reduction(+:sjNew)
        for (size_t i = 0; i < 2 * genome.sjdbNum; ++i) {
            for (size_t sjStart = 0; sjStart < config_.sjdbLength; ++sjStart) {
                size_t ind = i * config_.sjdbLength + sjStart;
                if (sjdbSeq_[i * config_.sjdbLength + sjStart] == '#' ||
                    sjdbSeq_[i * config_.sjdbLength + sjStart] == 'N') {
                    insertPos[ind].pos = -1;
                } else {
                    insertPos[ind].pos = genomeIndex.suffixArray.findInsertPosition(
                            sjdbSeq_.substr(i * config_.sjdbLength + sjStart),
                            genome.sequence_);
                    insertPos[ind].sjPos = i * config_.sjdbLength + sjStart;
                }

            }
        }
        int trueIndNum = 0;
        for (size_t i = 0; i < insertPos.size(); ++i) {
            if (insertPos[i].pos != -1) {
                insertPos[trueIndNum++] = insertPos[i];
            }
        }

        std::sort(insertPos.begin(), insertPos.begin() + trueIndNum, insertRecordComparator(sjdbSeq_.c_str()));

        insertPos[trueIndNum].pos = MARK_END;

        // modify genome and index


        // positive genome strand - negative genome strand - #### - positive sj - negative sj
        size_t originalGenomeLength = genomeIndex.genome->originalGenomeLength;

        int64_t nowInsertSjIndex = 0;
        int64_t nowSAIndex = 0;

        genomeIndex.suffixArrayFirstPass.sa_.buildFromPackedArrayExtendForward(
                genomeIndex.suffixArray.sa_, trueIndNum);
        genomeIndex.suffixArrayFirstPass.sa_.isOwner = false;


        for (size_t i = 0; i < genomeIndex.suffixArray.length(); ++i) {
            while (i == insertPos[nowInsertSjIndex].pos) {
                size_t sjPos = insertPos[nowInsertSjIndex].sjPos;

                sjPos += originalGenomeLength;

                genomeIndex.suffixArrayFirstPass.sa_.set(nowSAIndex, sjPos);
                ++nowInsertSjIndex;
                ++nowSAIndex;

            }
            size_t nowSAValue = genomeIndex.suffixArray.sa_.get(i);
            genomeIndex.suffixArrayFirstPass.sa_.set(nowSAIndex, nowSAValue);
            ++nowSAIndex;
        }

        for (; nowInsertSjIndex < trueIndNum; ++nowInsertSjIndex) {
            size_t sjPos = insertPos[nowInsertSjIndex].sjPos;

            sjPos += originalGenomeLength;

            genomeIndex.suffixArrayFirstPass.sa_.set(nowSAIndex, sjPos);
            ++nowSAIndex;
        }
        size_t saLength = genomeIndex.suffixArray.length() + trueIndNum;
        genomeIndex.suffixArray.sa_.buildFromPackedArrayExtendForward(genomeIndex.suffixArrayFirstPass.sa_, 0);

        genomeIndex.suffixArray.sa_.setLength(saLength);

        genome.sequence_ += sjdbSeq_;


        // recalculate LCP
        std::vector<uint16_t> &LCP = genomeIndex.longestCommonPrefix_;

        LCP.resize(saLength, 0);
        PackedArray rk;
        PackedArray &sa = genomeIndex.suffixArray.sa_;
        rk.init(genome.sequence_.length(), sa.getWordLengthBits());
        for (size_t i = 0; i < genome.sequence_.length(); ++i) rk.set(i, 0);
        for (size_t i = 0; i < saLength; ++i) rk.set(sa.get(i), i);
        size_t k = 0;
        for (int64_t i = 0; i < genome.sequence_.length(); ++i) {
            size_t j = rk.get(i);
            if (j == 0) {
                k = 0;
                continue;
            }
            if (k > 0) k--;
            while (genome.sequence_[i + k] == genome.sequence_[sa.get(j - 1) + k]) k++;
            if (k < std::numeric_limits<uint16_t>::max()) LCP[j] = k;
            else LCP[j] = std::numeric_limits<uint16_t>::max();
        }





        // update primary index and secondary index
        const int MER_LENGTH = genomeIndex.MER_LENGTH;
        std::vector<sjHash> sjHashRecord;
        sjHashRecord.resize(sjSeqLength * 2, {-1, 0});


        for (int i = 0; i < sjSeqLength * 2; ++i) {
            sjHashRecord[i] = encodeKMer(sjdbSeq_.substr(i), MER_LENGTH);

        }

        std::sort(sjHashRecord.begin(), sjHashRecord.end());

        int64_t nowShift = 0;
        int32_t prevHashInsert = -1;
        int32_t nowHashInsert = 0;

        for (size_t i = 0; i < sjSeqLength * 2; ++i) {
            if (sjHashRecord[i].hash == -1) continue;
            nowHashInsert = sjHashRecord[i].hash;
            int32_t length = sjHashRecord[i].length;
            if (length == 0) continue; // should not happen
            if (nowHashInsert != prevHashInsert) {
                if (prevHashInsert != -1) {
                    for (int32_t j = prevHashInsert + 1; j <= nowHashInsert; ++j) {
                        genomeIndex.patternMerMap_[j].leftSAIndex += nowShift;
                    }
                    //rebuild secondary index for prevHashInsert
                    genomeIndex.buildKMerMapSingle(prevHashInsert);

                }
                prevHashInsert = nowHashInsert;
            }
            ++nowShift;
            if (length == MER_LENGTH) {
                genomeIndex.patternMerMap_[nowHashInsert].upperRange++;
            }
            int32_t originalLength = genomeIndex.patternMerMap_[nowHashInsert].length;
            if (originalLength < length) {
                genomeIndex.patternMerMap_[nowHashInsert].length = length;
            }

        }

        if (prevHashInsert != -1) {
            for (int32_t j = prevHashInsert + 1; j < genomeIndex.MER_NUM; ++j) {
                genomeIndex.patternMerMap_[j].leftSAIndex += nowShift;
            }
            // rebuild secondary index for prevHashInsert
            genomeIndex.buildKMerMapSingle(prevHashInsert);
        }


    }


}