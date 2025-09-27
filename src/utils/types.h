#ifndef RNAALIGNMENT_TYPES_H
#define RNAALIGNMENT_TYPES_H

#include <cstdint>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#define READ_PAIRED 0x1
#define READ_MAPPED_IN_PROPER_PAIR 0x2
#define READ_UNMAPPED 0x4
#define READ_MATE_UNMAPPED 0x8
#define READ_REVERSE_STRAND 0x10
#define READ_MATE_REVERSE_STRAND 0x20
#define READ_FIRST_IN_PAIR 0x40
#define READ_SECOND_IN_PAIR 0x80

namespace rna {
    using GenomePos = int64_t;
    using ReadPos = int64_t;
    using Score = int64_t;
    using Length = int64_t;
    struct Read {
        std::string sequence[2]; // 0 for forward, 1 for complementary
        std::string quality;
        std::string name;
        int64_t length;
        int64_t mate1Length; // for paired-end reads
        int64_t mate2Length;
    };

    struct MatchResultStored {
        uint8_t length;
        int32_t upperRange;
        int64_t leftSAIndex;

        MatchResultStored() : length(0), upperRange(-1), leftSAIndex(-1) {};
    };

    struct MatchStats {
        int64_t maxLength;
        size_t count;
        size_t LongestMatchLeftPosInSA;
        size_t LongestMatchRightPosInSA;

        MatchStats(int64_t len = 0, size_t cnt = 0, size_t l = 0, size_t r = 0) : maxLength(len), count(cnt),
                                                                                  LongestMatchLeftPosInSA(l),
                                                                                  LongestMatchRightPosInSA(r) {}
    };

    struct Split {
        ReadPos readStart{0};
        size_t length{0};
        int iFragment{0};
    };

    struct Align {
        ReadPos readStart{0};
        int64_t length{0};
        size_t leftSAIndex{0};
        size_t rightSAIndex{0};
        size_t rep{0};
        int direction{0};
        bool isAnchor{false};
        int iFragment{0};

        bool operator<(const Align &other) const {
            return readStart < other.readStart;
        }
    };

    struct WindowAlign {
        ReadPos readStart{0};
        GenomePos genomeStart{0};
        int64_t length{0};
        int64_t score{0};
        bool isAnchor{false};
        int64_t isj{-1};// annotation index
        int iFragment{0};

        bool operator<(const WindowAlign &other) const {
            if (genomeStart != other.genomeStart) return genomeStart < other.genomeStart;
            return length < other.length;
        }
    };

    struct Window {
        int chrIndex{0};
        int direction{0};
        int numAnchors{0};
        int64_t startBin{0};
        int64_t endBin{0};
        WindowAlign *aligns{nullptr};
        int alignNum{0};
        int minLengthWhenFull{0};
    };


    struct SpliceJunction {

        int64_t type{0};//canon
        int64_t length{0};
        bool isAnnotated{false};
        int64_t shiftLeft{0};
        int64_t shiftRight{0};
        int64_t sjStrand{0};
    };

    // the change of the score, matches, mismatches, exons length of one single stitching step (stitching WA[j] to WA[i])
    struct StitchingRecord {
        enum Type {
            CANNOT_STITCH,
            PERFECT_MATCH,
            SAME_GAP,
            DELETION,
            INSERTION,
            SPLICE_JUNCTION,
            CROSS_FRAGMENTS
        } type{CANNOT_STITCH};
        int64_t score{0}; // including the gap penalty
        int matches{0};
        int mismatches{0};
        int64_t formerExonLengthShift{0};
        int64_t latterExonLengthShift{0};

        SpliceJunction stitchingType{0, 0, false};
    };

    struct ExtensionRecord {
        //score can be calculated by matches * MATCH_SCORE - i * MATCH_SCORE
        struct singleExtensionRecord {
            int length{0}; // extension length
            int matched{0}; // number of matched bases, in case of the appearance of 'N'
            int mismatches{
                    0}; // number of mismatches, since the best extension may not reach the max mismatch tolerance
        };
        // array of length [max_mismatch], to store the max extension length when allowing a certain number of mismatches
        // Note: must free the memory after use
        singleExtensionRecord *maxExtensionLengthWithMismatch{nullptr};
        int maxMismatch{0};
        int maxExtensionScore{0};

    };

    // selected pieces and corresponding stitching records
    // can be used to reconstruct the transcript
    // must ensure that the score is right (same as the transcript's final score)
    struct RawTranscript {
        int previousTranscriptId{-1};//id in the allWindowAligns_ array
        int newAlignId{-1}; //
        int64_t score{0}; // true score
        int mismatches{0}; // number of mismatches in the transcript
        int matches{0};
        int exonCount{0}; // number of exons in the transcript
        int extendedLengthForward{0}; //extension 5' end
        int extendedLengthBackward{0}; // extension 3' end
        int startAlignId{-1};
        bool crossFragments{false};
        int extendedLengthFragmentFormer{0};
        int extendedLengthFragmentLatter{0};
        int64_t firstFragmentMatchEnd{0};


    };

    struct SJDBOutput {
        std::string chr;
        int64_t left;
        int64_t right;
        bool annotated;
        int64_t type;
        int strand{0}; // 0 undefined, 1 forward, 2 reverse
    };


    struct Exon {
        GenomePos start{0};
        int64_t length{0};
        ReadPos readStart{0};
        Score score{0};
    };


    struct Transcript {
        std::string readName;
        std::string chr;
        std::string CIGAR;
        int strand{0};
        int64_t matched{0};
        int64_t unmatched{0};
        int64_t nIns{0};
        int64_t nDel{0};
        std::vector<Exon> exons;
        std::vector<SpliceJunction> sj;
        std::vector<WindowAlign> aligns;//Actually no need to store this, for debugging only
        int64_t readStart{0};
        int64_t genomeStart{0};
        int64_t posInChr{0};
        int64_t score{0};
        int64_t readLength{0};
        bool isPaired{false};
        int64_t chrStartPos{0};

        std::string getCIGAR() const {
            std::string cigar;
            if (readStart > 0) {
                cigar += std::to_string(readStart) + "S"; // soft clipping at the start
            }
            for (int64_t i = 0; i < exons.size(); ++i) {
                cigar += std::to_string(exons[i].length) + "M"; // exon
                if (i < exons.size() - 1) {
                    if (sj[i].type >= 0) {
                        cigar += std::to_string(sj[i].length) + "N"; // splice junction
                    } else if (sj[i].type == -1) {
                        cigar += std::to_string(sj[i].length) + "D"; // deletion
                    } else if (sj[i].type == -2) {
                        cigar += std::to_string(sj[i].length) + "I"; // insertion
                    } else if (sj[i].type == -3) {
                        cigar += "\n";
                        cigar += std::to_string(exons[i + 1].start - exons[i].start) + "d";
                        cigar += "\n"; //cross fragments
                    }
                }
            }
            if (exons.back().readStart + exons.back().length < readLength) {
                cigar += std::to_string(readLength - (exons.back().readStart + exons.back().length)) +
                         "S"; // soft clipping at the end
            }

            return cigar;
        }

        std::string outputSam(Read &read, int mulmapNum) const {
            int qualityScore = 255;
            if (mulmapNum == 2) qualityScore = 3;
            if (mulmapNum > 3) qualityScore = 1;
            if (!isPaired) {
                return read.name + "\t" +
                       std::to_string(strand == 0 ? 0 : 16) + "\t" +
                       chr + "\t" +
                       std::to_string(posInChr + 1) + "\t" + // SAM is 1-based
                       std::to_string(qualityScore) + "\t" +
                       CIGAR + "\t*\t0\t0\t" +
                       (strand == 0 ? read.sequence[0] : read.sequence[1]) + "\t" +
                       read.quality + "\n";
            } else {
                std::string read1output, read2output;
                int Flag1 = 0, Flag2 = 0;
                if (strand == 0) {
                    Flag1 = 99; // read1, mapped, forward, mate reverse
                    Flag2 = 147; // read2, mapped, reverse, mate forward
                } else {
                    Flag1 = 163; // read2, mapped, reverse, mate forward
                    Flag2 = 83; // read1, mapped, forward, mate reverse
                }
                std::string cigar1, cigar2;
                int64_t iExon = 0;
                int64_t mate2StartExon = 0;
                for (iExon = 0; iExon < exons.size(); ++iExon) {
                    cigar1 += std::to_string(exons[iExon].length) + "M"; // exon
                    if (iExon < exons.size() - 1) {
                        if (sj[iExon].type >= 0) {
                            cigar1 += std::to_string(sj[iExon].length) + "N"; // splice junction
                        } else if (sj[iExon].type == -1) {
                            cigar1 += std::to_string(sj[iExon].length) + "D"; // deletion
                        } else if (sj[iExon].type == -2) {
                            cigar1 += std::to_string(sj[iExon].length) + "I"; // insertion
                        } else if (sj[iExon].type == -3) {
                            if (strand == 0) {
                                if (exons[iExon].readStart + exons[iExon].length < read.mate1Length) {
                                    cigar1 += std::to_string(
                                            read.mate1Length - (exons[iExon].readStart + exons[iExon].length)) + "S";
                                }
                            } else {
                                if (exons[iExon].readStart + exons[iExon].length < read.mate2Length) {
                                    cigar1 += std::to_string(
                                            read.mate2Length - (exons[iExon].readStart + exons[iExon].length)) + "S";
                                }
                            }
                            mate2StartExon = iExon + 1;
                            break;
                        }
                    }
                }

                for (iExon = mate2StartExon;iExon < exons.size(); ++iExon) {
                    cigar2 += std::to_string(exons[iExon].length) + "M"; // exon
                    if (iExon < exons.size() - 1) {
                        if (sj[iExon].type >= 0) {
                            cigar2 += std::to_string(sj[iExon].length) + "N"; // splice junction
                        } else if (sj[iExon].type == -1) {
                            cigar2 += std::to_string(sj[iExon].length) + "D"; // deletion
                        } else if (sj[iExon].type == -2) {
                            cigar2 += std::to_string(sj[iExon].length) + "I"; // insertion
                        }
                    }
                }

                if (exons[0].readStart > 0){
                    cigar1 = std::to_string(exons[0].readStart) + "S" + cigar1;
                }
                if (exons.back().readStart + exons.back().length < read.length){
                    cigar2 += std::to_string(readLength - (exons.back().readStart + exons.back().length)) + "S";
                }
                if (strand == 0){
                    if (exons[mate2StartExon].readStart > read.mate1Length + 1){
                        cigar2 = std::to_string(exons[mate2StartExon].readStart - read.mate1Length - 1) + "S" + cigar2;
                    }
                }else{
                    if (exons[mate2StartExon].readStart > read.mate2Length + 1){
                        cigar2 = std::to_string(exons[mate2StartExon].readStart - read.mate2Length - 1) + "S" + cigar2;
                    }
                }
                int64_t distance = (exons.back().start + exons.back().length) - exons[0].start;
                std::string readSeq1 ,readSeq2, readQual1, readQual2;
                if (strand == 0) {
                    readSeq1 = read.sequence[0].substr(0, read.mate1Length);
                    readSeq2 = read.sequence[0].substr(read.mate1Length + 1, read.mate2Length);
                    readQual1 = read.quality.substr(0, read.mate1Length);
                    readQual2 = read.quality.substr(read.mate1Length + 1, read.mate2Length);
                } else {
                    readSeq1 = read.sequence[1].substr(0, read.mate2Length);
                    readSeq2 = read.sequence[1].substr(read.mate2Length + 1, read.mate1Length);
                    readQual1 = read.quality.substr(0, read.mate2Length);
                    readQual2 = read.quality.substr(read.mate2Length + 1, read.mate1Length);
                    std::reverse(readQual1.begin(),readQual1.end());
                    std::reverse(readQual2.begin(),readQual2.end());
                }
                int64_t mate2StartInChr = exons[mate2StartExon].start - chrStartPos + 1;

                read1output = read.name + "\t" + std::to_string(Flag1) + "\t" +
                              chr + "\t" +
                              std::to_string(posInChr + 1) + "\t" + // SAM is 1-based
                              std::to_string(qualityScore) + "\t"+ cigar1 + "\t = \t" +
                              std::to_string(mate2StartInChr) + "\t" +
                              std::to_string(distance) + "\t" + readSeq1 + "\t" + readQual1 + "\n";
                read2output = read.name + "\t" + std::to_string(Flag2) + "\t" +
                              chr + "\t" +
                              std::to_string(mate2StartInChr) + "\t" + // SAM is 1-based
                              std::to_string(qualityScore) + "\t"+ cigar2 + "\t = \t" +
                              std::to_string(posInChr + 1) + "\t" +
                              std::to_string(-distance) + "\t" + readSeq2 + "\t" + readQual2 + "\n";
                return read1output + read2output;

            }
        }
    };

    struct SAMEntry {
        std::string readName;
        int flag{0};
        std::string chr{"*"};
        int64_t pos{0};// 1-based
        int64_t mapq{0};
        std::string CIGAR;
        std::string RNEXT{"*"};
        int64_t PNEXT{0};
        int64_t TLEN{0};
        std::string seq;
        std::string qual;

        SAMEntry(Transcript &t, Read &read, int mapQ = 255) {
            readName = read.name;
            flag = t.strand == 0 ? 0 : 16; // 0 for forward, 16 for reverse
            chr = t.chr;
            pos = t.posInChr + 1; // SAM is 1-based
            mapq = mapQ;
            CIGAR = t.CIGAR;
            RNEXT = "*";
            PNEXT = 0;
            TLEN = 0; // not used in none paired-end reads
            seq = t.strand == 0 ? read.sequence[0] : read.sequence[1];
            qual = read.quality;
        }

        std::string toString() const {
            return readName + "\t" +
                   std::to_string(flag) + "\t" +
                   chr + "\t" +
                   std::to_string(pos) + "\t" +
                   std::to_string(mapq) + "\t" +
                   CIGAR + "\t" +
                   RNEXT + "\t" +
                   std::to_string(PNEXT) + "\t" +
                   std::to_string(TLEN) + "\t" +
                   seq + "\t" +
                   qual + "\n";
        }

        friend std::ostream &operator<<(std::ostream &os, const SAMEntry &entry) {
            os << entry.readName << "\t"
               << entry.flag << "\t"
               << entry.chr << "\t"
               << entry.pos << "\t"
               << entry.mapq << "\t"
               << entry.CIGAR << "\t"
               << entry.RNEXT << "\t"
               << entry.PNEXT << "\t"
               << entry.TLEN << "\t"
               << entry.seq << "\t"
               << entry.qual;
            return os;
        }

    };

    using ReadPtr = std::shared_ptr<Read>;
    using TranscriptPtr = std::shared_ptr<Transcript>;

}
#endif //RNAALIGNMENT_TYPES_H
