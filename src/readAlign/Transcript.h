#ifndef RNAALIGNREFACTORED_TRANSCRIPT_H
#define RNAALIGNREFACTORED_TRANSCRIPT_H
#include "StitchingUtils.h"
#include "../utils/types.h"
#include <cstdint>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
namespace rna {




    struct Exon {
        int64_t start{0};
        int64_t length{0};
        int64_t readStart{0};
        int64_t score{0};
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
}

#endif //RNAALIGNREFACTORED_TRANSCRIPT_H
