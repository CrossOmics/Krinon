#include <algorithm>
#include "Transcript.h"

namespace RefactorProcessing {
    std::string Transcript::getCIGAR() const {
        std::string cigar;
        if (readStart > 0) {
            cigar += std::to_string(readStart) + "S"; // soft clipping at the start
        }
        for (size_t i = 0; i < exons.size(); ++i) {
            cigar += std::to_string(exons[i].length) + "M"; // exon
            if (i < exons.size() - 1) {
                if (sj[i].type >= 0) {
                    cigar += std::to_string(sj[i].length) + "N"; // splice junction
                } else if (sj[i].type == -1) {
                    cigar += std::to_string(sj[i].length) + "D"; // deletion
                } else if (sj[i].type == -2) {
                    cigar += std::to_string(-sj[i].length) + "I"; // insertion
                }
            }
        }
        if (exons.back().readStart + exons.back().length < readLength) {
            cigar += std::to_string(readLength - (exons.back().readStart + exons.back().length)) +
                     "S"; // soft clipping at the end
        }

        return cigar;
    }

    std::string Transcript::convertToSAM(Read& r,bool isPaired, int mulmapNum) const {
        int qualityScore = 255;
        if (mulmapNum == 2) qualityScore = 3;
        if (mulmapNum == 3) qualityScore = 1;
        if (mulmapNum > 3) qualityScore = 0;
        if (!isPaired) {
            return r.name + "\t" +
                   std::to_string(strand == 0 ? 0 : 16) + "\t" +
                   chr + "\t" +
                   std::to_string(posInChr + 1) + "\t" + // SAM is 1-based
                   std::to_string(qualityScore) + "\t" +
                   CIGAR + "\t*\t0\t0\t" +
                   (strand == 0 ? r.sequence[0] : r.sequence[1]) + "\t" +
                   r.quality + "\n";
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
            size_t iExon = 0;
            int64_t mate2StartExon = 0;
            for (iExon = 0; iExon < exons.size(); ++iExon) {
                cigar1 += std::to_string(exons[iExon].length) + "M"; // exon
                if (iExon < exons.size() - 1) {
                    if (sj[iExon].type >= 0) {
                        cigar1 += std::to_string(sj[iExon].length) + "N"; // splice junction
                    } else if (sj[iExon].type == -1) {
                        cigar1 += std::to_string(sj[iExon].length) + "D"; // deletion
                    } else if (sj[iExon].type == -2) {
                        cigar1 += std::to_string(-sj[iExon].length) + "I"; // insertion
                    } else if (sj[iExon].type == -3) {
                        if (strand == 0) {
                            if (exons[iExon].readStart + exons[iExon].length < r.mate1Length) {
                                cigar1 += std::to_string(
                                        r.mate1Length - (exons[iExon].readStart + exons[iExon].length)) + "S";
                            }
                        } else {
                            if (exons[iExon].readStart + exons[iExon].length < r.mate2Length) {
                                cigar1 += std::to_string(
                                        r.mate2Length - (exons[iExon].readStart + exons[iExon].length)) + "S";
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
            if (exons.back().readStart + exons.back().length < r.length){
                cigar2 += std::to_string(readLength - (exons.back().readStart + exons.back().length)) + "S";
            }
            if (strand == 0){
                if (exons[mate2StartExon].readStart > r.mate1Length + 1){
                    cigar2 = std::to_string(exons[mate2StartExon].readStart - r.mate1Length - 1) + "S" + cigar2;
                }
            }else{
                if (exons[mate2StartExon].readStart > r.mate2Length + 1){
                    cigar2 = std::to_string(exons[mate2StartExon].readStart - r.mate2Length - 1) + "S" + cigar2;
                }
            }
            int64_t distance = (exons.back().genomeStart + exons.back().length) - exons[0].genomeStart;
            std::string readSeq1 ,readSeq2, readQual1, readQual2;
            if (strand == 0) {
                readSeq1 = r.sequence[0].substr(0, r.mate1Length);
                readSeq2 = r.sequence[0].substr(r.mate1Length + 1, r.mate2Length);
                readQual1 = r.quality.substr(0, r.mate1Length);
                readQual2 = r.quality.substr(r.mate1Length + 1, r.mate2Length);
            } else {
                readSeq1 = r.sequence[1].substr(0, r.mate2Length);
                readSeq2 = r.sequence[1].substr(r.mate2Length + 1, r.mate1Length);
                readQual1 = r.quality.substr(0, r.mate2Length);
                readQual2 = r.quality.substr(r.mate2Length + 1, r.mate1Length);
                std::reverse(readQual1.begin(),readQual1.end());
                std::reverse(readQual2.begin(),readQual2.end());
            }
            int64_t mate2StartInChr = exons[mate2StartExon].genomeStart - (exons[0].genomeStart - posInChr) + 1;

            read1output = r.name + "\t" + std::to_string(Flag1) + "\t" +
                          chr + "\t" +
                          std::to_string(posInChr + 1) + "\t" + // SAM is 1-based
                          std::to_string(qualityScore) + "\t"+ cigar1 + "\t = \t" +
                          std::to_string(mate2StartInChr) + "\t" +
                          std::to_string(distance) + "\t" + readSeq1 + "\t" + readQual1 + "\n";
            read2output = r.name + "\t" + std::to_string(Flag2) + "\t" +
                          chr + "\t" +
                          std::to_string(mate2StartInChr) + "\t" + // SAM is 1-based
                          std::to_string(qualityScore) + "\t"+ cigar2 + "\t = \t" +
                          std::to_string(posInChr + 1) + "\t" +
                          std::to_string(-distance) + "\t" + readSeq2 + "\t" + readQual2 + "\n";
            return read1output + read2output;

        }
    }

}
