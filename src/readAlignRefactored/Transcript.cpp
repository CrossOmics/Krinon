#include <charconv>
#include <algorithm>
#include "Transcript.h"
#include <iostream>

namespace RefactorProcessing {
    /**
     * Append a single integer onto the string, followed by given trailing
     * character.
     */
    inline void appendIntAndChar(std::string& str, int64_t val, char trailingChar) {
        char buf[32];
        auto [ptr, ec] = std::to_chars(buf, buf + sizeof(buf), val);
        *ptr = trailingChar;
        str.append(buf, ptr - buf + 1);
    }

    std::string Transcript::getCIGAR() const {
        std::string cigar;
        if (readStart > 0) {
            // soft clipping at the start
            appendIntAndChar(cigar, readStart, 'S');
        }
        for (size_t i = 0; i < exons.size(); ++i) {
            // exon
            appendIntAndChar(cigar, exons[i].length, 'M');
            if (i < exons.size() - 1) {
                if (sj[i].type >= 0) {
                    // splice junction
                    appendIntAndChar(cigar, sj[i].length, 'N');
                } else if (sj[i].type == -1) {
                    // deletion
                    appendIntAndChar(cigar, sj[i].length, 'D');
                } else if (sj[i].type == -2) {
                    // insertion
                    appendIntAndChar(cigar, -sj[i].length, 'I');
                }
            }
        }
        if (exons.back().readStart + exons.back().length < readLength) {
            // soft clipping at the end
            appendIntAndChar(cigar, readLength - (exons.back().readStart + exons.back().length), 'S');
        }

        return cigar;
    }

    void Transcript::convertToSAM(Read& r, bool isPaired, int mulmapNum, std::string& out) const {
        int qualityScore = 255;
        if (mulmapNum == 2) qualityScore = 3;
        if (mulmapNum == 3) qualityScore = 1;
        if (mulmapNum > 3) qualityScore = 0;

        if (!isPaired) {
            out.append(r.name).append("\t");
            appendIntAndChar(out, strand == 0 ? 0 : 16, '\t');
            out.append(chr).append("\t");
            appendIntAndChar(out, posInChr + 1, '\t');
            appendIntAndChar(out, qualityScore, '\t');
            out.append(CIGAR).append("\t*\t0\t0\t");
            out.append(strand == 0 ? r.sequence[0] : r.sequence[1]).append("\t");
            out.append(r.quality).append("\n");
            return;
        }

        // std::string read1output, read2output;
        int Flag1 = (strand == 0) ? 99 : 163;
        int Flag2 = (strand == 0) ? 147 : 83;

        // Once again, we pre-allocate these (plus TODO: for Arvin ...)
        std::string cigar1, cigar2;
        cigar1.reserve(1024);
        cigar2.reserve(1024);

        size_t mate2StartExon = 0;
        for (size_t iExon = 0; iExon < exons.size(); ++iExon) {
            // exon
            appendIntAndChar(cigar1, exons[iExon].length, 'M');
            if (iExon < exons.size() - 1) {
                if (sj[iExon].type >= 0) {
                    // splice junction
                    appendIntAndChar(cigar1, sj[iExon].length, 'N');
                } else if (sj[iExon].type == -1) {
                    // deletion
                    appendIntAndChar(cigar1, sj[iExon].length, 'D');
                } else if (sj[iExon].type == -2) {
                    // insertion
                    appendIntAndChar(cigar1, -sj[iExon].length, 'I');
                } else if (sj[iExon].type == -3) {
                    int64_t length = (strand == 0) ? r.mate1Length : r.mate2Length;
                    int64_t offset = exons[iExon].readStart + exons[iExon].length;
                    if (offset < length) {
                        appendIntAndChar(cigar1, length - offset, 'S');
                    }
                    mate2StartExon = iExon + 1;
                    break;
                }
            }
        }

        for (size_t iExon = mate2StartExon; iExon < exons.size(); ++iExon) {
            // exon
            appendIntAndChar(cigar2, exons[iExon].length, 'M');
            if (iExon < exons.size() - 1) {
                if (sj[iExon].type >= 0) {
                    // splice junction
                    appendIntAndChar(cigar2, exons[iExon].length, 'N');
                } else if (sj[iExon].type == -1) {
                    // deletion
                    appendIntAndChar(cigar2, exons[iExon].length, 'D');
                } else if (sj[iExon].type == -2) {
                    // insertion
                    appendIntAndChar(cigar2, exons[iExon].length, 'I');
                }
            }
        }

        if (exons[0].readStart > 0) {
            std::string prefix;
            appendIntAndChar(prefix, exons[0].readStart, 'S');
            cigar1.insert(0, prefix);
        }

        if (exons.back().readStart + exons.back().length < r.length){
            appendIntAndChar(cigar2, readLength - (exons.back().readStart + exons.back().length), 'S');
        }

        int64_t mateLength = (strand == 0) ? r.mate1Length + 1 : r.mate2Length + 1;
        int64_t exonStart2 = exons[mate2StartExon].readStart;
        if (exonStart2 > mateLength) {
            std::string prefix;
            appendIntAndChar(prefix, exonStart2 - mateLength, 'S');
            cigar2.insert(0, prefix);
        }
        
        int64_t distance = (exons.back().genomeStart + exons.back().length) - exons[0].genomeStart;
        size_t mate2StartInChr = exons[mate2StartExon].genomeStart - (exons[0].genomeStart - posInChr) + 1;

        std::string_view readSeq = (strand == 0) ? r.sequence[0] : r.sequence[1];
        std::string_view readQual = r.quality;
        std::string_view readSeq1, readSeq2, readQual1, readQual2;
        int64_t seq1SubstringEnd = (strand == 0) ? r.mate1Length : r.mate2Length;
        int64_t seq2SubstringEnd = (strand == 0) ? r.mate2Length : r.mate1Length;

        // Here, we defer reversing the read qualities until we want to acutally write!
        readSeq1 = readSeq.substr(0, seq1SubstringEnd);
        readSeq2 = readSeq.substr(seq1SubstringEnd + 1, seq2SubstringEnd);
        readQual1 = readQual.substr(0, seq1SubstringEnd);
        readQual2 = readQual.substr(seq1SubstringEnd + 1, seq2SubstringEnd);

        // Output for the first line ...
        out.append(r.name).append("\t");
        appendIntAndChar(out, Flag1, '\t');
        out.append(chr).append("\t");
        appendIntAndChar(out, posInChr + 1, '\t');
        appendIntAndChar(out, qualityScore, '\t');
        out.append(cigar1).append("\t = \t");
        appendIntAndChar(out, mate2StartInChr, '\t');
        out.append("\t").append(readSeq1).append("\t");
        // NOW, add the reversed quality reading if needed
        if (strand == 0) {
            out.append(readQual1);
        }
        else {
            out.append(readQual1.rbegin(), readQual1.rend());
        }
        out.append("\n");
        // Output for the second line ...
        out.append(r.name).append("\t");
        appendIntAndChar(out, Flag2, '\t');
        out.append(chr).append("\t");
        appendIntAndChar(out, mate2StartInChr, '\t');
        appendIntAndChar(out, qualityScore, '\t');
        out.append(cigar2).append("\t = \t");
        appendIntAndChar(out, posInChr + 1, '\t');
        appendIntAndChar(out, -distance, '\t');
        out.append(readSeq2).append("\t");
        if (strand == 0) {
            out.append(readQual2);
        } else {
            out.append(readQual2.rbegin(), readQual2.rend());
        }
        out.append("\n");

        return;
    }
}
