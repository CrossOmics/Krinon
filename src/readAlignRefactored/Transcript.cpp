#include "Transcript.h"

namespace RefactorProcessing {
std::string Transcript::getCIGAR() const {
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
            }
        }
    }
    if (exons.back().readStart + exons.back().length < readLength) {
        cigar += std::to_string(readLength - (exons.back().readStart + exons.back().length)) +
                 "S"; // soft clipping at the end
    }

    return cigar;
}
}
