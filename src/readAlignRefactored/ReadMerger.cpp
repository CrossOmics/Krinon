#include "ReadMerger.h"

namespace RefactorProcessing {
    void ReadMerger::setParam(const RefactorProcessing::Parameters &P) {
        isPaired_ = P.isPaired;
        maxMismatch_ = P.maxMismatch;
    }

    //brute-force merging
    int calcBestOverlap(const std::string& seq1, const std::string& seq2, int maxMismatch,double maxMismatchRate) {
        int bestOverlap = -1;
        int bestMatches = -1;
        for (size_t i = 0; i < seq1.length(); ++i){
            int mismatch = 0;
            int overlapLength = std::min(seq1.length() - i, seq2.length());
            int mismatchThreshold = std::min(maxMismatch, int(overlapLength * maxMismatchRate));

            for (size_t j = 0; j < seq2.length(); ++j){
                if (i + j >= seq1.length()) break;
                if (seq1[i + j] != seq2[j]) {
                    ++mismatch;
                    if (mismatch > mismatchThreshold) break;
                }
            }
            if (mismatch > mismatchThreshold) continue;
            if (bestMatches < overlapLength - mismatch) {
                bestOverlap = overlapLength;
                bestMatches = overlapLength - mismatch;
            }
        }
        return bestOverlap;
    }

    std::pair<bool,Read> ReadMerger::merge(const Read &r) {
        int overlap1 = calcBestOverlap(r.sequence[0], r.sequence[1], maxMismatch_, 0.1);
        int overlap2 = calcBestOverlap(r.sequence[1], r.sequence[0], maxMismatch_, 0.1);
        if (overlap1 == -1 && overlap2 == -1) {
            return {false, r};
        }
        if (overlap1 < minOverlapLength_ && overlap2 < minOverlapLength_) {
            // filter too short overlaps
            // too short overlap
            // todo: make this a parameter
            return {false, r};
        }
        Read mergedRead;
        mergedRead.name = r.name;
        if (overlap1 > overlap2) {
            mergedRead.sequence[0] = r.sequence[0].substr(0, r.sequence[0].length() - overlap1) + r.sequence[1];
            mergedRead.quality = r.quality.substr(0, r.quality.length() - overlap1) + ' ' + r.quality;
        } else {
            mergedRead.sequence[0] = r.sequence[1].substr(0, r.sequence[1].length() - overlap2) + r.sequence[0];
            mergedRead.quality = r.quality.substr(0, r.quality.length() - overlap2) + ' ' + r.quality;
        }

        mergedRead.length = mergedRead.sequence[0].length();
        mergedRead.sequence[1].resize(mergedRead.length);
        for (int i = 0; i < mergedRead.length; ++i) {
            char c = mergedRead.sequence[0][mergedRead.length - 1 - i];
            switch (c) {
                case 'A':
                    mergedRead.sequence[1][i] = 'T';
                    break;
                case 'T':
                    mergedRead.sequence[1][i] = 'A';
                    break;
                case 'C':
                    mergedRead.sequence[1][i] = 'G';
                    break;
                case 'G':
                    mergedRead.sequence[1][i] = 'C';
                    break;
                default:
                    mergedRead.sequence[1][i] = 'N';
                    break;
            }
        }
        return {true, mergedRead};

    }
}
