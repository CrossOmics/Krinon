#include "SeedMapping.h"
#include <string>

namespace RefactorProcessing {
    void SeedMapping::setParam(const Parameters &P) {
        minSplitLength_ = P.minSplitLength;
        maxSeedPerRead_ = P.maxSeedPerRead;
        splits_.reserve(maxSeedPerRead_/2);
    }

    void SeedMapping::alignRead() {
        for (Split &split: splits_) {
            genomeIndex_.find(split, aligns_);
        }
    }

    void SeedMapping::splitRead() {
        splits_.clear();

        int len = read_->length;
        const std::string &readSeq = read_->sequence[0];
        const std::string &readSeqComp = read_->sequence[1];


        int splitLen = 0;
        int maxLen = 0;

        int i;
        int iFragment = 0; // which fragment in paired-end
        for (i = 0; i < len; ++i) {
            if (readSeq[i] != 'N' && readSeq[i] != '#') {
                splitLen++;
            } else {

                if (splitLen >= minSplitLength_) {
                    //valid split
                    splits_.emplace_back(
                            Split{
                                    i - splitLen,
                                    splitLen,
                                    iFragment,
                                    readSeq.substr(i - splitLen, splitLen),
                                    readSeqComp.substr(len - i, splitLen),
                            }
                    );
                }
                if (readSeq[i] == '#') iFragment++;
                maxLen = std::max(maxLen, splitLen);
                splitLen = 0;
            }
        }

        // the last split
        if (splitLen >= minSplitLength_) {
            splits_.emplace_back(
                    Split{
                            i - splitLen,
                            splitLen,
                            iFragment,
                            readSeq.substr(i - splitLen, splitLen),
                            readSeqComp.substr(len - i, splitLen),
                    }
            );
        }
    }

    void SeedMapping::clear() {
        aligns_.clear();
        splits_.clear();
    }



    void SeedMapping::process(Read *r) {
        clear();
        read_ = r;
        splitRead();
        alignRead();
    }
}