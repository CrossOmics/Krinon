#include "SeedMapping.h"

namespace rna {

    // split into high quality splits
    void SeedMapping::splitRead() {
        splits.clear();

        size_t splitLen = 0;
        size_t maxLen = 0;

        size_t len = readSeq[0].length();
        ReadPos i;
        int iFragment = 0; // which fragment in paired-end
        for (i = 0; i < len; ++i) {
            if (readSeq[0][i] != 'N' && readSeq [0][i] != '#') {
                splitLen++;
            } else {

                if (splitLen >= config.minSplitLength) {
                    splits.emplace_back(Split{ReadPos(i - splitLen), splitLen, iFragment});
                }
                if (readSeq [0][i] == '#') iFragment++;
                maxLen = std::max(maxLen, splitLen);
                splitLen = 0;
            }
        }
        //last split
        if (splitLen >= config.minSplitLength) {
            splits.emplace_back(Split{ReadPos(i - splitLen), splitLen, iFragment});
        }


    }

    void SeedMapping::alignRead() {
        for (Split &split: splits) {

            genomeIndex.alignRead(readSeq, split, aligns, alignNum, config.maxSeedPerRead);

        }
    }

    void SeedMapping::processRead(ReadPtr &r) {
        clear();
        read = r;
        readSeq[0] = read->sequence[0];
        readSeq[1] = read->sequence[1]; // complement sequence

        splitRead();
        alignRead();
    }
}