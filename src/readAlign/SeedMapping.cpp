#include "SeedMapping.h"

namespace rna {

    // split into high quality splits
    void SeedMapping::splitRead() {
        splits.clear();
        size_t splitLen = 0;
        size_t maxLen = 0;
        for (int dir = 0; dir < 2; ++dir) {
            size_t len = readSeq[dir].length();
            ReadPos i;
            for (i = 0; i < len; ++i) {
                if (readSeq[dir][i] != 'N' ) {
                    splitLen++;
                } else {
                    if (splitLen >= config.minSplitLength) {
                        splits.emplace_back(Split{ReadPos(i - splitLen), splitLen, dir});
                    }
                    maxLen = std::max(maxLen, splitLen);
                    splitLen = 0;
                }
            }
            //last split
            if (splitLen >= config.minSplitLength) {
                splits.emplace_back(Split{ReadPos(i - splitLen), splitLen, dir});
            }
            splitLen = 0;

        }
    }

    void SeedMapping::alignRead() {
        aligns.reserve(readSeq->length()/5);
        for (Split &split: splits) {


            std::vector<Align> splitAligns=genomeIndexPrefix.alignRead(readSeq[split.direction], split);
             for (Align &align: splitAligns) {
                if(align.rep  != 0) //which means the seed is invalid (i.e. repeat too many times)
                    aligns.emplace_back(align);
            }
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