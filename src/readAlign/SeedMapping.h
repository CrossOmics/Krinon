#ifndef RNAALIGNREFACTORED_SEEDMAPPING_H
#define RNAALIGNREFACTORED_SEEDMAPPING_H
#include "../genome/GenomeIndex.h"
#include "../utils/types.h"
#include "../genome/GenomeIndexPrefix.h"
#include <vector>
#include <algorithm>

namespace rna {
    struct SeedMappingConfig {
        size_t minSplitLength = 20;

    };
    class SeedMapping {
    private:

        GenomeIndexPrefix& genomeIndexPrefix;
        SeedMappingConfig config;
        ReadPtr read;
        std::string readSeq[2];
        std::vector<Split> splits;
    public:
        std::vector<Align> aligns;
        SeedMapping(GenomeIndexPrefix& gInPre, const SeedMappingConfig& cfg) : genomeIndexPrefix(gInPre), config(cfg) {}
        void alignRead();
        void splitRead();
        void clear() {
            aligns.clear();
            splits.clear();
        }
        void processRead(ReadPtr& r);


    };
}
#endif //RNAALIGNREFACTORED_SEEDMAPPING_H
