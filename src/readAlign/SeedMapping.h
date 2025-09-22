#ifndef RNAALIGNREFACTORED_SEEDMAPPING_H
#define RNAALIGNREFACTORED_SEEDMAPPING_H

#include "../utils/types.h"
#include "../genome/GenomeIndex.h"
#include <vector>
#include <algorithm>

namespace rna {
    struct SeedMappingConfig {
        int minSplitLength = 20;
        int maxSeedPerRead = 1000;
    };
    class SeedMapping {
    private:

        const GenomeIndex& genomeIndex;
        SeedMappingConfig config;
        ReadPtr read;
        std::string readSeq[2];
        std::vector<Split> splits;
    public:
        Align* aligns;
        int alignNum;
        SeedMapping(const GenomeIndex& gInPre, const SeedMappingConfig& cfg) : genomeIndex(gInPre), config(cfg) {
            aligns = new Align[config.maxSeedPerRead];
            alignNum = 0;
        }
        ~SeedMapping(){
            delete []aligns;
        }
        void alignRead();
        void splitRead();
        void clear() {
            alignNum = 0;
            splits.clear();
        }
        void processRead(ReadPtr& r);
        void setConfig(const SeedMappingConfig& cfg) {
            config = cfg;
        }


    };
}
#endif //RNAALIGNREFACTORED_SEEDMAPPING_H
