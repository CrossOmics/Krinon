#ifndef RNAALIGNREFACTORED_SEEDMAPPING_H
#define RNAALIGNREFACTORED_SEEDMAPPING_H
#include "Read.h"
#include "../genomeRefactored/GenomeIndex.h"
#include "../io/Parameters.h"
#include <vector>
namespace RefactorProcessing {

    class SeedMapping {
    private:
        //configs
        int minSplitLength_;
        int maxSeedPerRead_;

        //data
        std::vector<Split> splits_;
        GenomeIndex genomeIndex_;
        Read* read_;


    public:
        std::vector<Align> aligns_;

        void setParam(const Parameters& P);

        void alignRead();

        void splitRead();

        void clear();

        void process(Read* r);

    };
}

#endif //RNAALIGNREFACTORED_SEEDMAPPING_H
