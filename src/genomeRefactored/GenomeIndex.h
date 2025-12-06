#ifndef RNAALIGNREFACTORED_GENOMEINDEX_H
#define RNAALIGNREFACTORED_GENOMEINDEX_H
#include "Genome.h"
#include "../io/Parameters.h"
namespace RefactorProcessing{
    class GenomeIndex{
    private:
        // configs

        int kMerSize_;
        // data

    public:
        GenomeIndex();
        ~GenomeIndex();

        // parameters
        void setParam(const Parameters& P);
    };
}

#endif //RNAALIGNREFACTORED_GENOMEINDEX_H
