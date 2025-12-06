#ifndef RNAALIGNREFACTORED_PARAMETERS_H
#define RNAALIGNREFACTORED_PARAMETERS_H
#include "argparse.hpp"
#include "../genome/Genome.h"
#include "../genome/GenomeIndex.h"
#include "../genome/SJDB.h"
#include "../genome/SuffixArray.h"
#include "../readAlign/ReadAligner.h"
#include "../readAlign/ReadAlignMultiThread.h"
#include "../readAlign/SeedMapping.h"
#include "../readAlign/Stitching.h"
namespace rna{
    class Parameters {
    public:
        argparse::ArgumentParser program;
        std::string mode;
        int genomeBinSize;
        std::string genomeFile;
        bool isPaired;
        std::string readFile;
        std::string readFile2;
        std::string gtfFile;
        int threads;
        std::string outPutDir;
        std::string genomeGenerateFileStoreDir;
        std::ofstream outLogFile;

        GenomeIndexConfig genomeIndexConfig;
        GTFConfig gtfConfig;
        StitchingConfig stitchConfig;
        StitchingScoreConfig stitchingScoreConfig;
        SeedMappingConfig seedMappingConfig;

        Parameters();
        int process(int argc, char* argv[]);

    };
}

#endif //RNAALIGNREFACTORED_PARAMETERS_H
