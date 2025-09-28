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
#include "../readAlign/StitchingManagement.h"
#include <yaml-cpp/yaml.h>
#include <linux/limits.h>

namespace rna{
    class Parameters {
    public:
        argparse::ArgumentParser program;

        /**
         * A global config node that contains the default values for all parameters
         * Inspect `config.yaml` to see what they are.
         */
        YAML::Node m_globalConfig;

        /**
         * This helps a lot to keep track of
         */
        std::filesystem::path m_binaryDir;

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
