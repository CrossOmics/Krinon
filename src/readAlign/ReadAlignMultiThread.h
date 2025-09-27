#ifndef RNAALIGNREFACTORED_READALIGNMULTITHREAD_H
#define RNAALIGNREFACTORED_READALIGNMULTITHREAD_H
#include <thread>
#include <mutex>
#include "ReadFile.h"
#include "ReadAligner.h"
#include "../utils/exceptions.h"
#include <plog/Log.h>

namespace rna{
    class Parameters;
    class ReadAligner;
    class ReadAlignMultiThread {
        // ReadFile readFile;
        // std::mutex outputLock;                  // mutex for output operations
        // std::mutex alignStatusLock;             // mutex for alignment status output
        std::mutex m_alignProgressLock;           // mutex for alignment progress output
        int threadNum;
        int readCount{0};
        int m_generationCounter{0};
        std::string m_inputFile1;
        std::string m_inputFile2;
        std::vector<std::thread> threads;
        std::vector<ReadAligner> readAligners;
        // char* outputAlignBuffer;                //buffer for output alignment
        // FILE* outFile;
        // std::string outDir;
        // std::ofstream logFile;
        // std::ofstream alignProgressFile;

        StitchingConfig stitchConfig;
        StitchingScoreConfig stitchingScoreConfig;
        SeedMappingConfig seedMappingConfig;

        void singleThreadProcess(ReadAligner&r, int threadIndex, std::string outputDir);
    public:
        ReadAlignMultiThread(Parameters& P);
        // ~ReadAlignMultiThread(){
        //     delete[] outputAlignBuffer;
        // };
        // ~ReadAlignMultiThread();
        void setConfig(const StitchingConfig& scfg, const StitchingScoreConfig& sscfg, const SeedMappingConfig& smcfg) {
            stitchConfig = scfg;
            stitchingScoreConfig = sscfg;
            seedMappingConfig = smcfg;
        }
        void processReadFile(rna::Parameters& params, GenomeIndex& gInPre, bool partialOutput );
    };
}

#endif //RNAALIGNREFACTORED_READALIGNMULTITHREAD_H
