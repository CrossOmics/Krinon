#ifndef RNAALIGNREFACTORED_READALIGNMULTITHREAD_H
#define RNAALIGNREFACTORED_READALIGNMULTITHREAD_H
#include <thread>
#include <mutex>
#include "ReadAligner.h"
#include "../utils/exceptions.h"
#include "ReadFile.h"

namespace rna{
    class Parameters;
    class ReadAlignMultiThread {
        ReadFile readFile;
        std::mutex outputLock; // mutex for output operations
        std::mutex alignStatusLock; // mutex for alignment status output
        std::mutex alignProgressLock; // mutex for alignment progress output
        int threadNum;
        std::vector<std::thread> threads;
        int readCount{0};
        std::vector<ReadAligner> readAligners;
        //char* outputAlignBuffer; //buffer for output alignment
        //FILE* outFile;
        int outFile;
        std::string outDir;
        std::ofstream logFile;
        std::ofstream alignProgressFile;

        StitchingConfig stitchConfig;
        StitchingScoreConfig stitchingScoreConfig;
        SeedMappingConfig seedMappingConfig;

        void singleThreadProcess(ReadAligner&r) {

            r.processReadFile(readFile, outFile, logFile, alignProgressFile, outputLock, alignStatusLock, alignProgressLock, readCount);
        }
    public:
        ReadAlignMultiThread(Parameters& P);
        ~ReadAlignMultiThread(){
            //delete[] outputAlignBuffer;
        };
        void setConfig(const StitchingConfig& scfg, const StitchingScoreConfig& sscfg, const SeedMappingConfig& smcfg) {
            stitchConfig = scfg;
            stitchingScoreConfig = sscfg;
            seedMappingConfig = smcfg;
        }
        void processReadFile(int tNum, GenomeIndex& gInPre, bool partialOutput );
    };
}

#endif //RNAALIGNREFACTORED_READALIGNMULTITHREAD_H
