#ifndef RNAALIGNREFACTORED_READALIGNMULTITHREAD_H
#define RNAALIGNREFACTORED_READALIGNMULTITHREAD_H
#include <thread>
#include <mutex>
#include "ReadFile.h"
#include "ReadAligner.h"
#include "../utils/exceptions.h"
#include "../utils/atomic_generation.hpp"
#include "../utils/memory_mapped_output.hpp"
#include <plog/Log.h>

namespace rna{
    class Parameters;
    class ReadAligner;
    class ReadAlignMultiThread {
        // ReadFile readFile;
        // std::mutex outputLock;                  // mutex for output operations
        // std::mutex alignStatusLock;             // mutex for alignment status output
        int threadNum;

        /**
         * Atomic counter for how many things we have read across all
         * readers combined.
         * Incremented periodically.
         */
        std::atomic_int64_t m_totalReadCount{0};

        /**
         * The current offset into the memory-mapped output file.
         */
        std::atomic_int64_t m_outputOffset{0};
        
        /**
         * An atomic latch that keeps track of who has updated its
         * read count recently.
         * When a reader sees that it has reset (i.e. becamse 0 then went
         * back to maximum value), then it print the periodic report.
         */
        AtomicGenerationCounter m_generationCounter{0};

        std::string m_inputFile1;
        std::string m_inputFile2;

        /**
         * Vector of `ReadAligner` objects spawned for each thread.
         */
        std::vector<ReadAligner> m_readAligners;
        /**
         * Vector of threads running the alginment procedure for each
         * object in `m_readAligners`.
         */
        std::vector<std::thread> m_threads;
        
        // char* outputAlignBuffer;                //buffer for output alignment
        // FILE* outFile;
        // std::string outDir;
        // std::ofstream logFile;
        // std::ofstream alignProgressFile;

        StitchingConfig stitchConfig;
        StitchingScoreConfig stitchingScoreConfig;
        SeedMappingConfig seedMappingConfig;

    public:
        ReadAlignMultiThread(Parameters& P);
        void setConfig(const StitchingConfig& scfg, const StitchingScoreConfig& sscfg, const SeedMappingConfig& smcfg) {
            stitchConfig = scfg;
            stitchingScoreConfig = sscfg;
            seedMappingConfig = smcfg;
        }
        /**
         * Create `ReadAligner` objects for each thread and set the starting
         * offset to begin reading from on the input file.
         */
        void initiateAligners(
            rna::Parameters& params, 
            rna::GenomeIndex &gInPre,
            int64_t totalBytesToRead, 
            bool partialOutput
        );
        /**
         * Spawn aligner threads and join
         */
        void doMultiThreadedAlignment(
            rna::Parameters& params,
            MemoryMappedFile& output
        );
        void processReadFile(
            rna::Parameters& params, 
            GenomeIndex& gInPre, 
            bool partialOutput
        );
        void singleThreadProcess(
            ReadAligner&r, 
            int threadIndex, 
            MemoryMappedFile& output
        );
    };
}

#endif //RNAALIGNREFACTORED_READALIGNMULTITHREAD_H
