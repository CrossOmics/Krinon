#ifndef RNAALIGNREFACTORED_READALIGNER_H
#define RNAALIGNREFACTORED_READALIGNER_H


#include "../genome/GenomeIndex.h"
#include "../utils/types.h"
#include <vector>
#include "SeedMapping.h"
#include "StitchingManagement.h"
#include <memory>
#include <fstream>
#include <sstream>
#include <mutex>
#include <queue>
#include "ReadFile.h"
#include "../utils/Parameters.h"

namespace rna {
    // class ReadFile;
    class Parameters;
    class ReadAligner {
    private:
        const GenomeIndex& genomeIndexPrefix;

        /**
         * UPDATE: Moved this from the multi-threaded implementation into here.
         *         Each aligner instance now has its own stream.
         */
        std::unique_ptr<ReadFile> m_readFile;
        int64_t m_beginFrom{0};
        int64_t m_expectedReadableBytes{0};

        std::unique_ptr<SeedMapping> seedMapping;                   // seedMapping, get alignments of the read
        std::unique_ptr<StitchingManagement> stitchingManagement;   // stitch the previously obtained alignments into transcripts
        ReadPtr read;                                               // the read being processed now
        std::vector<TranscriptPtr> bestTranscripts;                 // the transcripts generated from the read
        std::vector<SJDBOutput> sjdbCandidates;                     // splice junction candidates from the read

        int threadId{0};
        int m_oneMinuteCycleReadCount{0};
        int inputBufferSize = 30000;
        int nowReadInd;
        int nowQueueSize{0};
        bool queueEmpty{true};
        std::vector<ReadPtr> readBufferQueue;

    public:
        bool partialOutput{false};

        void setConfig(const StitchingConfig& scfg, const StitchingScoreConfig& sscfg, const SeedMappingConfig& smcfg) {
            seedMapping->setConfig(smcfg);
            stitchingManagement->setConfig(scfg,sscfg);
        }
        // ReadAligner(const GenomeIndex& gInPre, int tId = 0) : genomeIndexPrefix(gInPre), threadId(tId) {
        ReadAligner(const Parameters& P, const GenomeIndex& gInPre, int64_t beginFrom, int64_t expectedReadableBytes, int tId = 0);

        // bool loadReadFromFastq(ReadFile& file);
        bool loadReadFromFastq();
        // void processReadFile(ReadFile& file,FILE* outFile,std::ofstream& alignStatusFile,std::ofstream& alignProgressFile,std::mutex& outputLock,std::mutex& alignStatusLock,std::mutex& alignProgressLock,int& totalReadsProcessed);
        // void processReadFile(ReadFile& file, int& totalReadsProcessed, const int threadIndex = 0);
        void processReadFile(std::string outFilePath, int& totalReadsProcessed, const int numThreads, const int threadIndex, 
            std::mutex& progressLock, int& generationCounter);
    };
} // namespace rna
#endif //RNAALIGNREFACTORED_READALIGNER_H
