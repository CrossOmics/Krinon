#ifndef RNAALIGNREFACTORED_READALIGNER_H
#define RNAALIGNREFACTORED_READALIGNER_H

#include "../genome/GenomeIndex.h"
#include "../utils/types.h"
#include "../utils/atomic_generation.hpp"
#include <vector>
#include "SeedMapping.h"
#include "StitchingManagement.h"
#include <memory>
#include <fstream>
#include <sstream>
#include <queue>
#include <atomic>
#include "ReadFile.h"
#include "../utils/Parameters.h"
#include "../utils/memory_mapped_output.hpp"

#define RNAALIGNER_PROGRESS_REPORT_INTERVAL_SECONDS 10

namespace rna {
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

        int m_threadIndex;
        int64_t m_periodicReadCount{0};

        /**
         * UPDATE: Made these two members of the classs, since it helps for
         *         other things as well.
         */
        std::chrono::_V2::system_clock::time_point m_alignmentStartTime;
        std::chrono::_V2::system_clock::time_point m_previousProgressReportTime;

        int inputBufferSize = 30000;
        int nowReadInd;
        int nowQueueSize{0};
        bool queueEmpty{true};
        std::vector<ReadPtr> readBufferQueue;

    public:
        bool partialOutput{false};

        ReadAligner(
            const Parameters& P, 
            const GenomeIndex& gInPre, 
            int64_t beginFrom, 
            int64_t expectedReadableBytes, 
            int threadIndex = 0
        );

        void setConfig(
            const StitchingConfig& scfg, 
            const StitchingScoreConfig& sscfg, 
            const SeedMappingConfig& smcfg) 
        {
            seedMapping->setConfig(smcfg);
            stitchingManagement->setConfig(scfg,sscfg);
        }

        bool loadReadFromFastq();
        void reportProgress(
            std::atomic_int64_t& totalReadsProcessed,
            AtomicGenerationCounter& generationCounter
        );
        void processReadFile(
            MemoryMappedFile& outputFileMapped,
            std::atomic_int64_t& outputOffset,
            std::atomic_int64_t& totalReadsProcessed,  
            AtomicGenerationCounter& generationCounter
        );
    };
} // namespace rna
#endif //RNAALIGNREFACTORED_READALIGNER_H
