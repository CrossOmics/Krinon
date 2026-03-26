#ifndef RNAALIGNREFACTORED_READALIGNER_H
#define RNAALIGNREFACTORED_READALIGNER_H

//todo move out of readAlignRefactored later
#include "../io/ReadScanner.h"
#include "../genomeRefactored/GenomeIndex.h"
#include "SeedMapping.h"
#include "Stitching.h"
#include "../readAlignRefactored/Transcript.h"
#include "../io/OutputSAM.h"
#include "../log/TimeReport.h"

namespace RefactorProcessing {
    class ReadAlignerSingleThread {
    private:
        // configs
        bool isPairedEnd_;
        size_t readBufferSize_;
        size_t outputBufferSize_;

        // multi-thread
        int threadId_;
        int threadNum_;

        // current read
        Read r;
        
        // all readAligners share a read scanner
        ReadScanner *readScanner_;
        SeedMapping seedMapping_;
        Stitching stitchingManagement_;
        
        // Pointer to the output SAM wrapper (without ownership)
        const std::shared_ptr<OutputSAM>& outputSAM_;
        // Pointer to the time report object (without ownership)
        const std::shared_ptr<TimeReport>& timeReport_;
        // todo add time and performance report
    public:
        ReadAlignerSingleThread(const GenomeIndex& g, const std::shared_ptr<OutputSAM>& output, const std::shared_ptr<TimeReport>& report) : 
            seedMapping_(g), 
            stitchingManagement_(g,seedMapping_.aligns_),
            outputSAM_(output),
            timeReport_(report)
        {};
        // initialize data structures, called once
        void init(ReadScanner *rs, const Parameters &P, int threadId, int threadNum);

        void setParam(const Parameters &P);

        /**
         * Record the next read in the interval given in `block` and on offset
         * `offset`. The beginning of `block` _MUST_ be aligned at a valid
         * read header.
         * Return the number of bytes read.
         */
        size_t getRead(std::string_view block, size_t offset);

        void process(std::atomic<int>& activeThreads);
        // We use 8MB output buffers
        const size_t MAX_OUTPUT_BUFFER_SIZE = 8ULL * 1024 * 1024;
    };

    class ReadAligner {
    private:
        int threads;
        std::string version{"v0"};
        std::string gIndexDir_;
        GenomeIndex gIndex_;
        ReadScanner readScanner_;
        std::vector<std::unique_ptr<ReadAlignerSingleThread>> aligners_;
        std::shared_ptr<OutputSAM> outputSAM_{nullptr};
        std::shared_ptr<TimeReport> timeReport_{nullptr};
        std::atomic<int> activeThreads_{0};
    
    public:
        void setParam(const Parameters &P);
        void loadGenome();
        void init(const Parameters &P, int threadNum);
        void reset();
        void alignReads();
    };
}

#endif //RNAALIGNREFACTORED_READALIGNER_H
