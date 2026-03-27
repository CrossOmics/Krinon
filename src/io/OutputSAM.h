#ifndef RNAALIGNREFACTORED_OUTPUTSAM_H
#define RNAALIGNREFACTORED_OUTPUTSAM_H
#include <atomic>
#include <string>
#include <filesystem>
#include "Parameters.h"
#include "MemoryMappedInput.hpp"
#include "MPSC.hpp"


namespace RefactorProcessing {
    using SAMOutputQueue = MPSCQueue<std::string>;
    using BufferedReads = QueueItem<std::string>;
    using OutputBufferPool = BufferPool<std::string>;

    class OutputSAM {
    private:
        // Full path to the output file
        std::string outputFilePath_;
        // for now, keep it false, will implement later
        bool sortByCoordinate_;
        // Output file descriptor
        int outputFD_;
        // The queue used to gather buffered reads from aligners
        SAMOutputQueue queue_;
        // Pool of pre-allocated buffers, passed to the consumer thread when full to write
        OutputBufferPool bufferPool_;

    public:
        OutputSAM(const Parameters& P)
            // Each thread gets one rather large buffer to write to
            : bufferPool_(P.workingThreads)
        {
            outputFilePath_ = (std::filesystem::path(P.outPutDir) / "outAligned.out.sam").string();
            std::cout << "Output file path: " << outputFilePath_ << std::endl;
            outputFD_ = open(outputFilePath_.c_str(), O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
            if (outputFD_ == -1) {
                throw std::runtime_error("Failed get file handle for " + outputFilePath_);
            }
        };

        ~OutputSAM() {
            if (outputFD_ != -1) {
                close(outputFD_);
            }
        };

        void submitBufferedReads(BufferedReads* reads) {
            queue_.push(reads);
        }

        /**
         * Attempt to get an output buffer. This call may block the caller
         * in case the writer thread is slow and we don't have any buffers 
         * to write to.
         */
        BufferedReads* getOutputBuffer() {
            BufferedReads* reads = bufferPool_.pop();
            while (!reads) {
                std::this_thread::yield();
                reads = bufferPool_.pop();
            }
            return reads;
        }

        void relinquishBuffer(BufferedReads* reads) {
            bufferPool_.push(reads);
        }

        /**
         * The (single) thread that consumes buffered reads and writes
         * directly to the output SAM file.
         * Receives a reference to an atomic counter that shows how many
         * threads are still processing things. Once this hits zero, this
         * thread can exit safely.
         */
        void consumerThread(std::atomic<int>& runningProducers);
    };
}
#endif //RNAALIGNREFACTORED_OUTPUTSAM_H
