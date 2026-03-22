// #include <time.h>
#include <cassert>
#include <cstring>
#include <filesystem>
#include "ReadAligner.h"

namespace RefactorProcessing{
    void ReadAlignerSingleThread::setParam(const Parameters &P) {
        isPairedEnd_ = P.isPaired;
    }

    void ReadAlignerSingleThread::init(ReadScanner *rs, const Parameters &P, int threadId, int threadNum) {
        setParam(P);
        readScanner_ = rs;
        threadId_ = threadId;
        threadNum_ = threadNum;
        seedMapping_.setParam(P);
        stitchingManagement_.setParam(P);
        stitchingManagement_.init();
        readScanner_ = rs;
        seedMapping_.aligns_.reserve(P.maxSeedPerRead);
    }

    size_t ReadAlignerSingleThread::getRead(std::string_view block, size_t offset) {
        // The offset may be outside of the current buffer, then we are done with it!
        if (offset >= block.length()) {
            return 0;
        }
        // Otherwise, get the pointer to the current and parse it
        return readScanner_->parseRead(r, block.data() + offset, nullptr);
    }

    void ReadAlignerSingleThread::process(std::atomic<int>& activeThreads) {
        int readCount = 0;
        if (threadId_ == 0) {
            timeReport_->init(threadNum_, "TimeReport.out");
        }

        // Get our first output buffer instance
        BufferedReads* bufferedReads = outputSAM_->getOutputBuffer();
        while (true) {
            // Claim a read buffer properly aligned to start from a valid read
            std::string_view readBuffer = readScanner_->loadFromFastq();
            // There may be nothing to read anymore (i.e. the current read partition
            // is finished, there is nothing left in the mapped region). In this case,
            // we are done and can break out.
            if (readBuffer.length() == 0) {
                break;
            }
            // Read until the buffer is exhausted ...
            size_t offset = 0, bytes = 0;
            while (true) {
                bytes = getRead(readBuffer, offset);
                offset += bytes;
                if (bytes == 0) {
                    // The current read buffer has been exhausted, break out and
                    // get a new one!
                    break;
                }
                // Map and write directly to the buffer (we _own_ this, and since it is
                // a string, we don't have to manually resize it).
                seedMapping_.aligns_.clear();
                seedMapping_.process(&r);
                stitchingManagement_.process(&r, bufferedReads->data);
                // If the output buffer is filled, pass the result to the writer
                if (bufferedReads->data.length() >= MAX_OUTPUT_BUFFER_SIZE) {
                    // Submit all buffered reads, then attempt to grab a new buffer to write again
                    outputSAM_->submitBufferedReads(bufferedReads);
                    bufferedReads = outputSAM_->getOutputBuffer();
                }
                ++readCount;
                if (readCount % 1000 == 0) {
                    if (threadId_ == 0) {
                        timeReport_->tryActivateReport(threadNum_);
                    }
                    timeReport_->tryReportProgress(threadId_,readCount);
                }
            }
        }
        // In case there are leftovers to write, just submit!
        if (!bufferedReads->data.empty()) {
            outputSAM_->submitBufferedReads(bufferedReads);
        }
        // Announce that you are done with this partition.
        activeThreads.fetch_sub(1, std::memory_order_release);
    }

    void ReadAligner::setParam(const Parameters &P) {
        threads = P.workingThreads;
        readScanner_.setParam(P);
        gIndex_.setParam(P);
        gIndexDir_ = (std::filesystem::path(P.genomeGenerateFileStoreDir) / "GeIndex").string();
    }

    void ReadAligner::loadGenome() {
        /**
         * TODO: `loadFromFile` handles many error cases. I cannot tell which one of them
         *       is fatal and which one can be handled later, so for now, I am forcing it
         *       to be error free.
         *       Since this is before calling `mmap`, we can just abort without leaking.
         */
        gIndex_.loadFromFile(gIndexDir_);
    }

    void ReadAligner::init(const Parameters &P, int threadNum) {
        setParam(P);
        std::cout << "Loading genome file(s) at " << gIndexDir_ << std::endl;
        loadGenome();
        std::cout << "Mapping read file ... " << std::endl;
        outputSAM_ = std::make_shared<OutputSAM>(P);
        timeReport_ = std::make_shared<TimeReport>();
        aligners_.reserve(threadNum);
        for (int i = 0; i < threadNum; ++i) {
            aligners_.emplace_back(
                new ReadAlignerSingleThread(
                    gIndex_, outputSAM_, timeReport_
                )
            );
            std::cout << "Init on aligner thread " << i << std::endl;
            aligners_[i]->init(&readScanner_, P, i, threadNum);
        }
    }

    void ReadAligner::reset() {
        activeThreads_.store(threads, std::memory_order_release);
    }

    void ReadAligner::alignReads() {
        /**
         * TODO: (ARVIN) This is only for one partition. Add logic to move on
         *       to the next partition, if there is any!
         */
        // Reset number of active threads
        reset();
        // Add aligner threads
        std::vector<std::thread> t;
        for (int i = 0; i < threads; ++i) {
            t.emplace_back([this, i] { aligners_[i]->process(activeThreads_);});
        }
        // Assume the role of the consumer thread
        outputSAM_->consumerThread(activeThreads_);
        // Join for cleanup ...
        for (auto &thread: t) {
            if (thread.joinable()) {
                thread.join();
            }
        }
    }
}
