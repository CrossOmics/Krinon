// #include <time.h>
#include <cassert>
#include <cstring>
#include <filesystem>
#include "ReadAligner.h"

namespace RefactorProcessing{
    void ReadAlignerSingleThread::setParam(const Parameters &P) {
        isPairedEnd_ = P.isPaired;
        readBufferSize_ = P.readBufferSize;
        outputBufferSize_ = P.outputBufferSize;
        r1Length_ = 0;
        r2Length_ = 0;
        std::cout << "\t[isPaired]: " << P.isPaired << "\n" <<
            "\t[isPaired]: " << P.isPaired << "\n" <<
            "\t[readBufferSize]: " << P.readBufferSize << "\n" <<
            "\t[outputBufferSize]: " << P.outputBufferSize << std::endl;
    }

    void ReadAlignerSingleThread::init(ReadScanner *rs, OutputSAM *o,TimeReport* t,const Parameters &P,int threadId,int threadNum) {
        setParam(P);
        readScanner_ = rs;
        outputSAM_ = o;
        timeReport_ = t;
        threadId_ = threadId;
        threadNum_ = threadNum;

        readBuffer1_ = new char[readBufferSize_ + 20000];
        read1pos_ = 0;
        if (isPairedEnd_){
            readBuffer2_ = new char[readBufferSize_ + 20000];
            read2pos_ = 0;
        }
        outputBuffer_ = new char[outputBufferSize_ + 100000];
        outputPos_ = 0;
        r1Length_ = 0;
        r2Length_ = 0;

        seedMapping_.setParam(P);
        stitchingManagement_.setParam(P);
        std::cout << "Calling init on management" << std::endl;
        stitchingManagement_.init();

        readScanner_ = rs;
        outputSAM_ = o;

        seedMapping_.aligns_.reserve(P.maxSeedPerRead);
    }

    int ReadAlignerSingleThread::getRead() {
        size_t r1Length = 0, r2Length = 0;
        if (read1pos_ == r1Length_) {
            readScanner_->loadFromFastq(readBuffer1_, readBuffer2_, readBufferSize_, r1Length, r2Length);
            read1pos_ = 0;
            r1Length_ = r1Length;
            if (isPairedEnd_) {
                read2pos_ = 0;
                r2Length_ = r2Length;
            }
            if (r1Length == 0) {
                return -1; // no more reads
            }
        }

        readScanner_->parseRead(r, readBuffer1_ + read1pos_, isPairedEnd_ ? readBuffer2_ + read2pos_ : nullptr, r1Length, r2Length);
        read1pos_ += r1Length;
        if (isPairedEnd_) read2pos_ += r2Length;
        return 0;
    }

    void ReadAlignerSingleThread::process() {
        int readCount = 0;
        if (threadId_ == 0) {
            timeReport_->init(threadNum_,"TimeReport.out");
        }
        while (getRead() == 0) {

            seedMapping_.aligns_.clear();
            seedMapping_.process(&r);
            stitchingManagement_.process( &r);
            // std::cout << "Read: " << readCount << " | Transcript Length: " << stitchingManagement_.resultTranscriptLength_ << std::endl;
            // output SAM records
            std::memcpy(outputBuffer_ + outputPos_, stitchingManagement_.resultTranscriptBuffer_, stitchingManagement_.resultTranscriptLength_);
            outputPos_ += stitchingManagement_.resultTranscriptLength_;
            if (outputPos_ > outputBufferSize_) {
                outputSAM_->outputSAM(outputBuffer_, outputPos_);
                outputPos_ = 0;
            }
            ++readCount;
            if (readCount % 1000 == 0) {
                if (threadId_ == 0) {
                    timeReport_->tryActivateReport(threadNum_);
                }
                timeReport_->tryReportProgress(threadId_,readCount);
            }
            // usleep(1000000);
        }
        if (outputPos_ > 0) {
            outputSAM_->outputSAM(outputBuffer_, outputPos_);
            outputPos_ = 0;
        }
    }

    void ReadAligner::setParam(const Parameters &P) {
        threads = P.workingThreads;
        readScanner_.setParam(P);
        /**
         * TODO: Please correct me if I am wrong, but when I discussed with Chandra, we were
         *       under the impression that the output `SAM` file is roughly of the same size
         *       as the input read (i.e. at most a few times bigger).
         *       Knowing this upper bound can help make a completely lock free multi-threaded
         *       program which can achieve much higher throughput with many threads, and only
         *       needs to truncate/resize the output at most once.
         *       But I see that the initial file size here seems to be hardcoded and then expanded
         *       later when needed. Do we actually have to be this careful?
         *       We only need a rough upper bound for the above to work.
         */
        outputSAM_.setParam(P, 1e9);
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
        std::cout << "Current genome length: " << this->gIndex_.genome_.genomeLength_ << std::endl;
        std::cout << "Genome loading done" << std::endl;
        aligners_.reserve(threadNum);
        for (int i = 0; i < threadNum; ++i) {
            aligners_.emplace_back(new ReadAlignerSingleThread(gIndex_));
            std::cout << "Init on aligner thread " << i << std::endl;
            aligners_[i]->init(&readScanner_, &outputSAM_, &timeReport_, P, i, threadNum);
        }
    }

    void ReadAligner::alignReads() {

        // add threads
        std::vector<std::thread> t;
        for (int i = 0; i < threads; ++i) {
            t.emplace_back([this, i] { aligners_[i]->process(); });
        }

        for (auto &thread: t) {
            if (thread.joinable()) {
                thread.join();
            }
        }

        outputSAM_.close();
    }
}
