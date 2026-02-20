#include <cstring>
#include "ReadAligner.h"
namespace RefactorProcessing{
    void ReadAlignerSingleThread::setParam(const Parameters &P) {
        isPairedEnd_ = P.isPaired;
        readBufferSize_ = 5e7; // todo replace by parameter
        outputBufferSize_ = 5e7; // todo replace by parameter
        r1Length_ = 0;
        r2Length_ = 0;

    }

    void ReadAlignerSingleThread::init(ReadScanner *rs, OutputSAM *o,const Parameters &P,int threadId) {
        setParam(P);
        readScanner_ = rs;
        outputSAM_ = o;
        threadId_ = threadId;

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
        stitchingManagement_.init(&aligns_);

        readScanner_ = rs;
        outputSAM_ = o;

        aligns_.resize(P.maxSeedPerRead);
    }

    int ReadAlignerSingleThread::getRead() {
        size_t r1Length, r2Length;
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
        while (getRead() == 0) {
            aligns_.clear();
            seedMapping_.process(&r);
            stitchingManagement_.process(aligns_, &r);
            // output SAM records
            std::memcpy(outputBuffer_, stitchingManagement_.resultTranscriptBuffer_,stitchingManagement_.resultTranscriptLength_);
            outputPos_ += stitchingManagement_.resultTranscriptLength_;
            if (outputPos_ > outputBufferSize_) {
                outputSAM_->outputSAM(outputBuffer_, outputPos_);
                outputPos_ = 0;
            }
            ++readCount;
        }
    }

    void ReadAligner::setParam(const Parameters &P) {
        threads = P.workingThreads;
        readScanner_.setParam(P);
        outputSAM_.setParam(P, 1e9);
        gIndex_.setParam(P);
        gIndexDir_ = P.genomeGenerateFileStoreDir + "GeIndex";
    }

    void ReadAligner::loadGenome() {
        gIndex_.loadFromFile(gIndexDir_);
    }

    void ReadAligner::init(const Parameters &P, int threadNum) {
        setParam(P);
        aligners_.reserve(threadNum);
        for (int i = 0; i < threadNum; ++i) {
            aligners_.emplace_back(new ReadAlignerSingleThread(gIndex_));
            aligners_[i]->init(&readScanner_, &outputSAM_, P, i);
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
    }


}
