#include <cstring>
#include "OutputSAM.h"
#include "../log/ErrorRecord.h"
namespace RefactorProcessing{
    void OutputSAM::setParam(const Parameters &P, int initialFileSize) {
        outputFileName_ = P.outPutDir + "outAligned.out.sam";
        sortByCoordinate_ = false;
        // memory-map the output file
        try {
            outputFile_.reset(new MemoryMappedFile(outputFileName_, (long long)initialFileSize));
            mmapPos_ = 0;
            fileSize_ = outputFile_->size();
        } catch (const std::exception &e) {
            rna::ErrorRecord().reportError(std::string("Memory map failed: ") + e.what());
        }
    }

    void OutputSAM::outputSAM(char* samRecord, size_t size) {


        if (!outputFile_ || fileSize_ == 0) {
            rna::ErrorRecord().reportError("Output file is not available for writing");
            return;
        }
        outputFileLock_.lock();
        if (mmapPos_ + size > fileSize_) {
            outputFile_->ensureSize(fileSize_ * 2); // double the file size
            fileSize_ = outputFile_->size();
        }
        char* base = outputFile_->getMapPtr() + mmapPos_;
        mmapPos_ += size;
        outputFileLock_.unlock();
        std::memcpy(base, samRecord, size);


    }

    void OutputSAM::close() {
        if (outputFile_) {
            outputFile_->truncate(mmapPos_); // truncate the file to the actual size of the data written
            outputFile_->memClose();
        }
    }
}
