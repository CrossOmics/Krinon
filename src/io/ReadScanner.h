#ifndef RNAALIGNREFACTORED_READSCANNER_H
#define RNAALIGNREFACTORED_READSCANNER_H

#include <thread>
#include <mutex>
#include <fstream>
#include <memory>
#include "Parameters.h"
#include "MemoryMappedFile.hpp"
#include "../readAlignRefactored/Read.h"

namespace RefactorProcessing {
    // scan reads, handle multi-threading
    class ReadScanner {
    private:
        std::mutex readFileLock_;

        // configs
        bool isPairedEnd_;
        std::string readFileName1_;
        std::string readFileName2_;

        size_t readBufferSize_;
        //char* readBuffer1_;
        //char* readBuffer2_;


        //data
        //std::ifstream readFile1_;
        //std::ifstream readFile2_;

        // memory-mapped fast path
        std::unique_ptr<MemoryMappedFile> mmapFile1_;
        std::unique_ptr<MemoryMappedFile> mmapFile2_;
        size_t mmapPos1_ = 0;
        size_t mmapPos2_ = 0;
        size_t fileSize1_ = 0;
        size_t fileSize2_ = 0;

    public:
        void setParam(const Parameters &P);

        void loadFromFastq(char* targetBuffer1,char* targetBuffer2,const size_t bufferSize,size_t &buffer1Length,size_t & buffer2Length);

        void parseRead(Read &r,const char* readBuffer1,const char* readBuffer2,size_t& r1Length,size_t& r2Length);
    };
}
#endif //RNAALIGNREFACTORED_READSCANNER_H
