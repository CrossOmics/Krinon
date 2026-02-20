#ifndef RNAALIGNREFACTORED_OUTPUTSAM_H
#define RNAALIGNREFACTORED_OUTPUTSAM_H
#include <string>
#include <mutex>
#include "Parameters.h"
#include "MemoryMappedFile.hpp"


namespace RefactorProcessing {
    class OutputSAM{
        // output alignment results in SAM format
    private:
        // configs
        std::string outputFileName_;
        bool sortByCoordinate_; // for now, keep it false, will implement later

        //multi-threading
        std::mutex outputFileLock_; // lock for writing to output file, may optimize later

        //data
        std::unique_ptr<MemoryMappedFile> outputFile_;
        size_t mmapPos_; //offset for memory-mapped output file
        size_t fileSize_; //size of memory-mapped output file



    public:
        OutputSAM(){};
        ~OutputSAM(){};

        void setParam(const Parameters &P,int initialFileSize = 1e7);

        void outputSAM(char* samRecord,size_t size);
    };
}
#endif //RNAALIGNREFACTORED_OUTPUTSAM_H
