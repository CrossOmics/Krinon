#ifndef RNAALIGNREFACTORED_OUTPUTSAM_H
#define RNAALIGNREFACTORED_OUTPUTSAM_H
#include <string>
#include "Parameters.h"
#include "MemoryMappedFile.hpp"

namespace RefactorProcessing {
    class OutputSAM{
        // output alignment results in SAM format
    private:
        // configs
        std::string outputFileName_;
        bool sortByCoordinate_; // for now, keep it false, will implement later

        //data

    public:
        OutputSAM(){};
        ~OutputSAM(){};

        void setParam(const Parameters &P);

        void outputSingleSAM(std::string &samRecord);
    };
}
#endif //RNAALIGNREFACTORED_OUTPUTSAM_H
