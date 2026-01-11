#ifndef RNAALIGNREFACTORED_OUTPUTSAM_H
#define RNAALIGNREFACTORED_OUTPUTSAM_H
#include <string>
#include "Parameters.h"

namespace RefactorProcessing {
    class OutputSAM{
        // output alignment results in SAM format
    private:
        // configs
        std::string outputFileName_;
        bool sortByCoordinate_; // for now, keep it false, will implement later
    public:
        OutputSAM(){};
        ~OutputSAM(){};

        void setParam(const Parameters &P);

        void outputOneSAM(std::string &samRecord);
    };
}
#endif //RNAALIGNREFACTORED_OUTPUTSAM_H
