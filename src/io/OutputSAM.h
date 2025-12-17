#ifndef RNAALIGNREFACTORED_OUTPUTSAM_H
#define RNAALIGNREFACTORED_OUTPUTSAM_H
#include <string>
namespace RefactorProcessing {
    class OutputSAM{
        // output alignment results in SAM format
    private:
        std::string outputFileName_;
        bool sortByCoordinate_; // for now, keep it false, will implement later
    public:
        OutputSAM(){};
        ~OutputSAM(){};
    };
}
#endif //RNAALIGNREFACTORED_OUTPUTSAM_H
