#ifndef RNAALIGNREFACTORED_ERRORRECORD_H
#define RNAALIGNREFACTORED_ERRORRECORD_H

#include <string>
#include <fstream>
#include "../io/Parameters.h"
namespace rna {
    class ErrorRecord {
        // record errors during processing
        // not using try-catch for performance consideration
    private:
        std::string errorFileName_;
        std::ofstream errorFile_;
    public:
        void setParam(const RefactorProcessing::Parameters& P);
        void reportError(const std::string& errorMsg);
    };
    class WarningRecord {
        // record warnings during processing
    private:
        std::string warningFileName_;
        std::ofstream warningFile_;
    public:
        void setParam(const RefactorProcessing::Parameters& P);
        void reportWarning(const std::string& warningMsg);
    };
}
#endif //RNAALIGNREFACTORED_ERRORRECORD_H
