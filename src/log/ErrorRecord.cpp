#include "ErrorRecord.h"
#include <fstream>
#include "../utilsRefactored//timeFunctions.hpp"


namespace rna{
    void ErrorRecord::setParam(const RefactorProcessing::Parameters& P){
        errorFileName_ = P.outErrorFile;
        errorFile_ = std::ofstream(errorFileName_, std::ios::app);
    }

    void ErrorRecord::reportError(const std::string& errorMsg){
        errorFile_<<getCurrentTimeString()<<" "<< errorMsg << std::endl;
        //todo free resources

    }

    void WarningRecord::setParam(const RefactorProcessing::Parameters& P){
        warningFileName_ = P.outWarningFile;
        warningFile_ = std::ofstream(warningFileName_, std::ios::app);
    }

    void WarningRecord::reportWarning(const std::string& warningMsg){
        warningFile_<<getCurrentTimeString()<<" " << warningMsg << std::endl;
    }
}