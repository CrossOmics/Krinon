
#ifndef RNAALIGNREFACTORED_PARAMETERSR_H
#define RNAALIGNREFACTORED_PARAMETERSR_H
#include <string>
namespace RefactorProcessing {
    class Parameters {
        // decode and store parameters from command line or config file
    private:

    public:
        Parameters();
        ~Parameters();
        std::string outErrorFile;
        std::string outWarningFile;

    };
}
#endif //RNAALIGNREFACTORED_PARAMETERSR_H
