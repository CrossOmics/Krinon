#ifndef RNAALIGNREFACTORED_LOADREADFILER_H
#define RNAALIGNREFACTORED_LOADREADFILER_H
#include <string>
#include "Parameters.h"
#include "../readAlignRefactored/Read.h"
namespace RefactorProcessing {

    class ReadFileLoader {
        // get and decode read files
    private:
        // configs
        std::string readFileType_; // fastq, fasta, etc.
        std::string readType_; // single-end or paired-end, maybe long reads in the future
        std::string fileName1_;
        std::string fileName2_; // for paired-end

    public:
        ReadFileLoader(){};
        ~ReadFileLoader(){};
        void setParam(const Parameters& P);
        void loadRead(Read &read);
    };
}
#endif //RNAALIGNREFACTORED_LOADREADFILER_H
