#ifndef RNAALIGNREFACTORED_READFILE_H
#define RNAALIGNREFACTORED_READFILE_H
#include "../utils/types.h"
#include <fstream>
#include <mutex>
#include <memory>

namespace rna {
    class Parameters;
    class ReadFile {
    private:
        char* readFileBuffer{nullptr};
        char* readFileBuffer2{nullptr};
    public:
        enum ReadType{
            single,paired
        };
        std::string readFileName;
        std::string readFileName2;
        std::ifstream readFile;
        std::ifstream readFile2; // for paired-end
        std::mutex fileLock;
        ReadType readType{single};
        int64_t uniqueReadCount{0};
        int64_t multiReadCount{0};
        int64_t totalReadCount{0};


        ReadFile(std::string &readType);
        ReadFile(Parameters &P);

        ~ReadFile(){
            delete[] readFileBuffer;
            delete[] readFileBuffer2;
        };
        void openFiles(const std::string& filename1,const std::string& filename2="");
        bool loadReadFromFastq(ReadPtr& read);
        void closeFiles();
    };
}
#endif //RNAALIGNREFACTORED_READFILE_H
