#ifndef RNAALIGNREFACTORED_READFILE_H
#define RNAALIGNREFACTORED_READFILE_H
#include "../utils/types.h"
#include <fstream>
#include <mutex>
#include <memory>
#include <chrono>
#include <atomic>

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


        //todo move them to a new class
        std::atomic<int64_t> uniqueReadCount{0};
        std::atomic<int64_t> multiReadCount{0};
        std::atomic<int64_t> totalReadCount{0};

        int64_t* threadUniqueReadCount;
        int64_t* threadMultiReadCount;
        int64_t* threadReadCount;
        int threadNum;




        ReadFile(std::string &readType);
        ReadFile(Parameters &P);

        ~ReadFile(){
            delete[] readFileBuffer;
            delete[] readFileBuffer2;
            delete[] threadMultiReadCount;
            delete[] threadUniqueReadCount;
            delete[] threadReadCount;
        };
        void openFiles(const std::string& filename1,const std::string& filename2="");
        int loadReadFromFastq(std::vector<std::array<std::string, 4>>& lines1,std::vector<std::array<std::string, 4>>& lines2,int chunkSize);
        //int loadReadChunk(char* buffer, int64_t bufferSize);
        //int loadReadChunkFromFastq(std::vector<ReadPtr>& reads, int chunkSize);
        void closeFiles();
    };
}
#endif //RNAALIGNREFACTORED_READFILE_H
