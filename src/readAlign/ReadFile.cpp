#include "ReadFile.h"
#include "../utils/exceptions.h"
#include <sstream>
#include <algorithm>
#include "../utils/Parameters.h"

#define MAX_BUFFER_READ_NUM 10000
#define MAX_READ_CHUNK_SIZE 50000

namespace rna {
    ReadFile::ReadFile(std::string &readType) {

        if (readType == "single") {
            this->readType = single;
        } else if (readType == "paired") {
            this->readType = paired;
        }
        readFileName = "";
        readFileName2 = "";

    }

    ReadFile::ReadFile(Parameters &P) {

        if (P.isPaired) {
            this->readType = paired;
            readFileName = P.readFile;
            readFileName2 = P.readFile2;

        } else {
            this->readType = single;
            readFileName = P.readFile;
        }
    }

    void ReadFile::openFiles(const std::string &filename1, const std::string &filename2) {
        int64_t inputBufferSize = 32*1024*1024;
        readFileBuffer = new char[inputBufferSize];
        readFile = std::ifstream(filename1);
        readFile.rdbuf()->pubsetbuf(readFileBuffer, inputBufferSize);

        if (readType == paired) {
            readFileBuffer2 = new char[inputBufferSize];
            readFile2 = std::ifstream(filename2);
            readFile2.rdbuf()->pubsetbuf(readFileBuffer2, inputBufferSize);

        }
    }

    void ReadFile::closeFiles() {
        if(readFile.is_open()) readFile.close();
        if (readType == paired) {
            if(readFile2.is_open()) readFile2.close();
        }
    }




    int ReadFile::loadReadFromFastq(std::vector<std::array<std::string, 4>>& lines1,std::vector<std::array<std::string, 4>>& lines2,int chunkSize) {


        int r;
        fileLock.lock();
        
        bool isEnd = false;
        for (r = 0; r < chunkSize; ++r) {
            auto &l1 = lines1[r];
            auto &l2 = lines2[r];
            for (int i = 0; i < 4; ++i) {
                if (!(std::getline(readFile, l1[i]))) {
                    isEnd = true;
                    break;
                }
            }
            if (readType == paired) {
                for (int i = 0; i < 4; ++i) {
                    if (!(std::getline(readFile2, l2[i]))) {
                        isEnd = true;
                        break;
                    }
                }
            }
            if (isEnd) break;
        }
        fileLock.unlock();




        return r;


    }

    /*int ReadFile::loadReadChunkFromFastq(std::vector<ReadPtr>& reads, int chunkSize){
        int count = 0;
        fileLock.lock();
        for(int i = 0; i < chunkSize; ++i){
            if(!loadReadFromFastq(reads[i])){
                fileLock.unlock();
                return count;
            }
            ++count;
        }
        fileLock.unlock();
        return count;
    }*/
}