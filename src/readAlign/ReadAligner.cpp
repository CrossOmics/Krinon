#include "ReadAligner.h"
#include "../utils/exceptions.h"
#include <fstream>
#include <filesystem>
#include <chrono>
#include <mutex>
#include <sstream>
#include "ReadFile.h"
#include <iomanip>
namespace rna {
    bool ReadAligner::loadReadFromFastq(ReadFile& file) {
        if(queueEmpty) {
            nowReadInd = 0;
            int cnt = file.loadReadChunkFromFastq(readBufferQueue, inputBufferSize);

            if (cnt == 0) {
                return false;
            }
            nowQueueSize = cnt;
            queueEmpty = false;
        }
        read = readBufferQueue[nowReadInd];
        ++nowReadInd;
        if(nowReadInd >= nowQueueSize) queueEmpty = true;
        return true;

    }
    void ReadAligner::processReadFile(ReadFile& file,FILE* outFile,std::ofstream& logFile,std::ofstream& alignProgressFile,std::mutex& outputLock,std::mutex& alignStatusLock,std::mutex& alignProgressLock,int& totalReadsProcessed) {

        std::stringstream outputBuffer((std::string()));
        int outputBufferSize = 50000;
        int outputBufferCnt = 0;


        // count the number of reads processed in this call
        auto AligningStartTime = std::chrono::high_resolution_clock::now();
        auto previousProgressReportTime = AligningStartTime;


        while (loadReadFromFastq(file)) {

            bool isUnique = false;
            bool isMulti = false;
            ++readCount;
            auto startTime = std::chrono::high_resolution_clock::now();
            seedMapping->processRead(read);

            auto alignEndTime = std::chrono::high_resolution_clock::now();

            stitchingManagement->processAlignments(seedMapping->aligns,seedMapping->alignNum,read);

            if (stitchingManagement->numGoodTranscripts_ == 1) {
                isUnique = true;
            } else if (stitchingManagement->numGoodTranscripts_ > 1){
                isMulti = true;
            }




            if(!partialOutput || totalReadsProcessed < 10000){
                if (stitchingManagement->status == StitchingManagement::SUCCESS){

                    for (int j = 0; j < stitchingManagement->numGoodTranscripts_; ++j) {
                        auto &t = stitchingManagement->goodTranscripts_[j];
                        t.isPaired = file.readType == ReadFile::paired;
                        std::string s = t.outputSam(*read, stitchingManagement->numGoodTranscripts_);
                        //fprintf(outFile,s.c_str(),s.length());
                        outputBuffer << s;
                        ++outputBufferCnt;
                        if (outputBufferCnt >= outputBufferSize) {
                            outputLock.lock();
                            fprintf(outFile, outputBuffer.str().c_str(), outputBuffer.str().length());
                            outputLock.unlock();
                            outputBuffer.str(std::string());
                            outputBufferCnt = 0;
                        }
                    }
                }
            }


            stitchingManagement->clear();
            auto stitchEndTime = std::chrono::high_resolution_clock::now();



            seedMapping->clear();



            alignProgressLock.lock();
            totalReadsProcessed++;
            file.totalReadCount ++;
            if (isUnique) file.uniqueReadCount++;
            if (isMulti) file.multiReadCount++;
            auto nowTime = std::chrono::high_resolution_clock::now();
            if (std::chrono::duration_cast<std::chrono::seconds>(nowTime - previousProgressReportTime).count() >= 60){
                previousProgressReportTime = nowTime;
                auto totalTime = std::chrono::duration_cast<std::chrono::seconds>(nowTime - AligningStartTime).count();
                alignProgressFile << "Current time:" << totalTime << "s\t";
                alignProgressFile << "Total reads processed (all threads): " << totalReadsProcessed << '\t';
                alignProgressFile << "Speed:" <<std::fixed << std::setprecision(1)<< double (totalReadsProcessed) / double (totalTime) * 3600.0/1000000 << "Million r/h\n";
            }

            alignProgressLock.unlock();
        }

        if (outputBufferCnt > 0) {
            outputLock.lock();
            fprintf(outFile, outputBuffer.str().c_str(), outputBuffer.str().length());
            outputLock.unlock();
            outputBuffer.str(std::string());
            outputBufferCnt = 0;
        }







    }


}