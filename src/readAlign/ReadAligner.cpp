#include "ReadAligner.h"
#include "../utils/exceptions.h"
#include <fstream>
#include <filesystem>
#include <chrono>
#include <mutex>
#include <sstream>
#include "ReadFile.h"
#include <iomanip>
#include <csignal>

namespace rna {

    inline void parseRead(ReadPtr &read,std::array<std::string,4>& lines1,std::array<std::string,4>& lines2,bool isPaired){
        read = std::make_shared<Read>();
        for (int i = 0;i<4;++i){
            if (lines1[i][lines1[i].size() - 1] == '\r') {
                    lines1[i].pop_back(); // Windows file in linux
            }
            if (lines2[i][lines2[i].size() - 1] == '\r') {
                    lines2[i].pop_back(); // Windows file in linux
            }

        }
        std::istringstream nameStream(lines1[0].substr(1));
        nameStream >> read->name;
        if(isPaired){
            std::reverse(lines2[1].begin(),lines2[1].end());
            std::reverse(lines2[3].begin(),lines2[3].end());
            for (char &c: lines2[1]) {
                switch (c) {
                    case 'A':
                        c = 'T';
                        break;
                    case 'T':
                        c = 'A';
                        break;
                    case 'C':
                        c = 'G';
                        break;
                    case 'G':
                        c = 'C';
                        break;
                    default:
                        c = 'N';
                        break;
                }
            }
            read->sequence[0] = lines1[1] + '#' + lines2[1];
            read->sequence[1] = lines1[1] + '#' + lines2[1];
            read->length = lines1[1].length() + 1 + lines2[1].length();
            read->quality = lines1[3] + ' ' + lines2[3];
            read->mate1Length = lines1[1].length();
            read->mate2Length = lines2[1].length();
        }else {
            read->sequence[0] = lines1[1];
            read->sequence[1] = lines1[1];
            read->length = lines1[1].length();
            read->quality = lines1[3];
        }

        std::reverse(read->sequence[1].begin(), read->sequence[1].end());
        for (char &c: read->sequence[1]) {
            switch (c) {
                case 'A':
                    c = 'T';
                    break;
                case 'T':
                    c = 'A';
                    break;
                case 'C':
                    c = 'G';
                    break;
                case 'G':
                    c = 'C';
                    break;
                case '#':
                    break;
                default:
                    c = 'N';
                    break;
            }
        }
    }

    bool ReadAligner::loadReadFromFastq(ReadFile& file) {
        if(queueEmpty) {
            nowReadInd = 0;
            int cnt = file.loadReadFromFastq(lines1,lines2, inputBufferSize);

            if (cnt == 0) {
                return false;
            }
            nowQueueSize = cnt;
            queueEmpty = false;
        }
        read = std::make_shared<Read>();
        parseRead(read,lines1[nowReadInd],lines2[nowReadInd],file.readType == ReadFile::paired);
        ++nowReadInd;
        if(nowReadInd >= nowQueueSize) queueEmpty = true;
        return true;

    }
    void ReadAligner::processReadFile(ReadFile& file,FILE* outFile,std::ofstream& logFile,std::ofstream& alignProgressFile,std::mutex& outputLock,std::mutex& alignStatusLock,std::mutex& alignProgressLock,int& totalReadsProcessed) {

        std::stringstream outputBuffer((std::string()));
        int outputBufferSize = 10000;
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
                            std::string outputStr = outputBuffer.str();
                            int64_t outputLen = outputStr.length();
                            outputLock.lock();
                            fprintf(outFile, outputStr.c_str(), outputLen);
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



            

            
            if (isUnique) file.threadUniqueReadCount[threadId]++;
            if (isMulti) file.threadMultiReadCount[threadId]++;
            file.threadReadCount[threadId]++;
            if (threadId == 0){
                int64_t nowTotalRead = 0;
                for (int i = 0; i<file.threadNum; ++i) nowTotalRead += file.threadReadCount[threadId];
                auto nowTime = std::chrono::high_resolution_clock::now();
                if (std::chrono::duration_cast<std::chrono::seconds>(nowTime - previousProgressReportTime).count() >= 60){
                previousProgressReportTime = nowTime;
                auto totalTime = std::chrono::duration_cast<std::chrono::seconds>(nowTime - AligningStartTime).count();
                alignProgressFile << "Current time:" << totalTime << "s\t";
                alignProgressFile << "Total reads processed (all threads): " << nowTotalRead << '\t';
                alignProgressFile << "Speed:" <<std::fixed << std::setprecision(1)<< double (nowTotalRead) / double (totalTime) * 3600.0/1000000 << "Million r/h\n";
                alignProgressFile.flush();
                }
            }
            


        }

        if (outputBufferCnt > 0) {
            outputLock.lock();
            fprintf(outFile, outputBuffer.str().c_str(), outputBuffer.str().length());
            outputLock.unlock();
            outputBuffer.str(std::string());
            outputBufferCnt = 0;
        }

         if (threadId == 0){

                int64_t nowTotalRead = 0;
                for (int i = 0; i<file.threadNum; ++i){
                    nowTotalRead += file.threadReadCount[i];
                }
                auto nowTime = std::chrono::high_resolution_clock::now();
                
                previousProgressReportTime = nowTime;
                auto totalTime = std::chrono::duration_cast<std::chrono::seconds>(nowTime - AligningStartTime).count();
                alignProgressFile << "Current time:" << totalTime << "s\t";
                alignProgressFile << "Total reads processed (all threads): " << nowTotalRead << '\t';
                alignProgressFile << "Speed:" <<std::fixed << std::setprecision(1)<< double (nowTotalRead) / double (totalTime) * 3600.0/1000000 << "Million r/h\n";
                alignProgressFile.flush();
                
        }







    }


}