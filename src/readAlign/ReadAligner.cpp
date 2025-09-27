#include "ReadAligner.h"
#include "../utils/exceptions.h"
#include <fstream>
#include <filesystem>
#include <chrono>
#include <mutex>
#include <sstream>
#include "ReadFile.h"
#include <iomanip>
#include <plog/Log.h>

namespace rna {
    ReadAligner::ReadAligner(const Parameters& P, const GenomeIndex& gInPre, int64_t beginFrom, int64_t expectedReadableBytes, int tId) : 
        genomeIndexPrefix(gInPre), 
        m_beginFrom(beginFrom),
        m_expectedReadableBytes(expectedReadableBytes),
        threadId(tId) 
    {
        seedMapping = std::make_unique<SeedMapping>( genomeIndexPrefix,SeedMappingConfig());
        stitchingManagement = std::make_unique<StitchingManagement>(StitchingConfig(), gInPre);
        inputBufferSize = 30000;
        // inputBufferSize = 100000;
        readBufferQueue.resize(inputBufferSize);
        for (int i = 0; i < inputBufferSize; ++i) {
            readBufferQueue[i] = std::make_shared<Read>();
        }

        // Open the file and prepare for reading
        m_readFile = std::make_unique<ReadFile>(P, beginFrom, expectedReadableBytes, tId);
        m_readFile->openFiles(P.readFile, P.readFile2);
    }

    // bool ReadAligner::loadReadFromFastq(ReadFile& file) {
    bool ReadAligner::loadReadFromFastq() {
        if(queueEmpty) {
            nowReadInd = 0;
            int cnt = m_readFile->loadReadChunkFromFastq(readBufferQueue, inputBufferSize);

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

    // void ReadAligner::processReadFile(ReadFile& file,FILE* outFile,std::ofstream& logFile,std::ofstream& alignProgressFile,std::mutex& outputLock,std::mutex& alignStatusLock,std::mutex& alignProgressLock,int& totalReadsProcessed) {
    // void ReadAligner::processReadFile(ReadFile& file, int& totalReadsProcessed, const int threadIndex) {
    void ReadAligner::processReadFile(
        std::string outFilePath, int& totalReadsProcessed, const int numThreads, const int threadIndex, std::mutex& progressLock,
        int& generationCounter) 
    {
        PLOG_DEBUG << "Aligner thread " << threadIndex << " spawned (producing into " << outFilePath << ")";
        
        auto outFile = fopen(outFilePath.c_str(),"w");
        std::stringstream outputBuffer((std::string()));

        // TODO: Tune this!!!
        int outputBufferSize = 10000;
        int outputBufferCnt = 0;

        // count the number of reads processed in this call
        auto AligningStartTime = std::chrono::high_resolution_clock::now();
        auto previousProgressReportTime = AligningStartTime;

        auto chunkCounts = 0;

        while (loadReadFromFastq()) {
            bool isUnique = false;
            bool isMulti = false;
            ++m_oneMinuteCycleReadCount;
            // auto startTime = std::chrono::high_resolution_clock::now();
            seedMapping->processRead(read);

            // auto alignEndTime = std::chrono::high_resolution_clock::now();

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
                        t.isPaired = (m_readFile->readType == ReadFile::paired);
                        std::string s = t.outputSam(*read, stitchingManagement->numGoodTranscripts_);
                        outputBuffer << s;
                        ++outputBufferCnt;
                        if (outputBufferCnt >= outputBufferSize) {
                            fprintf(outFile, outputBuffer.str().c_str(), outputBuffer.str().length());
                            outputBuffer.str(std::string());
                            outputBufferCnt = 0;
                        }
                    }
                }
            }

            stitchingManagement->clear();
            // auto stitchEndTime = std::chrono::high_resolution_clock::now();

            seedMapping->clear();

            // file.totalReadCount ++;
            // if (isUnique) file.uniqueReadCount++;
            // if (isMulti) file.multiReadCount++;
            auto nowTime = std::chrono::high_resolution_clock::now();
            if (std::chrono::duration_cast<std::chrono::seconds>(nowTime - previousProgressReportTime).count() >= 10){
                progressLock.lock();
                // Decrement the generation and add the number of things you read
                // If generation has hit zero, report the numbers
                --generationCounter;
                totalReadsProcessed += m_oneMinuteCycleReadCount;
                m_oneMinuteCycleReadCount = 0;
                previousProgressReportTime = nowTime;
                if (generationCounter == 0) {
                    generationCounter = numThreads;
                    auto totalTime = std::chrono::duration_cast<std::chrono::seconds>(nowTime - AligningStartTime).count();
                    PLOG_INFO << "Current time:" << totalTime << "s";
                    PLOG_INFO << "Total reads processed (all threads): " << totalReadsProcessed;
                    // PLOG_INFO << "Speed:" <<std::fixed << std::setprecision(1)<< double (totalReadsProcessed) / double (totalTime) * 3600.0/1000000 << "Million r/h";
                    PLOG_INFO << "Speed:" <<std::fixed << std::setprecision(1)<< double (totalReadsProcessed) / double (totalTime) * 360.0/100000 << "Million r/h";
                }
                progressLock.unlock();
            }

            // auto nowTime = std::chrono::high_resolution_clock::now();
            // auto totalTime = std::chrono::duration_cast<std::chrono::seconds>(nowTime - AligningStartTime).count();
            chunkCounts += 1;
            // PLOG_INFO << "Aligner thread " << threadIndex << " finished chunk in: " << totalTime << " seconds";
        }

        if (outputBufferCnt > 0) {
            fprintf(outFile, outputBuffer.str().c_str(), outputBuffer.str().length());
            outputBuffer.str(std::string());
            outputBufferCnt = 0;
        }
        m_readFile->closeFiles();
        fclose(outFile);
        PLOG_DEBUG << "Aligner thread " << threadIndex << " finished after processing " << chunkCounts << " chunks";
    }
} // namespace rna