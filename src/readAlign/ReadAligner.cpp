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
    ReadAligner::ReadAligner(
        const Parameters& P, 
        const GenomeIndex& gInPre, 
        int64_t beginFrom, 
        int64_t expectedReadableBytes, 
        int threadIndex
    ) : 
        genomeIndexPrefix(gInPre), 
        m_beginFrom(beginFrom),
        m_expectedReadableBytes(expectedReadableBytes),
        m_threadIndex(threadIndex) 
    {
        seedMapping = std::make_unique<SeedMapping>( genomeIndexPrefix,SeedMappingConfig());
        stitchingManagement = std::make_unique<StitchingManagement>(StitchingConfig(), gInPre);

        /**
         * TODO: Inherit this from the given parameter struct
         */
        inputBufferSize = 30000;

        readBufferQueue.resize(inputBufferSize);
        for (int i = 0; i < inputBufferSize; ++i) {
            readBufferQueue[i] = std::make_shared<Read>();
        }

        // Open the file and prepare for reading
        m_readFile = std::make_unique<ReadFile>(P, beginFrom, expectedReadableBytes, threadIndex);
        m_readFile->openFiles(P.readFile, P.readFile2);
    }

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

    void ReadAligner::reportProgress(
        std::atomic_int64_t& totalReadsProcessed,
        AtomicGenerationCounter& generationCounter) 
    {
        auto nowTime = std::chrono::high_resolution_clock::now();
        if (
            std::chrono::duration_cast<std::chrono::seconds>(
                nowTime - m_previousProgressReportTime).count() >= RNAALIGNER_PROGRESS_REPORT_INTERVAL_SECONDS
            )
        {
            // Decrement the generation and add the number of things you read
            // If generation has hit zero, report the numbers
            m_previousProgressReportTime = nowTime;
            auto currentTotal = totalReadsProcessed.fetch_add(m_periodicReadCount) + m_periodicReadCount;
            m_periodicReadCount = 0;
            if (generationCounter.decrement()) {
                auto totalTime = std::chrono::duration_cast<std::chrono::seconds>(nowTime - m_alignmentStartTime).count();
                PLOG_INFO << "Current time: " << totalTime << " s";
                PLOG_INFO << "Total reads processed (all threads): " << currentTotal;
                PLOG_INFO << "Speed: " << std::fixed 
                          << std::setprecision(1) 
                          << (double (currentTotal) / double (totalTime)) * (3600 / 1e6)
                          << " Million r/h";
            }
        }
    }

    void ReadAligner::processReadFile(
        MemoryMappedFile& outputFileMapped,
        std::atomic_int64_t& outputOffset,
        std::atomic_int64_t& totalReadsProcessed, 
        AtomicGenerationCounter& generationCounter) 
    {
        PLOG_DEBUG << "Aligner thread " << m_threadIndex << " spawned"; 
        std::stringstream outputBuffer((std::string()));

        // TODO: Tune this!!!
        int outputBufferSize = 10000 + m_threadIndex * 100;
        int outputBufferCnt = 0;

        // count the number of reads processed in this call
        m_alignmentStartTime = std::chrono::high_resolution_clock::now();
        m_previousProgressReportTime = m_alignmentStartTime;

        while (loadReadFromFastq()) {
            bool isUnique = false;
            bool isMulti = false;

            seedMapping->processRead(read);
            stitchingManagement->processAlignments(seedMapping->aligns,seedMapping->alignNum,read);

            if (stitchingManagement->numGoodTranscripts_ == 1) {
                isUnique = true;
            } else if (stitchingManagement->numGoodTranscripts_ > 1){
                isMulti = true;
            }

            /**
             * TODO: Bring out that magic number!
             */
            
            if(!partialOutput || totalReadsProcessed.load() < 10000){
                if (stitchingManagement->status == StitchingManagement::SUCCESS){
                    for (int j = 0; j < stitchingManagement->numGoodTranscripts_; ++j) {
                        auto &t = stitchingManagement->goodTranscripts_[j];
                        t.isPaired = (m_readFile->readType == ReadFile::paired);
                        std::string s = t.outputSam(*read, stitchingManagement->numGoodTranscripts_);
                        outputBuffer << s;
                        ++outputBufferCnt;
                        if (outputBufferCnt >= outputBufferSize) {
                            auto currentString = outputBuffer.str();
                            auto currentOffset = outputOffset.fetch_add(currentString.length());
                            std::memcpy(outputFileMapped.getMapPtr() + currentOffset, currentString.c_str(), currentString.length());
                            outputBuffer.str(std::string());
                            outputBufferCnt = 0;
                        }
                    }
                }
            }

            stitchingManagement->clear();

            seedMapping->clear();
            ++m_periodicReadCount;

            // if (isUnique) file.uniqueReadCount++;
            // if (isMulti) file.multiReadCount++;
            reportProgress(totalReadsProcessed, generationCounter);
        }

        if (outputBufferCnt > 0) {
            auto currentString = outputBuffer.str();
            auto currentOffset = outputOffset.fetch_add(currentString.length());
            std::memcpy(outputFileMapped.getMapPtr() + currentOffset, currentString.c_str(), currentString.length());
        }
        m_readFile->closeFiles();
        PLOG_DEBUG << "Aligner thread " << m_threadIndex << " finished";
    }
} // namespace rna