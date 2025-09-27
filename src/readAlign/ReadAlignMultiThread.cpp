#include "ReadAlignMultiThread.h"
#include "../utils/Parameters.h"
#include <filesystem>
#include "../utils/utils.h"

namespace rna {
    ReadAlignMultiThread::ReadAlignMultiThread(Parameters& params):
        threadNum(params.threads),
        m_inputFile1(params.readFile),
        m_inputFile2(params.readFile2),
        m_generationCounter(params.threads)
    {
    }

    void ReadAlignMultiThread::singleThreadProcess(ReadAligner&r, int threadIndex, std::string outputDir) {
        // r.processReadFile(readFile, outFile, logFile, alignProgressFile, outputLock, alignStatusLock, alignProgressLock, readCount);
        // r.processReadFile(readFile, readCount, threadIndex);
        auto alignerOutputFileName = (std::filesystem::path(outputDir) / std::filesystem::path("alignerOutput-" + std::to_string(threadIndex) + ".tmp")).string();
        r.processReadFile(alignerOutputFileName, readCount, threadNum, threadIndex, m_alignProgressLock, m_generationCounter);
    }
          
    void ReadAlignMultiThread::processReadFile(rna::Parameters& params, rna::GenomeIndex &gInPre, bool partialOutput) {
        // int64_t outputBufferSize = 50000000 * threadNum;  
        // outputAlignBuffer = new char[outputBufferSize];
        // std::string filename = readFile.readFileName;
        // readFile.openFiles(readFile.readFileName, readFile.readFileName2);
        // readFile.openFiles(readFile.readFileName);

        // setvbuf(outFile,outputAlignBuffer,_IOFBF,outputBufferSize);
        // logFile = std::ofstream (outDir + "outLog.out");
        // alignProgressFile = std::ofstream (outDir + "outLog.progress.out");
        // alignProgressFile << "v8\n";
        // alignProgressFile << "Started at: " << getTime() << '\n';

        // TODO@arvin: Handle pairs later
        assert(m_inputFile2 == "");

        // First, check how large the file is
        int64_t totalBytesToRead = std::filesystem::file_size(m_inputFile1);

        // Now, chunk it up given the number of threads
        int64_t threadChunkSize = totalBytesToRead / threadNum;
        int64_t remainder = totalBytesToRead % threadNum;
        int64_t start = 0;
        int64_t currentChunkSize = 0;

        for (int i = 0; i < threadNum; ++i){
            if (i == 0)
                currentChunkSize = threadChunkSize + remainder;
            else
                currentChunkSize = threadChunkSize;
            readAligners.push_back(ReadAligner(params, gInPre, start, currentChunkSize, i));
            readAligners[i].partialOutput = partialOutput;
            readAligners[i].setConfig(stitchConfig, stitchingScoreConfig, seedMappingConfig);
            start += currentChunkSize;
        }

        assert(start == totalBytesToRead);

        auto alignmentStart = std::chrono::high_resolution_clock::now();
        
        for (int i = 0; i < threadNum; ++i) {
            auto &r = readAligners[i];
            threads.emplace_back([this, &r, &params, i] { singleThreadProcess(r, i, params.outPutDir); });
        }
        for (auto &t: threads) {
            if (t.joinable()) {
                t.join();
            }
        }
        threads.clear();

        auto now = std::chrono::high_resolution_clock::now();
        auto totalTime = std::chrono::duration_cast<std::chrono::seconds>(now - alignmentStart).count();
        PLOG_INFO << "Total time for alignment: " << totalTime << " seconds";
        // readFile.closeFiles();
        // fclose(outFile);
        // logFile << "Total reads processed: " << readFile.totalReadCount << '\n';
        // logFile << "Unique reads: " << readFile.uniqueReadCount << '\n';
        // logFile << "Multi-mapping reads: " << readFile.multiReadCount << '\n';
        // logFile << "Unique mapping rate: " << std::setprecision(4) << double(readFile.uniqueReadCount) / double(readFile.totalReadCount) * 100.0 << "%\n";
        // logFile << "Multi-mapping rate: " << std::setprecision(4) << double(readFile.multiReadCount) / double(readFile.totalReadCount) * 100.0 << "%\n";

        // logFile.close();
        // alignProgressFile.close();

        // auto concatStart = std::chrono::high_resolution_clock::now();
        // auto outputSamFile = (std::filesystem::path(params.outPutDir) / std::filesystem::path("alignerOutput.out.sam")).string();
        // auto alignerOutputsWildcard = (std::filesystem::path(params.outPutDir) / std::filesystem::path("alignerOutput-*.tmp")).string();
        // PLOG_INFO << "Concatenating the files into " << outputSamFile;
        // system(("cat " + alignerOutputsWildcard + " > " + outputSamFile).c_str());
        // now = std::chrono::high_resolution_clock::now();
        // auto concatTime = std::chrono::duration_cast<std::chrono::seconds>(now - concatStart).count();
        // PLOG_INFO << "Output concat. done in " << concatTime << " seconds.";

        // auto cleanupStart = std::chrono::high_resolution_clock::now();
        // PLOG_INFO << "Removing aligner .tmp files ";
        // system(("rm " + alignerOutputsWildcard).c_str());
        // now = std::chrono::high_resolution_clock::now();
        // auto cleanupTime = std::chrono::duration_cast<std::chrono::seconds>(now - cleanupStart).count();
        // PLOG_INFO << "Output cleanup done in " << cleanupTime << " seconds.";
    }
} // namespace rna