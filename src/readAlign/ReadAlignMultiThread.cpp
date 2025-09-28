#include "ReadAlignMultiThread.h"
#include "../utils/Parameters.h"
#include <filesystem>
#include <sys/mman.h>
#include "../utils/utils.h"

namespace rna {
    ReadAlignMultiThread::ReadAlignMultiThread(Parameters& params):
        threadNum(params.threads),
        m_inputFile1(params.readFile),
        m_inputFile2(params.readFile2),
        m_generationCounter(params.threads)
    {}

    void ReadAlignMultiThread::singleThreadProcess(ReadAligner& r, int threadIndex, MemoryMappedFile& output) {
        r.processReadFile(output, m_outputOffset, m_totalReadCount, m_generationCounter);
    }

    void ReadAlignMultiThread::initiateAligners(
        rna::Parameters& params, 
        rna::GenomeIndex &gInPre,
        int64_t totalBytesToRead, 
        bool partialOutput) 
    {
        int64_t start = 0;
        int64_t currentChunkSize = 0;
        int64_t threadChunkSize = totalBytesToRead / threadNum;
        int64_t remainder = totalBytesToRead % threadNum;

        for (int i = 0; i < threadNum; ++i){
            if (i == 0)
                currentChunkSize = threadChunkSize + remainder;
            else
                currentChunkSize = threadChunkSize;
            m_readAligners.push_back(ReadAligner(params, gInPre, start, currentChunkSize, i));
            m_readAligners[i].partialOutput = partialOutput;
            m_readAligners[i].setConfig(stitchConfig, stitchingScoreConfig, seedMappingConfig);
            start += currentChunkSize;
        }

        assert(start == totalBytesToRead);
    }

    void ReadAlignMultiThread::doMultiThreadedAlignment(
        rna::Parameters& params,
        MemoryMappedFile& output) 
    {
        auto alignmentStart = std::chrono::high_resolution_clock::now();
        
        for (int i = 0; i < threadNum; ++i) {
            auto &r = m_readAligners[i];
            m_threads.emplace_back([this, &r, &output, i] { singleThreadProcess(r, i, output); });
        }
        for (auto &t: m_threads) {
            if (t.joinable())
                t.join();
        }
        m_threads.clear();

        auto now = std::chrono::high_resolution_clock::now();
        auto totalTime = std::chrono::duration_cast<std::chrono::seconds>(now - alignmentStart).count();
        PLOG_INFO << "Total time for alignment: " << totalTime << " seconds";
    }
          
    void ReadAlignMultiThread::processReadFile(rna::Parameters& params, rna::GenomeIndex &gInPre, bool partialOutput) {
        /**
         * TODO: Handle pairs later
         */
        assert(m_inputFile2 == "");

        // First, check how large the file is
        int64_t totalBytesToRead = std::filesystem::file_size(m_inputFile1);

        // Allocate the memory-mapped output file
        // Currently, I have hardcoded it to be 1.5 times the input size
        int64_t maximumTotalBytesToWrite = 150 * totalBytesToRead / 100;
        auto outputFilePath = 
            std::filesystem::path(params.outPutDir) / 
            std::filesystem::path("alignerOutput.out.sam");
        MemoryMappedFile outputFile(outputFilePath, maximumTotalBytesToWrite);

        // Now, chunk the input up given the number of threads
        initiateAligners(params, gInPre, totalBytesToRead, partialOutput);
        // Being alignment
        doMultiThreadedAlignment(params, outputFile);
        // Close the output
        outputFile.memClose();

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