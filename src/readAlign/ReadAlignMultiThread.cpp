#include "ReadAlignMultiThread.h"
#include "../utils/Parameters.h"
#include <filesystem>
#include <fcntl.h>
#include "../utils/utils.h"
#include <chrono>
namespace rna {
    ReadAlignMultiThread::ReadAlignMultiThread(Parameters& P):
        readFile(P){
        if (!std::filesystem::exists(P.outPutDir)) {
            std::filesystem::create_directories(P.outPutDir);
        }
        outDir = P.outPutDir;
        //outFile = fopen((P.outPutDir+"outAligned.out.sam").c_str(),"w");
        outFile = open((P.outPutDir+"outAligned.out.sam").c_str(), O_WRONLY | O_CREAT, 0644);


    }
    void ReadAlignMultiThread::processReadFile(int tNum, rna::GenomeIndex &gInPre,
                                               bool partialOutput) {
        threadNum = tNum;
        int64_t outputBufferSize = 50*1024*1024;
        //outputAlignBuffer = new char[outputBufferSize];

        std::string filename = readFile.readFileName;
        readFile.openFiles(readFile.readFileName, readFile.readFileName2);
        readFile.threadReadCount = new int64_t[threadNum];
        readFile.threadUniqueReadCount = new int64_t[threadNum];
        readFile.threadMultiReadCount = new int64_t[threadNum];
        readFile.threadNum = threadNum;


        //setvbuf(outFile,outputAlignBuffer,_IOFBF,outputBufferSize);
        logFile = std::ofstream (outDir + "outLog.out");

        alignProgressFile = std::ofstream (outDir + "outLog.progress.out");
        alignProgressFile << "v10\n";
        alignProgressFile << "Started at: " << getTime() << '\n';


        for (int i = 0; i<threadNum;++i){
            readFile.threadReadCount[i] = 0;
            readFile.threadMultiReadCount[i] = 0;
            readFile.threadUniqueReadCount[i] = 0;

            readAligners.push_back(ReadAligner(gInPre,i));
            readAligners[i].partialOutput = partialOutput;
            readAligners[i].setConfig(stitchConfig, stitchingScoreConfig, seedMappingConfig);
        }
        for (int i = 0; i < threadNum; ++i) {
            auto &r = readAligners[i];
            threads.emplace_back([this, &r ] { singleThreadProcess(r); });
        }

        auto startTime = std::chrono::high_resolution_clock::now();
        for (auto &t: threads) {

            if (t.joinable()) {
                t.join();
            }
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        auto totalTime = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
        threads.clear();
        readFile.closeFiles();
        for (int i = 0; i<threadNum;++i){
            readFile.totalReadCount += readFile.threadReadCount[i];
            readFile.uniqueReadCount += readFile.threadUniqueReadCount[i];
            readFile.multiReadCount += readFile.threadMultiReadCount[i];
        }
        logFile << "Total reads processed: " << readFile.totalReadCount << '\n';
        logFile << "Unique reads: " << readFile.uniqueReadCount << '\n';
        logFile << "Multi-mapping reads: " << readFile.multiReadCount << '\n';
        logFile << "Unique mapping rate: " << std::setprecision(4) << double(readFile.uniqueReadCount) / double(readFile.totalReadCount) * 100.0 << "%\n";
        logFile << "Multi-mapping rate: " << std::setprecision(4) << double(readFile.multiReadCount) / double(readFile.totalReadCount) * 100.0 << "%\n";
        logFile << "Average Speed: " << std::fixed << std::setprecision(1)<< double (readFile.totalReadCount) / double (totalTime) * 3600.0/1000000 << "Million r/h\n";
        logFile.close();
        alignProgressFile.close();

    }
}