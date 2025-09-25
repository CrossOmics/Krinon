#include "ReadAlignMultiThread.h"
#include "../utils/Parameters.h"
#include <filesystem>
#include "../utils/utils.h"
namespace rna {
    ReadAlignMultiThread::ReadAlignMultiThread(Parameters& P):
        readFile(P){
        if (!std::filesystem::exists(P.outPutDir)) {
            std::filesystem::create_directories(P.outPutDir);
        }
        outDir = P.outPutDir;
        outFile = fopen((P.outPutDir+"outAligned.out.sam").c_str(),"w");
        
    }
    void ReadAlignMultiThread::processReadFile(int tNum, rna::GenomeIndex &gInPre,
                                               bool partialOutput) {
        threadNum = tNum;
        int64_t outputBufferSize = 50000000 * threadNum; 
        outputAlignBuffer = new char[outputBufferSize];
        std::string filename = readFile.readFileName;
        readFile.openFiles(readFile.readFileName, readFile.readFileName2);

        setvbuf(outFile,outputAlignBuffer,_IOFBF,outputBufferSize);
        logFile = std::ofstream (outDir + "outLog.out");
        alignProgressFile = std::ofstream (outDir + "outLog.progress.out");
        alignProgressFile << "v7\n";
        alignProgressFile << "Started at: " << getTime() << '\n';


        for (int i = 0; i<threadNum;++i){

            readAligners.push_back(ReadAligner(gInPre,i));
            readAligners[i].partialOutput = partialOutput;
            readAligners[i].setConfig(stitchConfig, stitchingScoreConfig, seedMappingConfig);
        }
        for (int i = 0; i < threadNum; ++i) {
            auto &r = readAligners[i];
            threads.emplace_back([this, &r ] { singleThreadProcess(r); });
        }

        for (auto &t: threads) {

            if (t.joinable()) {
                t.join();
            }
        }
        threads.clear();
        readFile.closeFiles();
        fclose(outFile);
        logFile << "Total reads processed: " << readFile.totalReadCount << '\n';
        logFile << "Unique reads: " << readFile.uniqueReadCount << '\n';
        logFile << "Multi-mapping reads: " << readFile.multiReadCount << '\n';
        logFile << "Unique mapping rate: " << std::setprecision(4) << double(readFile.uniqueReadCount) / double(readFile.totalReadCount) * 100.0 << "%\n";
        logFile << "Multi-mapping rate: " << std::setprecision(4) << double(readFile.multiReadCount) / double(readFile.totalReadCount) * 100.0 << "%\n";

        logFile.close();
        alignProgressFile.close();

    }
}