#ifndef RNAALIGNREFACTORED_READALIGNMULTITHREAD_H
#define RNAALIGNREFACTORED_READALIGNMULTITHREAD_H
#include <thread>
#include <mutex>
#include "ReadAligner.h"
#include "../utils/exceptions.h"

namespace rna{
    class ReadAlignMultiThread {
        std::mutex fileReadLock; // mutex for file reading
        std::mutex outputLock; // mutex for output operations
        std::mutex alignStatusLock; // mutex for alignment status output
        std::mutex alignProgressLock; // mutex for alignment progress output
        int threadNum;
        int readChunkSize{100};
        std::vector<std::thread> threads;
        int readCount{0};
        std::vector<ReadAligner> readAligners;
        std::ifstream file;
        std::ofstream outFile;
        std::ofstream alignStatusFile;
        std::ofstream alignProgressFile;

        void singleThreadProcess(ReadAligner&r) {

            r.processReadFile(file,outFile,alignStatusFile,alignProgressFile, outputLock, alignStatusLock,alignProgressLock,fileReadLock,readCount);
        }
    public:
        void processReadFile(const std::string& filename, int tNum,GenomeIndexPrefix& gInPre,bool partialOutput ) {
            threadNum = tNum;
            file= std::ifstream (filename);
            outFile = std::ofstream (filename+ ".sam");
            alignStatusFile = std::ofstream (filename + ".alignTime");
            alignProgressFile = std::ofstream (filename + ".progress");
            if (!file) {
                throw RNAException(
                        ExceptionType::FILE_NOT_FOUND,
                        "Cannot open read file: " + filename
                );
            }

            for (int i = 0; i<threadNum;++i){
                readAligners.emplace_back(gInPre,i);
                readAligners[i].partialOutput = partialOutput;

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
            file.close();
            outFile.close();
            alignStatusFile.close();
            alignProgressFile.close();

        }

    };
}

#endif //RNAALIGNREFACTORED_READALIGNMULTITHREAD_H
