#include "ReadAligner.h"
#include "../utils/exceptions.h"
#include <fstream>
#include <filesystem>
#include <chrono>
#include <mutex>
namespace rna {
    bool ReadAligner::loadReadFromFastq(std::ifstream& s,std::mutex& readLock) {
        std::array<std::string, 4> lines;

        readLock.lock();
        for (int i = 0; i < 4; ++i) {
            if (!(std::getline(s, lines[i]))) {
                readLock.unlock();
                return false;
            }
            if (lines[i][lines[i].size() - 1] == '\r') {
                lines[i].pop_back(); // Windows file in linux
            }
        }
        readLock.unlock();

        // Validate FASTQ format
        if (lines[0].empty() || lines[0][0] != '@' ||
            lines[2].empty() || lines[2][0] != '+') {
            throw RNAException(
                    ExceptionType::FORMAT_ERROR,
                    "Invalid FASTQ format"
            );
        }

        read = std::make_shared<Read>();
        read->name = lines[0].substr(1);  // Remove '@'
        read->sequence[0] = lines[1];
        read->sequence[1] = lines[1];
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
                default:
                    c = 'N';
                    break;
            }
        }
        read->length = lines[1].length();
        read->quality = lines[3];

        return true;
    }
    void ReadAligner::processReadFile(std::ifstream& file,std::ofstream& outFile,std::ofstream& alignStatusFile,std::ofstream& alignProgressFile,std::mutex& outputLock,std::mutex& alignStatusLock,std::mutex& alignProgressLock,std::mutex& readFileLock,int& totalReadsProcessed) {




        int64_t totalAlignTime = 0;
        int64_t totalStitchTime = 0;
        // count the number of reads processed in this call
        auto AligningStartTime = std::chrono::high_resolution_clock::now();

        while (loadReadFromFastq(file,readFileLock)) {

            ++readCount;
            auto startTime = std::chrono::high_resolution_clock::now();
            seedMapping->processRead(read);

            auto alignEndTime = std::chrono::high_resolution_clock::now();
            stitchingManagement->processAlignments(seedMapping->aligns,read);
            auto stitchEndTime = std::chrono::high_resolution_clock::now();



            outputLock.lock();
            totalReadsProcessed++;
            if(!partialOutput || totalReadsProcessed < 10000){
                if (stitchingManagement->status == StitchingManagement::SUCCESS){
                    for (int j = 0; j < stitchingManagement->numGoodTranscripts_; ++j) {
                        auto &t = stitchingManagement->goodTranscripts_[j];
                        outFile << SAMEntry(t, *read) << '\n';
                    }
                }
            }
            outputLock.unlock();

            stitchingManagement->clear();



            alignStatusLock.lock();
            if(!partialOutput || totalReadsProcessed < 10000){
                alignStatusFile << read->name << '\t' << "AlignTime:\t"
                                << std::chrono::duration_cast<std::chrono::microseconds>(alignEndTime - startTime).count() << " StitchTime:\t"
                                << std::chrono::duration_cast<std::chrono::microseconds>(stitchEndTime - alignEndTime).count() << '\n';
            }
            alignStatusLock.unlock();


            seedMapping->clear();




            totalAlignTime += std::chrono::duration_cast<std::chrono::microseconds>(alignEndTime - startTime).count();
            totalStitchTime += std::chrono::duration_cast<std::chrono::microseconds>(stitchEndTime - alignEndTime).count();
            alignProgressLock.lock();
            if ((!partialOutput) ||  totalReadsProcessed % 1000000 == 0){
                auto nowTime = std::chrono::high_resolution_clock::now();
                alignProgressFile << "Thread: " << threadId  << "\t";
                alignProgressFile << "Current time:" << std::chrono::duration_cast<std::chrono::seconds>(nowTime - AligningStartTime).count() << " s\t";
                alignProgressFile << "Current reads processed: " << readCount  << '\t';
                alignProgressFile << "Average align time: " << totalAlignTime/readCount << " us\t";
                alignProgressFile << "Average stitch time: " << totalStitchTime/readCount << " us\n";

            }
            alignProgressLock.unlock();
        }



    }


}