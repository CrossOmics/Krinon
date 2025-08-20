#include "ReadAligner.h"
#include "../utils/exceptions.h"
#include <fstream>
#include <filesystem>
#include <chrono>
namespace rna {
    bool ReadAligner::loadReadFromFastq(std::ifstream& s) {
        std::array<std::string, 4> lines;

        for (int i = 0; i < 4; ++i) {
            if (!(std::getline(s, lines[i]))) {
                return false;
            }
            if (lines[i][lines[i].size() - 1] == '\r') {
                lines[i].pop_back(); // Windows file in linux
            }
        }

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
    void ReadAligner::processReadFile(const std::string& filename) {
        std::ifstream file(filename);
        std::ofstream outFile(filename+ ".sam");
        std::ofstream alignStatusFile(filename + ".alignTime");
        std::ofstream alignProgressFile(filename + ".progress");
        //std::ofstream sjdbFile(filename + ".sj");
        if (!file) {
            throw RNAException(
                    ExceptionType::FILE_NOT_FOUND,
                    "Cannot open read file: " + filename
            );
        }
        int i = 0;
        int64_t totalAlignTime = 0;
        int64_t totalStitchTime = 0;
        auto AligningStartTime = std::chrono::high_resolution_clock::now();
        while (loadReadFromFastq(file)){
            ++i;

            auto startTime = std::chrono::high_resolution_clock::now();
            seedMapping->processRead(read);

            auto alignEndTime = std::chrono::high_resolution_clock::now();
            stitchingManagement->processAlignments(seedMapping->aligns,read);


            sjdbCandidates.insert(sjdbCandidates.end(),
                                  stitchingManagement->getSJDB().begin(),
                                  stitchingManagement->getSJDB().end());

            if(!partialOutput || i < 10000){


                outFile << SAMEntry(*stitchingManagement->getBestTranscript(), *read) << '\n';
            }

            stitchingManagement->clear();
            auto stitchEndTime = std::chrono::high_resolution_clock::now();

            if(!partialOutput || i < 10000){
                alignStatusFile << read->name << '\t' << "AlignTime:\t"
                                << std::chrono::duration_cast<std::chrono::microseconds>(alignEndTime - startTime).count() << " StitchTime:\t"
                                << std::chrono::duration_cast<std::chrono::microseconds>(stitchEndTime - alignEndTime).count() << '\n';
            }


            seedMapping->clear();

            outFile.flush();
            alignStatusFile.flush();


            totalAlignTime += std::chrono::duration_cast<std::chrono::microseconds>(alignEndTime - startTime).count();
            totalStitchTime += std::chrono::duration_cast<std::chrono::microseconds>(stitchEndTime - alignEndTime).count();
            if ((!partialOutput) ||  i % 1000000 == 0){
                auto nowTime = std::chrono::high_resolution_clock::now();
                alignProgressFile << "Total time:" << std::chrono::duration_cast<std::chrono::seconds>(nowTime - AligningStartTime).count() << " s\t";
                alignProgressFile << "Total reads processed: " << i << '\t';
                alignProgressFile << "Average align time: " << totalAlignTime/i << " us\t";
                alignProgressFile << "Average stitch time: " << totalStitchTime/i << " us\n";

            }
        }


        file.close();
        outFile.close();
    }


}