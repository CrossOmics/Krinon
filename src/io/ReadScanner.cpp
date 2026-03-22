#include "ReadScanner.h"
#include "../log/ErrorRecord.h"
#include <filesystem>
#include <cstring>
#include <algorithm>
#include <sstream>

namespace RefactorProcessing{
    void ReadScanner::setParam(const RefactorProcessing::Parameters &P) {
        isPairedEnd_ = P.isPaired;
        readFileName1_ = P.readFile;
        readFileName2_ = P.readFile2;

        //set stream (keep existing streams for fallback)
        //set stream
        /*if (!readFile1_.is_open()) {
            rna::ErrorRecord().reportError("Cannot open read file: " + readFileName1_);
        }
        if (isPairedEnd_) {
            readFile2_.open(readFileName2_);
            if (!readFile2_.is_open()) {
                rna::ErrorRecord().reportError("Cannot open read file: " + readFileName2_);
            }
        }*/

        mmapFile1_.reset(new MemoryMappedFile(readFileName1_, MAX_READ_CHUNK_SIZE + CHUNK_END_OVERHANG_));
        if (isPairedEnd_) {
            mmapFile2_.reset(new MemoryMappedFile(readFileName2_, MAX_READ_CHUNK_SIZE + CHUNK_END_OVERHANG_));
        }
    }

    std::string_view ReadScanner::loadFromFastq() {
        // Claim `READ_BUFFER_SIZE` amount of data to process
        auto claimed = mmapFile1_->claim(READ_BUFFER_SIZE);
        size_t length = std::min(claimed.second + READ_BUFFER_SIZE, MAX_READ_CHUNK_SIZE);
        // Locate a legit read header and get a view
        auto block = findNextReadStart(claimed.first, length);
        return std::string_view(block.first, block.second);
    }

    // size_t ReadScanner::parseRead(Read &r, const char* readBuffer1, const char* readBuffer2, size_t& r1Length, size_t& r2Length){
    //     const char* ptr1 = readBuffer1;
    //     int len = std::strcspn(ptr1, " ");
    //     r.name = std::string(ptr1, len);
    //     ptr1 += len + 1;
    //     len = std::strcspn(ptr1, "\n");
    //     ptr1 += len + 1;
    //     len = std::strcspn(ptr1, "\n");
    //     r.sequence[0] = std::string(ptr1, len);
    //     if (ptr1[len - 1] == '\r') r.sequence[0].pop_back(); // Windows file in linux
    //     ptr1 += len + 1;
    //     len = std::strcspn(ptr1, "\n");
    //     ptr1 += len + 1;
    //     len = std::strcspn(ptr1, "\n");
    //     r.quality = std::string(ptr1, len);
    //     if (ptr1[len - 1] == '\r') r.quality.pop_back(); // Windows file in linux
    //     ptr1 += len + 1;
    //     r1Length = ptr1 - readBuffer1;

    //     if (isPairedEnd_){
    //         std::cerr << "We are NOT prepared for paired-end reads." << std::endl;
    //         std::abort()
    //         // // parse read mate
    //         // r.mate1Length = r.sequence[0].length();
    //         // const char* ptr2 = readBuffer2;
    //         // len = std::strcspn(ptr2, "\n"); // no need to parse name again
    //         // ptr2 += len + 1;
    //         // len = std::strcspn(ptr2, "\n");
    //         // r.sequence[0] += '#' + std::string(ptr2, len);
    //         // if (ptr2[len - 1] == '\r') r.sequence[0].pop_back(); // Windows file in linux
    //         // ptr2 += len + 1;
    //         // len = std::strcspn(ptr2, "\n");
    //         // ptr2 += len + 1;
    //         // len = std::strcspn(ptr2, "\n");
    //         // r.quality += ' ' + std::string(ptr2, len);
    //         // if (ptr2[len - 1] == '\r') r.quality.pop_back(); // Windows file in linux
    //         // ptr2 += len + 1;
    //         // r2Length = ptr2 - readBuffer2;
    //         // r.length = r.sequence[0].length();
    //         // r.mate2Length = r.length - r.mate1Length - 1;
    //     } else {
    //         r.length = r.sequence[0].length();
    //     }
    //     // generate reverse complement
    //     r.sequence[1].resize(r.length);

    //     for (int i = 0; i < r.length; ++i) {
    //         char& c = r.sequence[0][r.length - 1 - i];
    //         switch (c) {
    //             case 'A':
    //                 r.sequence[1][i] = 'T';
    //                 break;
    //             case 'T':
    //                 r.sequence[1][i] = 'A';
    //                 break;
    //             case 'C':
    //                 r.sequence[1][i] = 'G';
    //                 break;
    //             case 'G':
    //                 r.sequence[1][i] = 'C';
    //                 break;
    //             case '#':
    //                 r.sequence[1][i] = '#';
    //                 break;
    //             case 'N':
    //                 r.sequence[1][i] = 'N';
    //                 break;
    //             default:
    //                 // invalid, treat as 'N'
    //                 r.sequence[1][i] = 'N';
    //                 c = 'N';
    //                 break;
    //         }
    //     }
    // }
    size_t ReadScanner::parseRead(Read &r, const char* readBuffer1, const char* readBuffer2){
        const char* ptr1 = readBuffer1;
        int len = std::strcspn(ptr1, " ");
        r.name = std::string(ptr1, len);
        ptr1 += len + 1;
        len = std::strcspn(ptr1, "\n");
        ptr1 += len + 1;
        len = std::strcspn(ptr1, "\n");
        r.sequence[0] = std::string(ptr1, len);
        if (ptr1[len - 1] == '\r') r.sequence[0].pop_back(); // Windows file in linux
        ptr1 += len + 1;
        len = std::strcspn(ptr1, "\n");
        ptr1 += len + 1;
        len = std::strcspn(ptr1, "\n");
        r.quality = std::string(ptr1, len);
        if (ptr1[len - 1] == '\r') r.quality.pop_back(); // Windows file in linux
        ptr1 += len + 1;
        size_t length = ptr1 - readBuffer1;

        if (isPairedEnd_){
            std::cerr << "We are NOT prepared for paired-end reads." << std::endl;
            std::abort();
        } else {
            r.length = r.sequence[0].length();
        }
        // generate reverse complement
        r.sequence[1].resize(r.length);

        for (int i = 0; i < r.length; ++i) {
            char& c = r.sequence[0][r.length - 1 - i];
            switch (c) {
                case 'A':
                    r.sequence[1][i] = 'T';
                    break;
                case 'T':
                    r.sequence[1][i] = 'A';
                    break;
                case 'C':
                    r.sequence[1][i] = 'G';
                    break;
                case 'G':
                    r.sequence[1][i] = 'C';
                    break;
                case '#':
                    r.sequence[1][i] = '#';
                    break;
                case 'N':
                    r.sequence[1][i] = 'N';
                    break;
                default:
                    // invalid, treat as 'N'
                    r.sequence[1][i] = 'N';
                    c = 'N';
                    break;
            }
        }
        return length;
    }
}