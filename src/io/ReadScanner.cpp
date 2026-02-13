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

        // try to memory-map the files
        try {
            if (std::filesystem::exists(readFileName1_)) {
                fileSize1_ = static_cast<size_t>(std::filesystem::file_size(readFileName1_));
                if (fileSize1_ > 0) {
                    mmapFile1_.reset(new MemoryMappedFile(readFileName1_, (long long)fileSize1_));
                    mmapPos1_ = 0;
                }
            }
            if (isPairedEnd_ && std::filesystem::exists(readFileName2_)) {
                fileSize2_ = static_cast<size_t>(std::filesystem::file_size(readFileName2_));
                if (fileSize2_ > 0) {
                    mmapFile2_.reset(new MemoryMappedFile(readFileName2_, (long long)fileSize2_));
                    mmapPos2_ = 0;
                }
            }
        } catch (const std::exception &e) {
            rna::ErrorRecord().reportError(std::string("Memory map failed: ") + e.what());
        }
    }


    void ReadScanner::loadFromFastq(char* targetBuffer1,char* targetBuffer2,const int bufferSize) {
        if (bufferSize <= 1) return;
        std::lock_guard<std::mutex> lock(readFileLock_);

        auto fillFromMmap = [&](std::unique_ptr<MemoryMappedFile> &mmapFile, size_t &mmapPos, size_t fileSize, char* target) {


            if (!mmapFile || fileSize == 0) {
                if (target) target[0] = '\0';
                return;
            }
            char* base = mmapFile->getMapPtr();
            size_t remaining = fileSize - mmapPos;
            if (remaining == 0) { // reached EOF
                if (target) target[0] = '\0';
                return;
            }
            size_t desired = std::min<size_t>(static_cast<size_t>(bufferSize - 1), remaining);

            // try to include a full read
            if (mmapPos + desired < fileSize ) {
                int lineCount = 0; // must read all 4 lines for each read
                // try to find a new read forward within remaining buffer capacity
                size_t forwardLimit = std::min<size_t>(remaining, static_cast<size_t>(bufferSize - 1));
                size_t i = desired;
                bool foundNewRead = false;
                for (; i < forwardLimit; ++i) {
                    if (base[mmapPos + i] == '@') foundNewRead = true;
                    if (foundNewRead && base[mmapPos + i] == '\n') {
                        ++lineCount;
                        if(lineCount %4==0){
                            ++i; break;
                        }
                    }
                }
                if (i < forwardLimit) {
                    desired = i; // include new read
                } else {
                    // couldn't find new read forward within buffer - try to move back to last new read
                    size_t j = desired;
                    while (j > 0 && base[mmapPos + j - 1] != '@') --j;
                    --j; // move back to previous line
                    if (j > 0) desired = j; // cut
                }
            }

            std::memcpy(target, base + mmapPos, desired);
            if(target) target[desired] = '\0';
            mmapPos += desired;
        };

        // prefer mmap path if available
        if (mmapFile1_) {
            fillFromMmap(mmapFile1_, mmapPos1_, fileSize1_, targetBuffer1);
        }

        if (isPairedEnd_) {
            if (mmapFile2_) {
                fillFromMmap(mmapFile2_, mmapPos2_, fileSize2_, targetBuffer2);
            }
        }
    }

    void ReadScanner::parseRead(Read &r,const char* readBuffer1,const char* readBuffer2,int& r1Length,int& r2Length){
        // parse read from the provided buffers
        // return the number of bytes consumed


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
        r1Length = ptr1 - readBuffer1;

        if (isPairedEnd_){
            // parse read mate
            r.mate1Length = r.sequence[0].length();
            const char* ptr2 = readBuffer2;
            len = std::strcspn(ptr2, "\n"); // no need to parse name again
            ptr2 += len + 1;
            len = std::strcspn(ptr2, "\n");
            r.sequence[0] += '#' + std::string(ptr2, len);
            if (ptr2[len - 1] == '\r') r.sequence[0].pop_back(); // Windows file in linux
            ptr2 += len + 1;
            len = std::strcspn(ptr2, "\n");
            ptr2 += len + 1;
            len = std::strcspn(ptr2, "\n");
            r.quality += ' ' + std::string(ptr2, len);
            if (ptr2[len - 1] == '\r') r.quality.pop_back(); // Windows file in linux
            ptr2 += len + 1;
            r2Length = ptr2 - readBuffer2;
            r.length = r.sequence[0].length();
            r.mate2Length = r.length - r.mate1Length - 1;
        }else {
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



    }

}