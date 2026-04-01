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

        mmapFile1_.reset(new MemoryMappedInput(readFileName1_));
        if (isPairedEnd_) {
            mmapFile2_.reset(new MemoryMappedInput(readFileName2_));
        }
    }

    const char* ReadScanner::findNextReadStart(const char* startPtr, const char* endPtr) {
        const char* ptr = startPtr;

        while (ptr < endPtr) {
            while (ptr < endPtr && (*ptr == '\n' || *ptr == '\r')) {
                ptr++;
            }
            if (ptr >= endPtr) return endPtr;

            // If this is `@`, it might be header start, or a quality score.
            if (*ptr == '@') {
                size_t len_l1, len_l2, len_l3, len_l4;
                
                const char* l2 = getNextLine(ptr, endPtr, len_l1);
                if (l2 >= endPtr) return endPtr;
                const char* l3 = getNextLine(l2, endPtr, len_l2);
                if (l3 >= endPtr) return endPtr;
                if (*l3 == '+') {
                    const char* l4 = getNextLine(l3, endPtr, len_l3);
                    if (l4 >= endPtr) return endPtr;
                    getNextLine(l4, endPtr, len_l4); 
                    // Safe to start from this candidate, as it is a legit read.
                    if (len_l2 == len_l4) return ptr;
                }
            }
            
            // Not a valid header. Skip to the start of the next line.
            const char* nl = static_cast<const char*>(std::memchr(ptr, '\n', endPtr - ptr));
            if (!nl) return endPtr;
            ptr = nl + 1;
        }

        // If we are here, it means that 

        return endPtr;
    }

    std::string_view ReadScanner::loadFromFastq() {
        // Claim `READ_BUFFER_SIZE` amount of data to process
        auto endPtr = mmapFile1_->getEndPtr();
        const char* claimed = mmapFile1_->claim(READ_BUFFER_SIZE);
        // Locate the next legit read header and get a view
        const char* next = findNextReadStart(claimed, endPtr);
        // If it is beyond the current buffer, return empty view
        if (next >= claimed + READ_BUFFER_SIZE) return std::string_view();
        return std::string_view(next, endPtr - next);
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
    size_t ReadScanner::parseRead(Read &r, const char* readBuffer1, const char* readBuffer2) {
        const char* ptr1 = readBuffer1;
        //if (*ptr1 == 0) return 0; // invalid input
        int len = std::strcspn(ptr1, " ");
        //if (len == 0) return 0; // invalid read ("\000\000\000\000"), is it expected to happen?
        r.name = std::string(ptr1 + 1 , len - 1 ); // skip '@'

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
    // size_t ReadScanner::parseRead(Read &r, const char* readBuffer1, const char* readBuffer2, size_t r1Length, size_t r2Length) {
    //     size_t lineLen;
    //     const char* currentStart = readBuffer1;
    //     const char* nextStart;
    //     const char* end1 = readBuffer1 + r1Length;
        
    //     // Line 1
    //     r.name = std::string(currentStart, std::strcspn(currentStart, " "));
    //     nextStart = getNextLine(currentStart, end1, lineLen);
    //     currentStart = nextStart;
    //     // Line 2
    //     nextStart = getNextLine(currentStart, end1, lineLen);
    //     r.sequence[0] = std::string(currentStart, lineLen);
    //     currentStart = nextStart;
    //     // Skip line 3
    //     nextStart = getNextLine(currentStart, end1, lineLen);
    //     currentStart = nextStart;
    //     // Line 4
    //     nextStart = getNextLine(currentStart, end1, lineLen);
    //     r.quality = std::string(currentStart, lineLen);
    //     currentStart = nextStart;
    //     size_t length = currentStart - readBuffer1 + lineLen;

    //     if (isPairedEnd_){
    //         std::cerr << "We are NOT prepared for paired-end reads." << std::endl;
    //         std::abort();
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
    //     return length;
    // }
}