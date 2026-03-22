#ifndef RNAALIGNREFACTORED_READSCANNER_H
#define RNAALIGNREFACTORED_READSCANNER_H

#include <mutex>
#include <memory>
#include <thread>
#include <cstring>
#include <fstream>
#include "Parameters.h"
#include "MemoryMappedFile.hpp"
#include "../readAlignRefactored/Read.h"

namespace RefactorProcessing {
    // scan reads, handle multi-threading
    class ReadScanner {
    private:
        // configs
        bool isPairedEnd_;

        std::string readFileName1_;
        std::string readFileName2_;

        // memory-mapped fast path
        std::unique_ptr<MemoryMappedFile> mmapFile1_;
        /**
         * TODO: (ARVIN AND YUXUAN)
         *       We should be really really careful with paired-end reads. The current implementation
         *       probably will output gibberish. I will postpone that for now ...
         */
        std::unique_ptr<MemoryMappedFile> mmapFile2_;

    public:
        // We use 32 GB chunks for now
        const size_t MAX_READ_CHUNK_SIZE = 32LLU * 1024 * 1024 * 1024;
        /**
         * In order to not bother with aligning the end of the partitions, we
         * map a slightly bigger value than `MAX_READ_CHUNK_SIZE`. This tolerance
         * should be enough to account for the longest possible read. We settle 
         * on a 1 KB tolerance for now.
         */
        const size_t CHUNK_END_OVERHANG_ = 1LLU * 1024;
        // const size_t CHUNK_END_OVERHANG_ = 1LLU;
        
        // We use 2MB read buffers
        const size_t READ_BUFFER_SIZE = 2LLU * 1024 * 1024;
        // const size_t READ_BUFFER_SIZE = 2LLU * 128;


        void setParam(const Parameters &P);

        // void loadFromFastq(char* targetBuffer1, char* targetBuffer2, const size_t bufferSize, size_t& buffer1Length, size_t& buffer2Length);
        /**
         * Attempt to claim at least `READ_BUFFER_SIZE` amount of input reads
         * to pass to the scanner.
         */
        std::string_view loadFromFastq();

        /**
         * Parse read from the provided buffers and return the number of bytes consumed
         */
        // void parseRead(Read &r, const char* readBuffer1, const char* readBuffer2, size_t& r1Length, size_t& r2Length);
        size_t parseRead(Read &r, const char* readBuffer1, const char* readBuffer2);

        /**
         * Advances to the end of the current line, calculates its length (ignoring Windows \r),
         * and skips any subsequent newlines/carriage returns to land on the next actual data character.
         */
        inline const char* getNextLine(const char* current, const char* endPtr, size_t& outLength) {
            const char* nl = static_cast<const char*>(std::memchr(current, '\n', endPtr - current));
            if (!nl) {
                outLength = endPtr - current;
                return endPtr; 
            }
            outLength = nl - current;
            
            // Step back one character if this is a Windows CRLF (\r\n) line ending
            if (outLength > 0 && *(nl - 1) == '\r') {
                outLength--;
            }

            const char* next_data = nl + 1;
            while (next_data < endPtr && (*next_data == '\n' || *next_data == '\r')) {
                next_data++;
            }

            return next_data;
        }

        /**
         * Given inclusive start and end pointers to a valid memory region that maps
         * into a valid FASTQ file, it iterates from `startPtr` forward until it reaches
         * a valid FASTQ header.
         * If we find none and we hit `endPtr`, we will just return that.
         * 
         * This function is only called once, for each buffer 
         *
         * Notes
         * -----
         * 
         * - We iterate _forward_ and only forward, since not doing so means that a read
         * that begins slightly before `startPtr` might be processed twice in a multi-threaded
         * environment.
         * - We can't just look for the next `@` character, as it is a valid character that can
         * appear in a quality score.
         */
        std::pair<const char*, size_t> findNextReadStart(const char* startPtr, size_t length) {
            const char* ptr = startPtr;
            const char* endPtr = startPtr + length;

            while (ptr < endPtr) {
                // Skip new lines or windows CRLF
                const char* nl = static_cast<const char*>(std::memchr(ptr, '\n', endPtr - ptr));
                if (!nl) return std::make_pair(endPtr, 0);
                const char* candidate = nl + 1;
                while (candidate < endPtr && (*candidate == '\n' || *candidate == '\r')) {
                    candidate++;
                }
                if (candidate >= endPtr) return std::make_pair(endPtr, 0);

                // If this is `@`, it might be header start, or a quality score.
                if (*candidate == '@') {
                    size_t len_l1, len_l2, len_l3, len_l4;
                    
                    const char* l2 = getNextLine(candidate, endPtr, len_l1);
                    if (l2 >= endPtr) return std::make_pair(endPtr, 0);
                    const char* l3 = getNextLine(l2, endPtr, len_l2);
                    if (l3 >= endPtr) return std::make_pair(endPtr, 0);
                    if (*l3 == '+') {
                        const char* l4 = getNextLine(l3, endPtr, len_l3);
                        if (l4 >= endPtr) return std::make_pair(endPtr, 0);
                        getNextLine(l4, endPtr, len_l4); 
                        if (len_l2 == len_l4) {
                            // Safe to start from this candidate, as it is a legit read.
                            return std::make_pair(candidate, endPtr - candidate);
                        }
                    }
                }
                
                // If we get here, the '@' was a false positive (likely a quality score).
                // Set our search pointer to the candidate '@'. The next iteration will 
                // use `memchr` to instantly jump to the end of this fake header line.
                ptr = candidate; 
            }

            return std::make_pair(endPtr, 0);
        }
    };
}
#endif //RNAALIGNREFACTORED_READSCANNER_H
