#ifndef RNAALIGNREFACTORED_READSCANNER_H
#define RNAALIGNREFACTORED_READSCANNER_H

#include <mutex>
#include <memory>
#include <thread>
#include <cstring>
#include <fstream>
#include "Parameters.h"
#include "MemoryMappedInput.hpp"
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
        std::unique_ptr<MemoryMappedInput> mmapFile1_;
        /**
         * TODO: (ARVIN AND YUXUAN)
         *       We should be really really careful with paired-end reads. The current implementation
         *       probably will output gibberish. I will postpone that for now ...
         */
        std::unique_ptr<MemoryMappedInput> mmapFile2_;

    public:
        // We use 2MB read buffers
        const size_t READ_BUFFER_SIZE = 2LLU * 1024 * 1024;

        size_t readFileSize1() {
            return mmapFile1_->size();
        }

        size_t readFileSize2() {
            return mmapFile2_->size();
        }

        bool isEOF(const char* ptr) {
            return ptr >= mmapFile1_->getEndPtr();
        }

        void setParam(const Parameters &P);

        /**
         * Attempt to claim at least `READ_BUFFER_SIZE` amount of input reads
         * to pass to the scanner.
         */
        std::string_view loadFromFastq();

        /**
         * Parse read from the provided buffers and return the number of bytes consumed
         */
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
         * Given inclusive start and exclusive end pointers to a valid memory region that maps
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
        const char* findNextReadStart(const char* startPtr, const char* endPtr);
    };
}
#endif //RNAALIGNREFACTORED_READSCANNER_H
