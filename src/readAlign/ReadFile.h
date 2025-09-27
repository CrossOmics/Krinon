#ifndef RNAALIGNREFACTORED_READFILE_H
#define RNAALIGNREFACTORED_READFILE_H
#include "../utils/types.h"
#include <fstream>
// #include <mutex>
#include <memory>
#include <chrono>

namespace rna {
    class Parameters;

    /**
     * Each aligner thread will have its own instance of this class.
     */
    class ReadFile {
    private:
        char* readFileBuffer{nullptr};
        // char* readFileBuffer2{nullptr};

        /**
         * This is the _expected_ number of bytes that we are allowed to read.
         * Becomes important in the multi-threaded case.
         * A value of `0` means no restrictions.
         */
        int64_t expectedReadableBytes{0};

        /**
         * When chunking the FASTQ input file, it is possible to chunk in the middle of a single
         * read, thus leaving it incomplete such that parsing it will results in an error.
         * To prevent this, we move the initial file pointer until we see an `@` character. This
         * member records how many bytes we moved.
         * 
         * The _real_ maximum allowed reads thus becomes:
         * 
         *      `expectedReadableBytes - initialReadAdjust`
         */
        int64_t initialReadAdjust{0};


        /**
         * Traces how many bytes we have read.
         * We _MUST_ stop reading the moment:
         * 
         *      `(expectedReadableBytes - initialReadAdjust) <= totalBytesRead`
         * 
         * If we go beyond this, it is possible to read another entry that was meant for another
         * thread, getting a duplicate on the output.
         */
        int64_t totalBytesRead{0};

        /**
         * We will process a file from this byte address.
         */
        int64_t beginFrom{0};

        /**
         * This is the unique index of this reader instance.
         * When multi-threaded, this will be set to the thread ID.
         */
        int readerIndex{0};

    public:
        enum ReadType{
            single, paired
        };
        std::string readFileName;
        // std::string readFileName2;

        /**
         * The stream to the underlying FASTQ file.
         * Multiple streams can co-exist, we will align them such that they do not
         * get into any conflict.
         */
        std::ifstream readFile;
        // std::ifstream readFile2; // for paired-end
        // std::mutex fileLock;
        ReadType readType{single};

        //todo move them to a new class
        // int64_t uniqueReadCount{0};
        // int64_t multiReadCount{0};
        // int64_t totalReadCount{0};

        ReadFile(std::string &readType, int64_t start = 0, int64_t expected = 0, int index = 0);
        ReadFile(const Parameters &P, int64_t start = 0, int64_t expected = 0, int index = 0);

        ~ReadFile(){
            delete[] readFileBuffer;
            // delete[] readFileBuffer2;
        };
        void openFiles(const std::string& filename1,const std::string& filename2="");
        bool loadReadFromFastq(ReadPtr& read);
        int loadReadChunkFromFastq(std::vector<ReadPtr>& reads, int chunkSize);
        void closeFiles();
    };
}
#endif //RNAALIGNREFACTORED_READFILE_H
