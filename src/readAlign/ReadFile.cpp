#include "ReadFile.h"
#include "../utils/exceptions.h"
#include <sstream>
#include <algorithm>
#include "../utils/Parameters.h"
#include <cassert>
#include <plog/Log.h>

namespace rna {
    ReadFile::ReadFile(std::string &readType, int64_t start, int64_t expected, int index) {
        assert(readType == "single");
        this->readType = single;
        // if (readType == "single") {
        //     this->readType = single;
        // } else if (readType == "paired") {
        //     this->readType = paired;
        // }
        readFileName = "";
        beginFrom = start;
        expectedReadableBytes = expected;
        readerIndex = index;
        // readFileName2 = "";
    }

    ReadFile::ReadFile(const Parameters &P, int64_t start, int64_t expected, int index) {
        assert(!P.isPaired);
        // if (P.isPaired) {
        //     this->readType = paired;
        //     readFileName = P.readFile;
        //     // readFileName2 = P.readFile2;
        // } else {
        //     this->readType = single;
        //     readFileName = P.readFile;
        // }
        this->readType = single;
        readFileName = P.readFile;
        beginFrom = start;
        expectedReadableBytes = expected;
        readerIndex = index;
    }

    void ReadFile::openFiles(const std::string &filename1, const std::string &filename2) {
        // Open the input stream and set the input buffer size
        assert(readType == single);
        // TODO: Tune this!!!
        int64_t inputBufferSize = 30000000; 
        readFileBuffer = new char[inputBufferSize];
        readFile.rdbuf()->pubsetbuf(readFileBuffer, inputBufferSize);
        readFile.open(filename1);

        // Move to where we are told to start from
        PLOG_DEBUG << "[Reader " << readerIndex << "] Will begin reading from address " << beginFrom;
        readFile.seekg(beginFrom);

        // Move up until you see an `@` character
        while (true) {
            if (readFile.peek() == '@')
                break;
            readFile.get();
            ++initialReadAdjust;
        }
        PLOG_DEBUG << "[Reader " << readerIndex << "] Adjusted start address by " << initialReadAdjust;

        assert(readFile.peek() == '@');

        // if (readType == paired) {
        //     readFileBuffer2 = new char[inputBufferSize];
        //     readFile2.rdbuf()->pubsetbuf(readFileBuffer2, inputBufferSize);
        //     readFile2.open(filename2);
        // }
    }

    void ReadFile::closeFiles() {
        assert(readType == single);
        if(readFile.is_open()) readFile.close();
        // if (readType == paired) {
        //     if(readFile2.is_open()) readFile2.close();
        // }
    }

    /**
     * This reads the 4 lines from a FASTQ read instance and puts them into the
     * given `read` shared pointer.
     * UDPATE: This now explicitly assumes that the underlying `readFile` has been
     *         alligned such that it is seeking the `@` character that begins the
     *         read.
     * This will return `false` if there is nothing left to read, or we have hit our
     * maximum number of bytes that we are allowed to read.
     */
    bool ReadFile::loadReadFromFastq(rna::ReadPtr &read) {
        assert(readType == single);

        // First, are we allowed to continue?
        if ((expectedReadableBytes > 0) && ((expectedReadableBytes - initialReadAdjust) <= totalBytesRead))
            return false;

        // This MUST be aligned
        char current = readFile.peek();
        if (current == EOF)
            return false;
        assert(current == '@');

        std::array<std::string, 4> lines1;
        std::array<std::string, 4> lines2;

        for (int i = 0; i < 4; ++i) {
            if (!(std::getline(readFile, lines1[i]))) {
                return false;
            }

            // Adjust how many bytes we read
            totalBytesRead += lines1[i].size();

            if (lines1[i][lines1[i].size() - 1] == '\r') {
                lines1[i].pop_back(); // Windows file in linux
            }
        }
        // if (readType == paired) {
        //     for (int i = 0; i < 4; ++i) {
        //         if (!(std::getline(readFile2, lines2[i]))) {
        //             return false;
        //         }
        //         if (lines2[i][lines2[i].size() - 1] == '\r') {
        //             lines2[i].pop_back(); // Windows file in linux
        //         }
        //     }
        // }

        read = std::make_shared<Read>();
        std::istringstream nameStream(lines1[0].substr(1));
        nameStream >>read->name;
        if(readType == paired){
            std::reverse(lines2[1].begin(),lines2[1].end());
            std::reverse(lines2[3].begin(),lines2[3].end());
            for (char &c: lines2[1]) {
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
            read->sequence[0] = lines1[1] + '#' + lines2[1];
            read->sequence[1] = lines1[1] + '#' + lines2[1];
            read->length = lines1[1].length() + 1 + lines2[1].length();
            read->quality = lines1[3] + ' ' + lines2[3];
            read->mate1Length = lines1[1].length();
            read->mate2Length = lines2[1].length();
        }else {
            read->sequence[0] = lines1[1];
            read->sequence[1] = lines1[1];
            read->length = lines1[1].length();
            read->quality = lines1[3];
        }

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
                case '#':
                    break;
                default:
                    c = 'N';
                    break;
            }
        }
        return true;
    }

    int ReadFile::loadReadChunkFromFastq(std::vector<ReadPtr>& reads, int chunkSize){
        int count = 0;
        // fileLock.lock();
        for(int i = 0; i < chunkSize; ++i){
            if(!loadReadFromFastq(reads[i])){
                // fileLock.unlock();
                return count;
            }
            ++count;
        }
        // fileLock.unlock();
        return count;
    }
}