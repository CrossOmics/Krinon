#include <cstring>
#include <filesystem>
#include "OutputSAM.h"
#include "../log/ErrorRecord.h"

namespace RefactorProcessing{
    void OutputSAM::consumerThread(std::atomic<int>& runningProducers) {
        /**
         * TODO: (ARVIN) Put a global atomic flag to check if the user has interrupted 
         *       us or not so we can exit gracefully.
         */
        size_t total = 0;
        while (true) {
            // Grab one prepared buffer (if there is any!)
            BufferedReads* reads = queue_.pop();
            if (reads) {
                // Write the encapsulated string ...
                size_t toWrite = reads->data.length();
                const char* ptr = reads->data.data();
                while (toWrite > 0) {
                    ssize_t written = write(outputFD_, ptr, toWrite);
                    if (written > 0) {
                        toWrite -= written;
                        ptr += written;
                        total += written;
                    }
                }
                // Clear and return the buffer so that another thread can pick it up
                reads->data.clear();
                bufferPool_.push(reads);
            }
            else {
                // The queue is empty. Either the aligners are done, or we have to wait a bit more ...
                if (runningProducers.load(std::memory_order_acquire) == 0) {
                    break;
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
        // if (ftruncate(outputFD_, total) == -1) {
        //     throw std::runtime_error("Unable to truncate the output SAM file.");
        // }
    }
}
