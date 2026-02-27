#ifndef RNAALIGNREFACTORED_TIMEREPORT_H
#define RNAALIGNREFACTORED_TIMEREPORT_H
#include <atomic>
#include <chrono>
#include <fstream>
#include <iomanip>
namespace RefactorProcessing {

    class TimeReport {
    public:
        std::atomic<int64_t> totalReadsProcessed;
        std::atomic<bool>* reportStages;
        std::atomic<int> threadsUnfinished;
        std::ofstream outputFile;
        std::chrono::system_clock::time_point startTime;
        std::chrono::system_clock::time_point previousReportTime;

        ~TimeReport() {
            delete[] reportStages;
        }

        void init(int threadNum,const std::string outputFileName);

        void tryActivateReport(int threadNum);

        void tryReportProgress(int threadId,int64_t readsProcessed);

    };
}
#endif //RNAALIGNREFACTORED_TIMEREPORT_H
