#include "TimeReport.h"


namespace RefactorProcessing {
    void TimeReport::init(int threadNum, const std::string outputFileName) {
        totalReadsProcessed = 0;
        reportStages = new std::atomic<bool>[threadNum];
        for (int i = 0; i < threadNum; ++i) {
            reportStages[i] = false;
        }
        threadsUnfinished = threadNum;
        outputFile.open(outputFileName);
        previousReportTime = std::chrono::system_clock::now();
        startTime = previousReportTime;
    }


    void TimeReport::tryActivateReport(int threadNum) {
        //thread 0
        auto now = std::chrono::system_clock::now();
        if (std::chrono::duration_cast<std::chrono::seconds>(now - previousReportTime).count() >= 60) {
            // activate report
            previousReportTime = now;
            for (int i = 0; i < threadsUnfinished.load(); ++i) {
                reportStages[i] = true;
            }
            totalReadsProcessed = 0;
            threadsUnfinished = threadNum;

        }
    }

    void TimeReport::tryReportProgress(int threadId, int64_t readsProcessed) {
        if (reportStages[threadId]){
            reportStages[threadId] = false;
            totalReadsProcessed += readsProcessed;
            if(--threadsUnfinished == 0){
                // all threads have reported, output progress
                auto now = std::chrono::system_clock::now();
                int totalTime = std::chrono::duration_cast<std::chrono::seconds>(now - startTime).count();
                int64_t reads = totalReadsProcessed;
                outputFile << "Total reads processed:\t" << reads << "\t Now time:\t"
                           << totalTime << " s\t"
                           << "Now speed: \t"
                           <<std::fixed << std::setprecision(1)
                           << double (reads) * 3600.0 / 1000000*double (totalTime) << "M reads/hour\n";
            }
        }



    }
}
