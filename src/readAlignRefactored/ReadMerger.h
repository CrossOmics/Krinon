#ifndef RNAALIGNREFACTORED_READMERGER_H
#define RNAALIGNREFACTORED_READMERGER_H
#include "Read.h"
#include "Transcript.h"
#include "../io/Parameters.h"
namespace RefactorProcessing {
    class ReadMerger {
    private:
        // configs
        bool isPaired_;
        int maxMismatch_;
        int minOverlapLength_{10};

        // record the merging status of the read
        int mergedLength; // length of the merged part, 0 for unmerged reads , minus for read1 ahead of read0
    public:
        void setParam(const RefactorProcessing::Parameters &P);

        // merge paired-end reads into a single read for alignment, return the merged read (if mergeable) or the original read (if not mergeable)
        std::pair<bool,Read> merge(const Read& r);

        // split the transcript of the merged single read back into paired-end, return the split transcripts
        Transcript splitTranscript(const Transcript& t);
    };
}
#endif //RNAALIGNREFACTORED_READMERGER_H
