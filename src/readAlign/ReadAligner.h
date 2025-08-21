#ifndef RNAALIGNREFACTORED_READALIGNER_H
#define RNAALIGNREFACTORED_READALIGNER_H

#include "../genome/GenomeIndex.h"
#include "../genome/GenomeIndexPrefix.h"
#include "../utils/types.h"
#include <vector>
#include "SeedMapping.h"
#include "StitchingManagement.h"
#include <memory>
#include <fstream>
#include <sstream>
#include <mutex>

namespace rna {

    class ReadAligner {
    private:

            const GenomeIndexPrefix& genomeIndexPrefix;
            std::unique_ptr<SeedMapping> seedMapping; // seedMapping, get alignments of the read
            std::unique_ptr<StitchingManagement> stitchingManagement;// stitch the previously obtained alignments into transcripts

            ReadPtr read; // now processing Read
            std::vector<TranscriptPtr> bestTranscripts;// the transcripts generated from the read
            std::vector<SJDBOutput> sjdbCandidates; // splice junction candidates from the read

            int threadId{0};
            int readCount{0};

    public:
        bool partialOutput{false};
        ReadAligner(const GenomeIndexPrefix& gInPre,int tId = 0) : genomeIndexPrefix(gInPre),threadId(tId) {
            seedMapping = std::make_unique<SeedMapping>( genomeIndexPrefix,SeedMappingConfig());
            stitchingManagement = std::make_unique<StitchingManagement>(stitchingConfig(),gInPre);
        }
        bool loadReadFromFastq(std::ifstream& s,std::mutex& readLock);
        void processReadFile(std::ifstream& file,std::ofstream& outFile,std::ofstream& alignStatusFile,std::ofstream& alignProgressFile,std::mutex& outputLock,std::mutex& alignStatusLock,std::mutex& alignProgressLock,std::mutex& readFileLock,int& totalReadsProcessed);
    };
}
#endif //RNAALIGNREFACTORED_READALIGNER_H
