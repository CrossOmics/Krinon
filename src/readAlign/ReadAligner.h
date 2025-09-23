#ifndef RNAALIGNREFACTORED_READALIGNER_H
#define RNAALIGNREFACTORED_READALIGNER_H


#include "../genome/GenomeIndex.h"
#include "../utils/types.h"
#include <vector>
#include "SeedMapping.h"
#include "StitchingManagement.h"
#include <memory>
#include <fstream>
#include <sstream>
#include <mutex>
#include <queue>


namespace rna {
    class ReadFile;
    class ReadAligner {
    private:

            const GenomeIndex& genomeIndexPrefix;
            std::unique_ptr<SeedMapping> seedMapping; // seedMapping, get alignments of the read
            std::unique_ptr<StitchingManagement> stitchingManagement;// stitch the previously obtained alignments into transcripts

            ReadPtr read; // now processing Read
            std::vector<TranscriptPtr> bestTranscripts;// the transcripts generated from the read
            std::vector<SJDBOutput> sjdbCandidates; // splice junction candidates from the read



            int threadId{0};
            int readCount{0};
            int inputBufferSize = 30000;
            int nowReadInd;
            int nowQueueSize{0};
            bool queueEmpty{true};
            std::vector<ReadPtr> readBufferQueue;



    public:


        bool partialOutput{false};

        void setConfig(const StitchingConfig& scfg, const StitchingScoreConfig& sscfg, const SeedMappingConfig& smcfg) {
            seedMapping->setConfig(smcfg);
            stitchingManagement->setConfig(scfg,sscfg);
        }
        ReadAligner(const GenomeIndex& gInPre, int tId = 0) : genomeIndexPrefix(gInPre), threadId(tId) {
            seedMapping = std::make_unique<SeedMapping>( genomeIndexPrefix,SeedMappingConfig());
            stitchingManagement = std::make_unique<StitchingManagement>(StitchingConfig(), gInPre);
            inputBufferSize = 30000;
            readBufferQueue.resize(inputBufferSize);
            for (int i = 0; i < inputBufferSize; ++i) {
                readBufferQueue[i] = std::make_shared<Read>();
            }
        }

        bool loadReadFromFastq(ReadFile& file);
        void processReadFile(ReadFile& file,FILE* outFile,std::ofstream& alignStatusFile,std::ofstream& alignProgressFile,std::mutex& outputLock,std::mutex& alignStatusLock,std::mutex& alignProgressLock,int& totalReadsProcessed);
    };
}
#endif //RNAALIGNREFACTORED_READALIGNER_H
