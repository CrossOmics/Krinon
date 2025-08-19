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

namespace rna {

    class ReadAligner {
    private:

            GenomeIndexPrefix& genomeIndexPrefix;
            std::unique_ptr<SeedMapping> seedMapping; // seedMapping, get alignments of the read
            std::unique_ptr<StitchingManagement> stitchingManagement;// stitch the previously obtained alignments into transcripts

            ReadPtr read; // now processing Read
            std::vector<TranscriptPtr> bestTranscripts;// the transcripts generated from the read
            std::vector<SJDBOutput> sjdbCandidates; // splice junction candidates from the read


    public:
        bool partialOutput{false};
        ReadAligner(GenomeIndexPrefix& gInPre) : genomeIndexPrefix(gInPre) {
            seedMapping = std::make_unique<SeedMapping>( genomeIndexPrefix,SeedMappingConfig());
            stitchingManagement = std::make_unique<StitchingManagement>(stitchingConfig(),gInPre);
        }
        bool loadReadFromFastq(std::ifstream& s);
        void processReadFile(const std::string& filename);
    };
}
#endif //RNAALIGNREFACTORED_READALIGNER_H
