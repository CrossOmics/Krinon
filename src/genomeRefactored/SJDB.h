#ifndef RNAALIGNREFACTORED_SJDBR_H
#define RNAALIGNREFACTORED_SJDBR_H
// todo
#include <string>
#include <vector>
#include "../io/Parameters.h"
namespace RefactorProcessing {
    // sj data required for stitching
    struct sjData {

    };
    class SJDB {
        // splice junction database
        // decode and store essential information of GTF files
    private:
        // configs
        std::string sjdbGTFfeatureExon{"exon"};
        std::string sjdbGTFTagExonParentTranscriptId{"transcript_id"};
        std::string sjdbGTFChrPrefix;
        std::string sjdbGTFTagExonParentGene{"gene_id"};
        std::vector<std::string> sjdbGTFTagExonParentGeneName{"gene_name"};
        std::vector<std::string> sjdbGTFTagExonParentGeneType{"gene_type", "gene_biotype"};

        int sjdbOverhang{100};
        int sjdbLength{220};
        int limitSjdbInsertN{1000000};

    public:
        void setParam();
    };
}

#endif //RNAALIGNREFACTORED_SJDBR_H
