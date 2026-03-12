#ifndef RNAALIGNREFACTORED_SJDBR_H
#define RNAALIGNREFACTORED_SJDBR_H
// todo
#include <string>
#include <vector>
#include "../io/Parameters.h"
#include "Genome.h"

namespace RefactorProcessing {
    // sj data required for stitching
    struct exon {
        int64_t trID;
        int64_t start;
        int64_t end;
        int64_t geneId;

        bool operator<(const exon &b) const {
            if (trID != b.trID) return trID < b.trID;
            if (start != b.start) return start < b.start;
            return false;
        }
    };

    struct exonTrLoci {
        int64_t trStart;
        int64_t trEnd;
        int64_t trID;
        int64_t exonStart;
        int64_t exonEnd;
        int64_t geneId;

        exonTrLoci() = default;


        bool operator<(const exonTrLoci &b) const {
            if (trID != b.trID) return trID < b.trID;
            if (trStart != b.trStart) return trStart < b.trStart;
            if (trEnd != b.trEnd) return trEnd < b.trEnd;
            if (exonStart != b.exonStart) return exonStart < b.exonStart;
            if (exonEnd != b.exonEnd) return exonEnd < b.exonEnd;
            return false;
        }
    };

    struct sjStride {
        int64_t sjStart;
        int64_t sjEnd;
        int64_t sjStrand; // 0,1,2
        int64_t geneId;
        sjStride() = default;
        sjStride(int64_t sjStart, int64_t sjEnd, int64_t sjStrand, int64_t geneId)
                : sjStart(sjStart), sjEnd(sjEnd), sjStrand(sjStrand), geneId(geneId) {};
        bool operator<(const sjStride &b) const {
            if (sjStart != b.sjStart) return sjStart < b.sjStart;
            if (sjEnd != b.sjEnd) return sjEnd < b.sjEnd;
            if (sjStrand != b.sjStrand) return sjStrand < b.sjStrand;
            return false;
        }

        bool operator==(const sjStride &b) const {
            return sjStart == b.sjStart && sjEnd == b.sjEnd && sjStrand == b.sjStrand;
        }

        bool operator!=(const sjStride &b) const {
            return !(*this == b);
        }

    };

    struct sjdbSortRecord{
        size_t start;
        size_t end;
        size_t originalIndex;
        bool operator<(const sjdbSortRecord &b) const {
            if (start != b.start) return start < b.start;
            if (end != b.end) return end < b.end;
            return false;
        }
    };


    struct sjdbLoci {
        std::string chr;
        int64_t start;
        int64_t end;
        int8_t strand;
        uint8_t priority;
        std::set<size_t> gene;

        sjdbLoci() = default;

        sjdbLoci(std::string chr, int64_t start, int64_t end, int8_t strand, uint8_t priority = 0)
                : chr(std::move(chr)), start(start), end(end), strand(strand), priority(priority) {};

    };

    struct sjHash{
        int32_t hash;
        unsigned int length;
        bool operator<(const sjHash &b) const {
            if (hash != b.hash) return hash < b.hash;
            return length > b.length; // AATT# < AAT##
        }
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

        std::ofstream logFile_;

        size_t sjdbOverhang{100};
        size_t sjdbLength{220};
        int limitSjdbInsertN{1000000};



        std::string sjdbSeq_;

        std::map<std::string,int64_t> chromosomeNameToIndex_;

        std::vector<std::string> transcriptID_;
        std::vector<int> transcriptStrand_;
        std::vector<std::string> geneID_;
        std::vector<std::pair<std::string, std::string>> geneAttr_;
        std::vector<exon> exonLoci_;
        std::vector<sjdbLoci> sjdbLoci_;

    public:
        void setParam(const Parameters& P);

        void loadGTF(const std::string &gtfFile, const Genome &genome);

        int fillSjdbLoci(const std::string &dirOut, Genome &genome);

        friend GenomeIndex;

        friend Genome;
    };
}

#endif //RNAALIGNREFACTORED_SJDBR_H
