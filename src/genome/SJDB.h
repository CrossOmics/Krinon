#ifndef RNAALIGNREFACTORED_SJDB_H
#define RNAALIGNREFACTORED_SJDB_H

#include <string>
#include <utility>
#include <vector>
#include <set>


namespace rna {
    class Genome;
    class GenomeIndex;




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
        int32_t length;
        bool operator<(const sjHash &b) const {
            if (hash != b.hash) return hash < b.hash;
            return length > b.length; // AATT# < AAT##
        }
    };

    struct GTFConfig {
        // tags for processing GTF file
        std::string sjdbGTFfeatureExon{"exon"};
        std::string sjdbGTFTagExonParentTranscriptId{"transcript_id"};
        std::string sjdbGTFChrPrefix;
        std::string sjdbGTFTagExonParentGene{"gene_id"};
        std::vector<std::string> sjdbGTFTagExonParentGeneName{"gene_name"};
        std::vector<std::string> sjdbGTFTagExonParentGeneType{"gene_type", "gene_biotype"};


        int sjdbOverhang{100};
        int sjdbLength{220}; // = sjdbOverhang*2+SJDB_PADDING_LENGTH
        int limitSjdbInsertN{1000000};
    };



    class GTF {
    private:


        std::ofstream &logFile_;
    public:
        GTFConfig config_;
        std::vector<std::string> transcriptID_;
        std::vector<int> transcriptStrand_;
        std::vector<std::string> geneID_;
        std::vector<std::pair<std::string, std::string>> geneAttr_; // geneName, geneType

        std::vector<exon> exonLoci_;
        std::vector<sjdbLoci> sjdbLoci_;

        std::string sjdbSeq_;
        int64_t sjdbSeqLengthSingleStrand_;

        GTF(std::ofstream &outLogFile) : logFile_(outLogFile) {};

        void loadGTF(const std::string &gtfFile, Genome &genome);

        int fillSjdbLoci(const std::string& dirOut, Genome &genome);

        void insertJunctions(Genome &genome, GenomeIndex &genomeIndex);

        void modifyIndex(Genome &genome, GenomeIndex& genomeIndex);

        void clear();

        inline constexpr int charToIndex(char c) {
            switch (c) {
                case 'A': return 0;
                case 'C': return 1;
                case 'G': return 2;
                case 'T': return 3;
                default: return -1; // For 'N' or any other character
            }
        }

        inline sjHash encodeKMer(const std::string_view &pattern, const int mer_length){
            if (pattern.length() < mer_length) {
                return {-1, 0};
            }
            int32_t code = 0;
            int i = 0;
            for (char c: pattern.substr(0, mer_length)) {
                if (charToIndex(c) < 0) {
                    if (i == 0) return {-1, 0};
                    return {((code + 1) << (2 * (mer_length - i))) - 1 , i};
                }
                ++i;
                code = (code << 2) | charToIndex(c);
            }
            return {code , mer_length};
        }

    };
}


#endif //RNAALIGNREFACTORED_SJDB_H
