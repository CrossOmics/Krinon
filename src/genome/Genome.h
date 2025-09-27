#ifndef RNAALIGNMENT_GENOME_H
#define RNAALIGNMENT_GENOME_H

#include <string>
#include <string_view>
#include <vector>
#include <map>
#include <cstdint>
#include "../utils/types.h"

namespace rna{
    struct GenomeSjdbRecord {
        size_t start;
        size_t end;
        uint8_t motif;
        uint8_t shiftLeft;
        uint8_t shiftRight;
        uint8_t strand;
    };

    class Genome{
    private:

        void buildPosToChromosomeMap();

        static constexpr char SPACING_CHAR = '#';               // Character used for spacing in the sequence

    public:
        int binSize = 18;

        struct Chromosome {
            std::string name;
            GenomePos start;
            GenomePos length;
        };
        std::string sequence_;
        size_t originalGenomeLength;
        std::vector<Chromosome> chromosomes_;
        std::map<GenomePos,int64_t> chromosomeStartMap_;
        std::map<std::string,int64_t> chromosomeNameToIndex_;//todo move into class GTF

        int64_t sjdbNum{0};
        std::vector<GenomeSjdbRecord> sjdb;
        std::vector<int64_t> sjDonorStart_;
        std::vector<int64_t> sjAcceptorStart_;

        Genome();
        explicit Genome(const std::string& filename);
        void loadFromFasta(const std::string& filename);
        std::pair<std::string,int64_t> getPosChromosome(int64_t pos) const;
        int64_t getPosChrIndex(int64_t pos) const;
        std::string_view getSequence(size_t pos,size_t length) const;
        void writeChrInfo(const std::string& dirOut) const;
    };
} // namespace rna

#endif //RNAALIGNMENT_GENOME_H
