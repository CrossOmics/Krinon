#ifndef RNAALIGNMENT_GENOME_H
#define RNAALIGNMENT_GENOME_H
#include <string>
#include <string_view>
#include <vector>
#include <map>
#include <cstdint>
#include "../utils/types.h"
namespace rna{
    class Genome{
    private:

        void buildPosToChromosomeMap();

        static constexpr size_t INDEX_SHIFT = 18;
        static constexpr char SPACING_CHAR = '#';   // Character used for spacing in the sequence

    public:
        struct Chromosome {
            std::string name;
            GenomePos start;
            GenomePos length;
        };
        std::string sequence_;
        std::vector<Chromosome> chromosomes_;
        std::map<GenomePos,int64_t> chromosomeMap_;
        Genome();
        explicit Genome(const std::string& filename);
        void loadFromFasta(const std::string& filename);
        std::pair<std::string,int64_t> getPosChromosome(int64_t pos) const;
        int64_t getPosChrIndex(int64_t pos) const;
        std::string_view getSequence(size_t pos,size_t length) const;
        void writeChrInfo(const std::string& dirOut) const;
    };
}

#endif //RNAALIGNMENT_GENOME_H
