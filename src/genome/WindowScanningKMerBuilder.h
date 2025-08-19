#ifndef RNAALIGNREFACTORED_WINDOWSCANNINGKMERBUILDER_H
#define RNAALIGNREFACTORED_WINDOWSCANNINGKMERBUILDER_H
#include "KMerBuilder.h"
#include "../utils/types.h"
namespace rna{
    class WindowScanningKMerBuilder : public KMerIndexBuilder {
    private:
        static constexpr int32_t charToIndex(char c) {
            switch(c) {
                case 'A': return 0;
                case 'C': return 1;
                case 'G': return 2;
                case 'T': return 3;
                default: return -1; // non-ACGT characters are ignored
            }
        }
        std::vector<std::vector<int8_t>> appearance_flag; // records the appearance of each k-mer in the genome
        void fillKMerIndex();
        void outputDebugData() const;
        void scanGenome();
        void buildKMerIndex();
        inline void setFlag(int32_t hash,int32_t length);
        inline bool getFlag(int32_t hash,int32_t length);
        std::vector<int32_t> nowWindowHash;

        int64_t max_14_mer_hash{-1};
        int64_t max_14_mer_match{0};
        int64_t matched_14_mer_num{0};


    public:
        void build(GenomeIndex* index) override;
    };
}
#endif //RNAALIGNREFACTORED_WINDOWSCANNINGKMERBUILDER_H
