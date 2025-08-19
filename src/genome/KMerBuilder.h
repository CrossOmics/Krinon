#ifndef RNAALIGNREFACTORED_KMERBUILDER_H
#define RNAALIGNREFACTORED_KMERBUILDER_H
#include <memory>
#include <string_view>
namespace rna{
    class GenomeIndex;
    class KMerIndexBuilder {
    public:
        virtual void build(GenomeIndex* index) = 0;

        GenomeIndex* genomeIndex;

        static constexpr int32_t charToIndex(char c) {
            switch(c) {
                case 'A': return 0;
                case 'C': return 1;
                case 'G': return 2;
                case 'T': return 3;
                default: return -1;
            }
        }

        static int64_t encodeKMer(const std::string_view & pattern,const int mer_length) {
            if (pattern.length() < mer_length) return -1;
            uint32_t code = 0;
            for (char c : pattern.substr(0, mer_length)) {
                int32_t i = charToIndex(c);
                if(i == -1) return -1; // 非ACGT字符
                code = (code << 2) | i;
            }
            return code;
        }
    };

    class SuffixAutomatonKMerBuilder;
    class WindowScanningKMerBuilder;
}
#endif //RNAALIGNREFACTORED_KMERBUILDER_H
