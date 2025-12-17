#ifndef RNAALIGNREFACTORED_READ_H
#define RNAALIGNREFACTORED_READ_H
#include <string>
namespace RefactorProcessing {
    struct Read {
        std::string sequence[2]; // 0 for forward, 1 for complementary
        std::string quality;
        std::string name;
        int64_t length;
        int64_t mate1Length; // for paired-end reads
        int64_t mate2Length;
    };
}
#endif //RNAALIGNREFACTORED_READ_H
