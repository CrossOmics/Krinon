#ifndef RNAALIGNREFACTORED_STITCHINGWINDOW_H
#define RNAALIGNREFACTORED_STITCHINGWINDOW_H

#include <cstdint>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

namespace rna {
    // an align in a window
    class WindowAlign {
    public:
        int64_t readStart{0};
        int64_t genomeStart{0};
        int64_t length{0};
        int64_t score{0};
        bool isAnchor{false};
        int64_t isj{-1};// annotation index
        int iFragment{0};

        bool operator<(const WindowAlign &other) const {
            if (genomeStart != other.genomeStart) return genomeStart < other.genomeStart;
            return length < other.length;
        }
    };

    class Window {
    public:
        int chrIndex{0};
        int direction{0};
        int numAnchors{0};
        int64_t startBin{0};
        int64_t endBin{0};
        WindowAlign *aligns{nullptr};
        int alignNum{0};
        int minLengthWhenFull{0};
        int firstFragmentAlignNum{0};
        int secondFragmentAlignNum{0};
    };
}

#endif //RNAALIGNREFACTORED_STITCHINGWINDOW_H
