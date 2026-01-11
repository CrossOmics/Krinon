#ifndef RNAALIGNREFACTORED_STITCHINGR_H
#define RNAALIGNREFACTORED_STITCHINGR_H
//todo
#include <vector>
#include "../genomeRefactored/GenomeIndex.h"

namespace RefactorProcessing {

    class Stitching {
        // stitch alignment fragments
    private:
        // configs

        // data
        std::vector<Align> alignments_;//  alignments to be stitched


    public:
        Stitching(){};
        ~Stitching(){};

        // parameters
        void setParam();

        // stitching process


        void stitchAlignments();

        void createWindows(const std::vector<Align>& alignments);

        void assignAlignmentsToWindows();

        void convertAlignment();// convert alignment to window-alignment

        void assignSingleAlignments();

        void stitchWindow_SingleEnd();

        void stitchWindow_PairedEnd();

    };
};

#endif //RNAALIGNREFACTORED_STITCHINGR_H
