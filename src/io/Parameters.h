#ifndef RNAALIGNREFACTORED_PARAMETERSR_H
#define RNAALIGNREFACTORED_PARAMETERSR_H

#include <string>
#include "../utils/argparse.hpp"

namespace RefactorProcessing {
    class Parameters {
        // decode and store parameters from command line or config file
    private:
        argparse::ArgumentParser program;
    public:


        std::string mode;
        int workingThreads;

        //Genome Generate Parameters
        int genomeBinSize;
        std::string genomeFile;

        int kMerSize = 10;
        int extendAlternativeByte = 64;// 16 * 4 byte

        std::string gtfFile;
        bool insertSJ{false};
        int sjdbOverhang{100};
        int sjdbLength{220};
        int limitSjdbInsertN{20000};// the default in STAR is 1000000, 20000 is for smaller test.

        std::string sjdbGTFfeatureExon{"exon"};
        std::string sjdbGTFTagExonParentTranscriptId{"transcript_id"};
        std::string sjdbGTFChrPrefix;
        std::string sjdbGTFTagExonParentGene{"gene_id"};
        std::vector<std::string> sjdbGTFTagExonParentGeneName{"gene_name"};
        std::vector<std::string> sjdbGTFTagExonParentGeneType{"gene_type", "gene_biotype"};

        //Read Align Parameters - basic
        bool isPaired;
        std::string readFile;
        std::string readFile2;

        //Read Align Parameters - seed search

        int minSplitLength{20};
        int maxSeedPerRead{1000};

        //Read Align Parameters - stitching
        int maxAnchorRep{50};
        int winBinSizeLog{16}; //the size of each window is k* 2^winBinSizeLog
        int winAnchorDistBins{9};
        int flankSize{4};
        int maxWindows{10000};
        int maxSeedPerWindows{50};
        int maxRep{10000};
        int maxExons{20}; // max exons number in a transcript
        int transcriptStoredMax{100}; // maximum number of transcripts to store

        int outFilterMultimapMax{10};
        int maxMismatch{10};
        int multimapScoreRange{1};
        double outFilterScoreMinOverLRead{0.66};
        double outFilterMatchMinOverLRead{0.66};

        //Read Align Parameters - score

        int MATCH_SCORE{1};
        int MISMATCH_PENALTY{-1};
        int GAP_OPEN_PENALTY{0};
        int DEL_OPEN_PENALTY{-2};
        int DEL_EXTEND_PENALTY{-2};
        int INS_OPEN_PENALTY{-2};
        int INS_EXTEND_PENALTY{-2};
        int SCORE_STITCH_SJ_SHIFT{1};// maximum score reduction while searching for SJ boundaries
        int SCORE_GAP_GCAG{-4}; // GC/AG and CT/GC junction penalty
        int SCORE_GAP_ATAC{-8}; // AT/AC and GT/AT junction penalty
        int SCORE_GAP_NON_CANONICAL{-8};
        int SCORE_ANNOTATED_SJ{2};
        int MAX_SJ_REPEAT_SEARCH{255};
        int MIN_INTRON_LENGTH{21};
        int MAX_INTRON_LENGTH{2147483647};
        // max mismatch allowed for different SJ types: non-canonical, GT-AG, GC-AG, AT-AC, -1 means no limit
        int MAX_MISMATCH_FOR_SJ[4]{0, -1, 0, 0};


        //Output Parameters
        std::string genomeGenerateFileStoreDir;
        std::string outPutDir;
        std::string outLogFile;
        std::string outErrorFile;
        std::string outWarningFile;
        std::string outProgressFile;


        Parameters();

        int process(int argc, char *argv[]);
    };
}
#endif //RNAALIGNREFACTORED_PARAMETERSR_H
