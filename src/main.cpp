#include "genome/GenomeIndex.h"
#include "readAlign/ReadAligner.h"
#include "readAlign/ReadAlignMultiThread.h"

int main(int argc, char* argv[]){
    // usage: ./RNAAlignRefactored <reference_genome_file> <read_file>
    // or you can also change their names to test.fa and test.fastq
    // do not use the default genomeIndexPrefixConfig, it is not suitable for large genome
    // kMerSize is the basic k-mer size
    // when a k-mer have more than minExtendRep repetitions, extend it locally, (see GenomeIndexPrefixConfig.cpp)
    // extendLength is the newly added length of the k-mer
    // if the extended k-mer is still repeated more than minExtendRep times, extend it again

    // it is a rough version mainly for testing the index building and aligning process
    // just have made the stitching process work, but probably there are still some bugs and wrong results
    // the output files are:
    // <read_file>.sam: the main result file, in SAM format. The 1,3,4 columns are the read name,
    // the reference name and the position of the read in the reference genome
    // these are the most important colums and are expected to be right
    // <read_file>.alignTime:
    // a timestamp file, each line is the time (us) taken to align and stitch a read
    std::string genomeFile = "test.fa";
    std::string readFile = "test.fastq";
    std::string mode = "both";
    std::string outputMode = "all"; // "all" or "test" in test mode, only output the result and time of first 10000 reads
    //mode = "both";
    mode = "readAlign";
    if (argc == 3){
        genomeFile = argv[1];
        readFile = argv[2];
    }

    if (argc == 5){
        genomeFile = argv[1];
        readFile = argv[2];
        mode = argv[3];
        outputMode = argv[4];
    }
    if(mode == "both" || mode == "genomeGenerate"){
        rna::Genome g(genomeFile);
        rna::GenomeIndexPrefix genomeIndexPrefix(g);
        //These parameter controls the size of the k-mer index

        //Must modify the parameters below !!!
        genomeIndexPrefix.setConfig(rna::GenomeIndexPrefixConfig{
                .kMerSize = 10, // for large genome, suggest using 14-mer
                .extendLength = 2, // for large genome, suggest using 4
                .extendIndexSize = 17, // (1 << (2 * extendLength)) + 1 , suggest using 257
                .minExtendRep = 1000000, // minimum number of repetitions for extension ,
                // suggest using 100 or more to prevent from creating too many extend mer
                .maxLayer = 3, // max number of the depth of the extension, to control the size of the index
                .twoDirections = true
        });
        genomeIndexPrefix.build();
        genomeIndexPrefix.write("test");
    }
    if(mode == "both" || mode == "readAlign"){
        rna::GenomeIndexPrefix genomeIndexPrefix;
        genomeIndexPrefix.load("test");
        /*rna::ReadAligner aligner(genomeIndexPrefix);
        if (outputMode == "test") aligner.partialOutput = true;
        aligner.processReadFile(readFile); // the read file in fastq format*/
        rna::ReadAlignMultiThread readAlignMultiThread;
        readAlignMultiThread.processReadFile(readFile,4,genomeIndexPrefix,outputMode == "test");
    }
    if (mode == "test"){
        //without write & load
        rna::Genome g(genomeFile);
        rna::GenomeIndexPrefix genomeIndexPrefix(g);
        //These parameter controls the size of the k-mer index

        //Must modify the parameters below !!!
        genomeIndexPrefix.setConfig(rna::GenomeIndexPrefixConfig{
                .kMerSize = 10, // for large genome, suggest using 14-mer
                .extendLength = 2, // for large genome, suggest using 4
                .extendIndexSize = 17, // (1 << (2 * extendLength)) + 1 , suggest using 257
                .minExtendRep = 1000000, // minimum number of repetitions for extension ,
                // suggest using 100 or more to prevent from creating too many extend mer
                .maxLayer = 3, // max number of the depth of the extension, to control the size of the index
                .twoDirections = true
        });
        genomeIndexPrefix.build();
        rna::ReadAlignMultiThread readAlignMultiThread;
        readAlignMultiThread.processReadFile(readFile,4,genomeIndexPrefix,outputMode == "test");
    }


};