
#include "readAlign/ReadAligner.h"
#include "readAlign/ReadAlignMultiThread.h"
#include "utils/Parameters.h"


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
    std::string outputMode = "all"; // "all" or "test" in test mode, only output the result and time of first 10000 reads

    rna::Parameters P;
    P.process(argc,argv);
    std::string mode = P.mode;
    if (mode.empty()) {
        std::cout << "Please specify the mode: GenomeGenerate, ReadAlign, both or test" << std::endl;

        std::cout << P.program;
        exit(0);
    }
    if (mode == "GenomeGenerate" || mode == "both") {
        rna::Genome g;
        g.binSize = P.genomeBinSize;
        g.loadFromFasta(P.genomeFile);
        rna::GenomeIndex genomeIndex(g);
        genomeIndex.setConfig(P.genomeIndexConfig);
        genomeIndex.build();

        if (!P.gtfFile.empty()) {
            rna::GTF gtf(P.outLogFile);
            gtf.loadGTF(P.gtfFile, *genomeIndex.genome);
            gtf.fillSjdbLoci("./", *genomeIndex.genome);
            gtf.insertJunctions(*genomeIndex.genome, genomeIndex);
        }
        genomeIndex.write(P.genomeGenerateFileStoreDir);
    }
    if (mode == "ReadAlign" || mode == "both") {
        rna::GenomeIndex genomeIndex;
        genomeIndex.setConfig(P.genomeIndexConfig);
        genomeIndex.load(P.genomeGenerateFileStoreDir);

        rna::ReadAlignMultiThread readAlignMultiThread(P);
        readAlignMultiThread.processReadFile( P.threads, genomeIndex, false);
    }


};