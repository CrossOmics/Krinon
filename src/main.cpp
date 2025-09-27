
#include "readAlign/ReadAligner.h"
#include "readAlign/ReadAlignMultiThread.h"
#include "utils/Parameters.h"

#include <plog/Log.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Initializers/ConsoleInitializer.h>


void DoGenomeGenerate(rna::Parameters& params) {
    rna::Genome g;
    g.binSize = params.genomeBinSize;
    g.loadFromFasta(params.genomeFile);
    rna::GenomeIndex genomeIndex(g);
    genomeIndex.setConfig(params.genomeIndexConfig);
    genomeIndex.build();

    if (!params.gtfFile.empty()) {
        rna::GTF gtf(params.outLogFile);
        gtf.loadGTF(params.gtfFile, *genomeIndex.genome);
        gtf.fillSjdbLoci("./", *genomeIndex.genome);
        gtf.insertJunctions(*genomeIndex.genome, genomeIndex);
    }
    genomeIndex.write(params.genomeGenerateFileStoreDir);
}


void DoReadAlign(rna::Parameters& params) {
    rna::GenomeIndex genomeIndex;
    genomeIndex.setConfig(params.genomeIndexConfig);
    genomeIndex.load(params.genomeGenerateFileStoreDir);

    PLOG_INFO << "Genome index loading done. Will begin alignment.";

    rna::ReadAlignMultiThread readAlignMultiThread(params);
    readAlignMultiThread.processReadFile(params, genomeIndex, false);
}


int main(int argc, char* argv[]){

    plog::init<plog::TxtFormatter>(plog::info, plog::streamStdOut);

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
        PLOG_ERROR << "Please specify the mode: GenomeGenerate, ReadAlign, both or test";
        exit(-1);
    }
    else if (mode == "GenomeGenerate") {
        DoGenomeGenerate(P);
    }
    else if (mode == "ReadAlign") {
        DoReadAlign(P);
    }
    else if (mode == "both") {
        DoGenomeGenerate(P);
        DoReadAlign(P);
    }
    else {
        PLOG_ERROR << "Unknown mode option `" << mode << "`";
        exit(-1);
    }
};