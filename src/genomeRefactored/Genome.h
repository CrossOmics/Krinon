#ifndef RNAALIGNREFACTORED_GENOME_H
#define RNAALIGNREFACTORED_GENOME_H

#include <string>
#include <vector>
#include <map>
#include "../io/Parameters.h"
// todo a placeholder, change the namespace after finishing refactoring
namespace RefactorProcessing {
    class AdditionalGenomeSegment;
    class GenomeIndex;
    class Genome{
        struct Chromosome {
            std::string name;
            int64_t start;
            int64_t length;
        };
        // ONLY store basic genome information
    private:

        // configs
        int binSizeLog_{18}; //default 18
        std::string genomeFileName_;
        std::string outPutDir_;


        // data
        std::string sequence_;
        size_t genomeLength_;// length of one strand, original genome length. That is, no additional sequence such as sjdb;
        std::vector<Chromosome> chromosomes_;
        std::map<int64_t,int64_t> chromosomeMapStartToIndex_;

        //additional genome sequence data such as sjdb
        int64_t additionalSegmentNum_;
        std::vector<AdditionalGenomeSegment*> additionalSegments_;

        void buildChromosomeMap();


    public:

        friend GenomeIndex;

        Genome(){};
        ~Genome(){};

        // parameters
        void setParam(const Parameters& P);

        // genome generate
        void loadFromFasta(const std::string &file);

        void addComplementaryStrand();

        // seed search

        int64_t getPosChrIndex(int64_t pos) const;


        // modify data
        void modifyGenome();


        // input/output
        void writeChrInfo(const std::string& dirOut) const;
        int writeToFile(const std::string& dirName);
        int loadFromFile(const std::string& dirName);


    };
}

#endif //RNAALIGNREFACTORED_GENOME_H
