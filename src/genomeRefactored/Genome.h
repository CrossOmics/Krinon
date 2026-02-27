#ifndef RNAALIGNREFACTORED_GENOME_H
#define RNAALIGNREFACTORED_GENOME_H

#include <string>
#include <vector>
#include <map>
#include "../io/Parameters.h"
// todo a placeholder, change the namespace after finishing refactoring
namespace RefactorProcessing {

    class GenomeIndex;

    class SJDB;

    class sjDataBasePiece {
    public:
        size_t start;
        size_t end;
        uint8_t motif;
        uint8_t shiftLeft;
        uint8_t shiftRight;
        uint8_t strand;

    };

    struct Chromosome {
        std::string name;
        int64_t start;
        int64_t length;
    };

    class Genome {
        // ONLY store basic genome information
    private:

        // configs
        int binSizeLog_{18}; //default 18
        std::string genomeFileName_;
        //std::string outPutDir_;


        // data
        std::map<int64_t, int64_t> chromosomeMapStartToIndex_;



        void buildChromosomeMap();


    public:

        std::vector<Chromosome> chromosomes_;
        std::string sequence_;
        size_t genomeLength_;// length of one strand, original genome length. That is, no additional sequence such as sjdb;


        size_t sjdbSeqLength_;// length of sjdb sequence for one strand, that is, sjdb.sjdbLength*sjdb.sjdbNum

        int64_t sjdbStart_;// start position of sjdb sequence in genome sequence

        //additional genome sequence data such as sjdb

        int64_t sjdbNum_;
        std::vector<sjDataBasePiece> sjDataBase_;
        std::vector<int64_t> sjDonorStart_;
        std::vector<int64_t> sjAcceptorStart_;
        Genome() {};

        ~Genome() {};

        // parameters
        void setParam(const Parameters &P);

        // genome generate
        void loadFromFasta(const std::string &file);

        void addComplementaryStrand();

        // seed search

        int64_t getPosChrIndex(int64_t pos) const;


        // modify data
        void modifyGenome(SJDB &sjdb);


        // input/output
        int writeChrInfo(const std::string &fileName) const;

        int writeSjdbInfo(const std::string &fileName) const;

        int writeToFile(const std::string &fileName) const;

        int loadChrInfo(const std::string &fileName);

        int loadSjdbInfo(const std::string &fileName);

        int loadFromFile(const std::string &fileName);

        friend GenomeIndex;
    };
}

#endif //RNAALIGNREFACTORED_GENOME_H
