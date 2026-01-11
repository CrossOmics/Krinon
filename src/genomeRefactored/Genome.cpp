#include "Genome.h"
#include "../log/ErrorRecord.h"
#include "../utilsRefactored/defines.h"
#include <fstream>
namespace RefactorProcessing {
    void Genome::setParam(const Parameters &P) {
        binSizeLog_ = P.genomeBinSize;
        genomeFileName_ = P.genomeFile;
    }

    void Genome::addComplementaryStrand() {
        std::string revSeq;
        revSeq.reserve(sequence_.length());
        for (int64_t i = genomeLength_ - 1; i >= 0; --i) {
            char c = sequence_[i];
            if (c == 'A') revSeq += 'T';
            else if (c == 'T') revSeq += 'A';
            else if (c == 'C') revSeq += 'G';
            else if (c == 'G') revSeq += 'C';
            else revSeq += c; // keep N or other characters unchanged
        }
        sequence_ = sequence_ + "####################" + revSeq + "####################";

    }




    void Genome::loadFromFasta(const std::string &filename) {
        //todo optimize: load the whole genome in one go, use O_DIRECT flag
        std::ifstream file(filename);
        if (!file) rna::ErrorRecord().reportError("Cannot open genome file: " + filename);
        std::string line;
        std::stringstream currentSeq;
        std::string currentChr;
        int64_t currentStart = 0;

        while (std::getline(file, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (!currentChr.empty()) {
                    int64_t currentPos = static_cast<int64_t>(currentSeq.tellp());
                    int64_t chrLen = currentPos - currentStart;

                    int64_t paddingLen = ((chrLen + (1 << binSizeLog_)) & ~((1 << binSizeLog_) - 1)) - chrLen;

                    for (int64_t i = 0; i < paddingLen; ++i) {
                        currentSeq << SPACING_CHAR;
                    }

                    chromosomes_.push_back(Chromosome{
                            currentChr,
                            currentStart,
                            chrLen
                    });
                    currentStart = static_cast<int64_t>(currentSeq.tellp());
                } else {
                    currentStart = 0;
                }

                size_t nameEnd = line.find_first_of(" \t", 1);
                currentChr = line.substr(1, nameEnd - 1);
                if (currentChr[currentChr.size() - 1] == '\r') {
                    currentChr.pop_back(); // Windows file in linux
                }

            } else {
                line.erase(
                        std::remove_if(line.begin(), line.end(), ::isspace),
                        line.end()
                );
                std::transform(line.begin(), line.end(), line.begin(), ::toupper);
                currentSeq << line;
            }
        }


        if (!currentChr.empty()) {

            int64_t currentPos = static_cast<int64_t>(currentSeq.tellp());
            int64_t chrLen = currentPos - currentStart;
            int64_t paddingLen = ((chrLen + (1 << binSizeLog_) - 1) & ~((1 << binSizeLog_) - 1)) - chrLen;
            for (int64_t i = 0; i < paddingLen; ++i) {
                currentSeq << SPACING_CHAR;
            }
            chromosomes_.push_back(Chromosome{
                    currentChr,
                    currentStart,
                    chrLen
            });
        }
        file.close();
        sequence_ = currentSeq.str();
        genomeLength_ = sequence_.length();
        buildChromosomeMap();
        if (sequence_.empty()) {
            rna::ErrorRecord().reportError("No sequence loaded from genome file: " + filename);
        }

        addComplementaryStrand();
    }

    void Genome::buildChromosomeMap() {
        chromosomeMapStartToIndex_.clear();
        for (size_t i = 0; i < chromosomes_.size(); ++i) {
            chromosomeMapStartToIndex_[chromosomes_[i].start] = i;
        }
    }

    int64_t Genome::getPosChrIndex(int64_t pos) const {
        return (--chromosomeMapStartToIndex_.upper_bound(pos))->second;
    }


}