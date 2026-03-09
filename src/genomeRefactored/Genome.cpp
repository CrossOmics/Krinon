#include "Genome.h"
#include "../log/ErrorRecord.h"
#include "../utilsRefactored/defines.h"
#include "../utilsRefactored/seqFunctions.h"
#include <fstream>
namespace RefactorProcessing {
    void Genome::setParam(const Parameters &P) {
        binSizeLog_ = P.genomeBinSize;
        genomeFileName_ = P.genomeFile;
        sjdbNum_ = 0;
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

        genomeLength_ += 10; // half of the spacing characters, for safety

    }




    void Genome::loadFromFasta(const std::string &filename) {
        //todo optimize: load the whole genome in one go, use O_DIRECT flag
        std::ifstream file(filename);
        if (!file) rna::ErrorRecord().reportError("Cannot open genome file: " + filename);
        std::string line;
        std::stringstream currentSeq;
        std::string currentChr;
        size_t currentStart = 0;

        while (std::getline(file, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (!currentChr.empty()) {
                    size_t currentPos = static_cast<size_t>(currentSeq.tellp());
                    size_t chrLen = currentPos - currentStart;

                    size_t paddingLen = ((chrLen + (1 << binSizeLog_)) & ~((1 << binSizeLog_) - 1)) - chrLen;

                    for (size_t i = 0; i < paddingLen; ++i) {
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

            size_t currentPos = static_cast<size_t>(currentSeq.tellp());
            size_t chrLen = currentPos - currentStart;
            size_t paddingLen = ((chrLen + (1 << binSizeLog_) - 1) & ~((1 << binSizeLog_) - 1)) - chrLen;
            for (size_t i = 0; i < paddingLen; ++i) {
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

    int Genome::writeChrInfo(const std::string &fileName) const {
        std::ofstream ofs(fileName, std::ios::binary);
        if (!ofs.is_open()) return -1;
        uint64_t chrNum = chromosomes_.size();
        ofs.write(reinterpret_cast<const char*>(&chrNum), sizeof(chrNum));
        for (const auto& chr : chromosomes_) {
            if (!writeString(ofs, chr.name)) return -2;
            ofs.write(reinterpret_cast<const char*>(&chr.start), sizeof(chr.start));
            ofs.write(reinterpret_cast<const char*>(&chr.length), sizeof(chr.length));
        }
        return 0;
    }

    int Genome::writeToFile(const std::string &fileName) const {
        std::ofstream ofs(fileName, std::ios::binary);
        if (!ofs.is_open()) return -1;
        if(writeChrInfo(fileName+".chr") != 0) return -2;
        if(writeSjdbInfo(fileName+".sjdb") != 0) return -3;


        ofs.write(reinterpret_cast<const char*>(&genomeLength_), sizeof(genomeLength_));
        // write sequence
        writeString(ofs, sequence_);
        return 0;
    }

    int Genome::loadChrInfo(const std::string &fileName) {
        std::ifstream ifs(fileName, std::ios::binary);
        if (!ifs.is_open()) return -1;
        uint64_t chrNum;
        ifs.read(reinterpret_cast<char*>(&chrNum), sizeof(chrNum));
        if (!ifs) return -2;
        chromosomes_.clear();
        for (size_t i = 0; i < chrNum; ++i) {
            Chromosome chr;
            if (!readString(ifs, chr.name)) return -3;
            ifs.read(reinterpret_cast<char*>(&chr.start), sizeof(chr.start));
            ifs.read(reinterpret_cast<char*>(&chr.length), sizeof(chr.length));
            if (!ifs) return -4;
            chromosomes_.push_back(chr);
        }
        buildChromosomeMap();
        return 0;
    }

    int Genome::loadFromFile(const std::string &fileName) {
        std::ifstream ifs(fileName, std::ios::binary);
        if (!ifs.is_open()) return -1;
        ifs.read(reinterpret_cast<char*>(&genomeLength_), sizeof(genomeLength_));
        if (!ifs) return -2;
        if (!readString(ifs, sequence_)) return -3;
        if (loadChrInfo(fileName + ".chr") != 0) return -4;
        if (loadSjdbInfo(fileName + ".sjdb") != 0) return -5;
        return 0;
    }

    int Genome::writeSjdbInfo(const std::string &fileName) const {
        std::ofstream ofs(fileName, std::ios::binary);
        if (!ofs.is_open()) return -1;
        ofs.write(reinterpret_cast<const char*>(&sjdbNum_), sizeof(sjdbNum_));
        ofs.write(reinterpret_cast<const char*>(&sjdbStart_), sizeof(sjdbStart_));
        ofs.write(reinterpret_cast<const char*>(&sjdbSeqLength_), sizeof(sjdbSeqLength_));
        //write vectors
        if (sjdbNum_ > 0){
            ofs.write(reinterpret_cast<const char*>(sjDataBase_.data()), sjDataBase_.size() * sizeof(sjDataBasePiece));
            ofs.write(reinterpret_cast<const char*>(sjDonorStart_.data()), sjDonorStart_.size() * sizeof(int64_t));
            ofs.write(reinterpret_cast<const char*>(sjAcceptorStart_.data()), sjAcceptorStart_.size() * sizeof(int64_t));
        }

        return 0;
    }

    int Genome::loadSjdbInfo(const std::string &fileName) {
        std::ifstream ifs(fileName, std::ios::binary);
        if (!ifs.is_open()) return -1;
        ifs.read(reinterpret_cast<char*>(&sjdbNum_), sizeof(sjdbNum_));
        ifs.read(reinterpret_cast<char*>(&sjdbStart_), sizeof(sjdbStart_));
        ifs.read(reinterpret_cast<char*>(&sjdbSeqLength_), sizeof(sjdbSeqLength_));
        //read vectors
        if (sjdbNum_ > 0){
            sjDataBase_.resize(sjdbNum_);
            sjDonorStart_.resize(sjdbNum_);
            sjAcceptorStart_.resize(sjdbNum_);
            ifs.read(reinterpret_cast<char*>(sjDataBase_.data()), sjDataBase_.size() * sizeof(sjDataBasePiece));
            ifs.read(reinterpret_cast<char*>(sjDonorStart_.data()), sjDonorStart_.size() * sizeof(int64_t));
            ifs.read(reinterpret_cast<char*>(sjAcceptorStart_.data()), sjAcceptorStart_.size() * sizeof(int64_t));
        }

        return 0;
    }


}