#include "Genome.h"
#include "../utils/exceptions.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>

namespace rna{
    Genome::Genome() = default;
    Genome::Genome(const std::string& filename) {
        loadFromFasta(filename);
    }
    void Genome::loadFromFasta(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) {
            throw RNAException(ExceptionType::FILE_NOT_FOUND,
                               "Cannot open genome file: " + filename);
        }
        std::string line;
        std::stringstream currentSeq;
        std::string currentChr;
        GenomePos currentStart = 0;
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (!currentChr.empty()) {
                    GenomePos currentPos = static_cast<GenomePos>(currentSeq.tellp());
                    GenomePos chrLen = currentPos - currentStart;
                    // 计算需要补充的长度，使总长为2^INDEX_SHIFT的整数倍
                    GenomePos paddingLen = ((chrLen + (1 << INDEX_SHIFT)) & ~((1 << INDEX_SHIFT) - 1)) - chrLen;

                    for (GenomePos i = 0; i < paddingLen; ++i) {
                        currentSeq << SPACING_CHAR;
                    }

                    chromosomes_.push_back(Chromosome{
                            currentChr,
                            currentStart,
                            chrLen
                    });
                    currentStart = static_cast<GenomePos>(currentSeq.tellp());
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

            GenomePos currentPos = static_cast<GenomePos>(currentSeq.tellp());
            GenomePos chrLen = currentPos - currentStart;
            GenomePos paddingLen = ((chrLen + (1 << INDEX_SHIFT) - 1) & ~((1 << INDEX_SHIFT) - 1)) - chrLen;
            for (GenomePos i = 0; i < paddingLen; ++i) {
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
        buildPosToChromosomeMap();
        if (sequence_.empty()) {
            throw RNAException(ExceptionType::INVALID_INPUT,
                               "Empty genome sequence in file: " + filename);
        }
    }
    std::string_view Genome::getSequence(size_t pos,size_t length) const {
        if (pos + length > sequence_.length()) {
            throw RNAException(ExceptionType::INVALID_INPUT,
                               "Genome sequence access out of range");
        }
        return std::string_view(sequence_).substr(pos, length);
    }
    void Genome::buildPosToChromosomeMap() {
        for(int64_t i = 0; i < chromosomes_.size(); ++i) {
            chromosomeMap_[chromosomes_[i].start] = i;
        }
    }
    // return 0-based position in chromosome
    std::pair<std::string,int64_t> Genome::getPosChromosome(int64_t pos) const {
        Chromosome c = chromosomes_[(--chromosomeMap_.upper_bound(pos))->second ];
        return {c.name, pos - c.start};
    }

    int64_t Genome::getPosChrIndex(int64_t pos) const {
        return (--chromosomeMap_.upper_bound(pos))->second;
    }

    void Genome::writeChrInfo(const std::string& dirOut) const {

        std::filesystem::create_directories(dirOut);

        {
            std::ofstream chrNameFile(dirOut + "/chrName.txt");
            if (!chrNameFile) {
                throw RNAException(ExceptionType::FILE_NOT_FOUND,
                                   "Cannot create file: " + dirOut + "/chrName.txt");
            }

            for (const auto& chr : chromosomes_) {
                chrNameFile << chr.name << '\n';
            }
            chrNameFile.close();
        }

        {
            std::ofstream chrLengthFile(dirOut + "/chrLength.txt");
            if (!chrLengthFile) {
                throw RNAException(ExceptionType::FILE_NOT_FOUND,
                                   "Cannot create file: " + dirOut + "/chrLength.txt");
            }

            for (const auto& chr : chromosomes_) {
                chrLengthFile << chr.length << '\n';
            }
            chrLengthFile.close();
        }

        {
            std::ofstream chrStartFile(dirOut + "/chrStart.txt");
            if (!chrStartFile) {
                throw RNAException(ExceptionType::FILE_NOT_FOUND,
                                   "Cannot create file: " + dirOut + "/chrStart.txt");
            }

            for (const auto& chr : chromosomes_) {
                chrStartFile << chr.start << '\n';
            }
            chrStartFile.close();
        }
    }
}