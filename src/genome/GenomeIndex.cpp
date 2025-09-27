#include <cassert>
#include "GenomeIndex.h"
#include <chrono>
#include <fstream>
#include <filesystem>
#include <plog/Log.h>

namespace rna {
    GenomeIndex::GenomeIndex(Genome &g) : genome(std::make_unique<Genome>(g)) {
        setConfig(GenomeIndexConfig());
    }

    int64_t GenomeIndex::encodeKMer(const std::string_view &pattern, const int mer_length) {
        if (pattern.length() < mer_length) return -1;
        uint32_t code = 0;
        for (char c: pattern.substr(0, mer_length)) {
            if (charToIndex(c) < 0) {
                return -1; // Invalid character, return -1
            }
            code = (code << 2) | charToIndex(c);
        }
        return code;
    }

    void GenomeIndex::setConfig(const GenomeIndexConfig &cfg) {
        config = cfg;
        MER_LENGTH = config.kMerSize;
        MER_NUM = (1 << (MER_LENGTH * 2)) + 1;
    }

    void GenomeIndex::build() {
        std::string revSeq;
        revSeq.reserve(genome->sequence_.length());
        for (int64_t i = genome->sequence_.length() - 1; i >= 0; --i) {
            char c = genome->sequence_[i];
            if (c == 'A') revSeq += 'T';
            else if (c == 'T') revSeq += 'A';
            else if (c == 'C') revSeq += 'G';
            else if (c == 'G') revSeq += 'C';
            else revSeq += c; // keep N or other characters unchanged
        }
        genomeLength = genome->sequence_.length();
        genome->sequence_ = genome->sequence_ + revSeq + "#########################"; // add padding characters
        size_t reservedLength = 0;
        if (config.insertSJ){
            reservedLength += config.limitSjdbInsertN * config.sjdbLength * 2 + 1;
        }
        if (config.insertSecondPass){
            reservedLength += config.limitSjdbInsertN * config.sjdbLength * 2 + 1;
        }
        genome->originalGenomeLength = genome->sequence_.length();

        suffixArray.build(genome->sequence_, true, reservedLength);


        PLOG_INFO << "finished building suffix array";
        patternMerMap_.resize(MER_NUM);
        longestCommonPrefix_.resize(suffixArray.length() + reservedLength, 0);
        buildLCP();
        PLOG_INFO << "finished building LCP";
        buildKMerIndex();
        fillKMerIndex();
        PLOG_INFO << "finished building kMerIndex";
        appearance_flag.clear();
        buildKMerMap();
        PLOG_INFO << "finished building index";
    }

    void GenomeIndex::buildLCP() {
        PackedArray rk;
        PackedArray& sa = suffixArray.sa_;
        rk.init(suffixArray.fullLength_, sa.getWordLengthBits());
        for (size_t i = 0; i < suffixArray.fullLength_; ++i) rk.set(i, 0);
        for (size_t i = 0; i < sa.length(); ++i) rk.set(sa.get(i), i);
        size_t k = 0;
        for (int64_t i = 0; i < suffixArray.fullLength_; ++i) {
            size_t j = rk.get(i);
            if (j == 0) {
                k = 0;
                continue;
            }
            if (k > 0) k--;
            while (genome->sequence_[i + k] == genome->sequence_[sa.get(j - 1) + k]) k++;
            if (k < std::numeric_limits<uint16_t>::max()) longestCommonPrefix_[j] = k;
            else longestCommonPrefix_[j] = std::numeric_limits<uint16_t>::max();
        }
    }

    inline void GenomeIndex::setFlag(int32_t hash, int32_t length) {
        if (length == 1) {
            appearance_flag[length][0] |= (1 << hash);
        } else {
            int32_t pos = hash >> 3;
            int32_t bit = hash & 0x7;
            appearance_flag[length][pos] |= (1 << bit);
        }
    }

    inline bool GenomeIndex::getFlag(int32_t hash, int32_t length) {
        if (length == 1) {
            return (appearance_flag[length][0] >> hash) & 1;
        } else {
            int32_t pos = hash >> 3;
            int32_t bit = hash & 0x7;
            return (appearance_flag[length][pos] >> bit) & 1;
        }
    }

    void GenomeIndex::scanGenome() {
        // initialize appearance_flag
        const std::string &seq = genome->sequence_;
        appearance_flag.resize(MER_LENGTH + 1); // 0 is empty, for convenience
        for (int i = 1; i <= MER_LENGTH; ++i) {
            appearance_flag[i].resize(1 << (std::max(2 * i - 3, 0)), 0);
        }
        nowWindowHash.resize(MER_LENGTH + 1, 0); // Also, 0 is empty
        int nowAvailableLength = 0;
        for (size_t i = 0; i < seq.length(); ++i) {
            charToIndex(seq[i]) >= 0 ? ++nowAvailableLength : nowAvailableLength = 0;
            if (nowAvailableLength == 0) continue;
            if (nowAvailableLength <= MER_LENGTH) {
                nowWindowHash[nowAvailableLength] = (nowWindowHash[nowAvailableLength - 1] << 2) | charToIndex(seq[i]);
                for (int j = 1; j < nowAvailableLength; ++j) {
                    nowWindowHash[j] = (nowWindowHash[j] << 2) | charToIndex(seq[i]);
                    nowWindowHash[j] &= (1 << (2 * j)) - 1; // keep only the last 2*j bits
                }
                for (int j = 1; j <= nowAvailableLength; ++j) {
                    setFlag(nowWindowHash[j], j);
                }
            } else {
                for (int j = 1; j <= MER_LENGTH; ++j) {
                    nowWindowHash[j] = (nowWindowHash[j] << 2) | charToIndex(seq[i]);
                    nowWindowHash[j] &= (1 << (2 * j)) - 1; // keep only the last 2*j bits
                    setFlag(nowWindowHash[j], j);
                }
            }
        }
    }

    void GenomeIndex::buildKMerIndex() {
        //build k-mer index by window scanning
        scanGenome();
        for (size_t i = 0; i < MER_NUM - 1; ++i) {
            int32_t hash = i;
            int32_t length = MER_LENGTH;
            while (!getFlag(hash, length) && length > 0) {
                --length;
                hash >>= 2; // remove the last two bits
            }
            patternMerMap_[i].length = length;
        }
    }

    void GenomeIndex::fillKMerIndex() {
        const std::string &text_ = genome->sequence_;
        const PackedArray &sa_ = suffixArray.sa_;
        int64_t prevHash = -1;
        for (int64_t i = 0; i < sa_.length(); ++i) {
            int64_t hash = encodeKMer(text_.substr(sa_.get(i), MER_LENGTH), MER_LENGTH);
            if (hash > prevHash) {
                if (prevHash != -1 && patternMerMap_[prevHash].upperRange == -1)
                    patternMerMap_[prevHash].upperRange = (int32_t) (i - 1 - patternMerMap_[prevHash].leftSAIndex);
                prevHash = hash;
            }
            if (hash == -1) {
                if (prevHash != -1 && patternMerMap_[prevHash].upperRange == -1)
                    patternMerMap_[prevHash].upperRange = (int32_t) (i - 1 - patternMerMap_[prevHash].leftSAIndex);
            }
            if (patternMerMap_[hash].leftSAIndex == -1)
                patternMerMap_[hash].leftSAIndex = i;
        }

        if (prevHash != -1 && patternMerMap_[prevHash].upperRange == -1)
            patternMerMap_[prevHash].upperRange = (int32_t) (sa_.length() - 1 - patternMerMap_[prevHash].leftSAIndex);

        int64_t prevRightBound = 0;

        for (int64_t i = 0; i < MER_NUM - 1; --i) {
            if (patternMerMap_[i].length < MER_LENGTH) {
                patternMerMap_[i].leftSAIndex = prevRightBound + 1;
            } else {
                prevRightBound = patternMerMap_[i].leftSAIndex + patternMerMap_[i].upperRange;
            }
        }
    }


    std::pair<size_t, bool>
    GenomeIndex::matchGenomeSeq(size_t matchedLength, size_t pos, const std::string_view &pattern) const {
        size_t l;
        bool greater = false;
        for (l = matchedLength; l < pattern.length(); ++l) {
            if (pattern[l] != genome->sequence_[pos + l]) {
                if (charToIndex(genome->sequence_[pos + l]) < 0)
                    greater = true;// 'N' and '#' are regarded as the largest
                else greater = pattern[l] < genome->sequence_[pos + l];
                break;
            }
        }
        return {l, greater};
    }

    // Note: [lBound, rBound)


    MatchStats GenomeIndex::getMatchStatsLCP(const std::string_view &pattern, size_t rangeLeft, size_t rangeRight,
                                             size_t matchedLength) const {
        //first ,get the longest match
        if (rangeLeft >= rangeRight) {
            return {0, 0, 0, 0};
        }
        rangeRight--;
        if (rangeLeft == rangeRight) {
            auto [longestLength, greater] = matchGenomeSeq(matchedLength, suffixArray.sa_.get(rangeLeft),
                                                           pattern);
            return {int64_t(longestLength), 1, rangeLeft, rangeLeft};
        }
        size_t leftMaxLength = matchedLength;
        size_t rightMaxLength = matchedLength;
        size_t cnt = 0;
        while (rangeRight - rangeLeft > 1) {
            size_t mid = (rangeLeft + rangeRight) / 2;
            size_t midPos = suffixArray.sa_.get(mid);
            auto [midLength, midGreater] = matchGenomeSeq(matchedLength, midPos, pattern);
            if (midLength == pattern.length()) {
                size_t leftPos = mid;
                size_t rightPos = mid;
                // left range
                int64_t leftCommonPrefix = longestCommonPrefix_[leftPos];
                while (leftCommonPrefix >= midLength) {
                    --leftPos;
                    leftCommonPrefix = longestCommonPrefix_[leftPos];
                }
                // right range
                int64_t rightCommonPrefix = longestCommonPrefix_[rightPos + 1];
                while (rightCommonPrefix >= midLength) {
                    ++rightPos;
                    rightCommonPrefix = longestCommonPrefix_[rightPos + 1];
                }
                return {(int64_t) midLength, rightPos - leftPos + 1, leftPos, rightPos};
            }
            if (midGreater) {
                rangeRight = mid;
                rightMaxLength = midLength;
            } else {
                rangeLeft = mid;
                leftMaxLength = midLength;
            }
            matchedLength = std::min(leftMaxLength, rightMaxLength);
        }
        //now, rangeLeft == rangeRight - 1
        auto [lLength, leftGreater] = matchGenomeSeq(leftMaxLength, suffixArray.sa_.get(rangeLeft), pattern);
        auto [rLength, rightGreater] = matchGenomeSeq(rightMaxLength, suffixArray.sa_.get(rangeRight), pattern);
        size_t longestLength = std::max(lLength, rLength);
        size_t leftPos = rangeLeft;
        size_t rightPos = rangeRight;
        //left range
        if (lLength == longestLength) {
            leftPos = rangeLeft;
            int64_t leftCommonPrefix = longestCommonPrefix_[leftPos];
            ++cnt;
            while (leftCommonPrefix >= longestLength) {
                --leftPos;
                ++cnt;
                if (cnt > 10000) return {(int64_t) longestLength, 0, 0, 0};
                leftCommonPrefix = longestCommonPrefix_[leftPos];
            }
        } else leftPos = rangeLeft + 1;

        //right range
        if (rLength == longestLength) {
            rightPos = rangeRight;
            int64_t rightCommonPrefix = longestCommonPrefix_[rightPos + 1];
            ++cnt;
            while (rightCommonPrefix >= longestLength) {
                ++rightPos;
                ++cnt;
                if (cnt > 10000) return {(int64_t) longestLength, 0, 0, 0};
                rightCommonPrefix = longestCommonPrefix_[rightPos + 1];
            }
        } else rightPos = rangeRight - 1;

        return {(int64_t) longestLength, rightPos - leftPos + 1, leftPos, rightPos};
    }


    int GenomeIndex::alignRead(const std::string *readSeq, const Split &split, Align *aligns, int &alignNum,
                               int maxAlignPerRead) const {
        // align the read to the genome, return the alignments

        int64_t length;
        int64_t iStart = 1 + (split.length / 50);
        int64_t lStart = split.length / iStart;

        bool fullMatch = false;
        for(int dir = 0; dir <= 1; ++dir){
            int64_t splitReadStart = dir == 0 ? split.readStart : (readSeq[0].length() - split.readStart - split.length);
            for (int i = 0; i < iStart; ++i) {
                int64_t nowMappedLength = i * lStart;
                do {
                    MatchStats matchStats;
                    if (split.length - nowMappedLength < MER_LENGTH && split.length - nowMappedLength != 0) {
                        if (split.length -nowMappedLength <= 5) break; // too short to align
                        matchStats = find(readSeq[dir].substr(splitReadStart + split.length - MER_LENGTH, MER_LENGTH));
                        if (matchStats.count > 0) {
                            int posToInsert = -1;
                            int64_t readStart = splitReadStart + split.length - MER_LENGTH;
                            int64_t nowLength = matchStats.maxLength;
                            for (int j = 0; j < alignNum; ++j) {
                                if (aligns[j].readStart < readStart) continue;
                                if (aligns[j].readStart == readStart) {
                                    if (aligns[j].length < nowLength) continue;
                                    if (aligns[j].length == nowLength) {
                                        if (dir != aligns[j].direction) continue;
                                        else {
                                            posToInsert = -2;
                                            break; // already exists
                                        }
                                    }
                                }
                                posToInsert = j;
                            }
                            if (posToInsert == -1) posToInsert = alignNum;
                            if (posToInsert != -2 && alignNum < maxAlignPerRead) {
                                for (int j = alignNum; j > posToInsert; --j) {
                                    aligns[j] = aligns[j - 1];
                                }
                                aligns[posToInsert] = Align{
                                        .readStart = readStart,
                                        .length = matchStats.maxLength,
                                        .leftSAIndex = matchStats.LongestMatchLeftPosInSA,
                                        .rightSAIndex = matchStats.LongestMatchRightPosInSA,
                                        .rep = matchStats.count,
                                        .direction = dir,
                                        .isAnchor = false,
                                        .iFragment = split.iFragment
                                };
                                alignNum++;
                            }
                        }
                    } else {
                        matchStats = find(
                                readSeq[dir].substr(splitReadStart + nowMappedLength, split.length - nowMappedLength));
                        length = matchStats.maxLength;
                        int64_t readStart = splitReadStart + nowMappedLength;

                        int posToInsert = -1;
                        for (int j = 0; j < alignNum; ++j) {
                            if (aligns[j].readStart < readStart) continue;
                            if (aligns[j].readStart == readStart) {
                                if (aligns[j].length < length) continue;
                                if (aligns[j].length == length) {
                                    if (dir != aligns[j].direction) continue;
                                    else {
                                        posToInsert = -2;
                                        break; // already exists
                                    }
                                }
                            }
                            posToInsert = j;
                        }
                        if (posToInsert == -1) posToInsert = alignNum;
                        if (posToInsert != -2 && alignNum < maxAlignPerRead) {
                            for (int j = alignNum; j > posToInsert; --j) {
                                aligns[j] = aligns[j - 1];
                            }
                            aligns[posToInsert] = Align{
                                    .readStart = readStart,
                                    .length = matchStats.maxLength,
                                    .leftSAIndex = matchStats.LongestMatchLeftPosInSA,
                                    .rightSAIndex = matchStats.LongestMatchRightPosInSA,
                                    .rep = matchStats.count,
                                    .direction = dir,
                                    .isAnchor = false,
                                    .iFragment = split.iFragment
                            };
                            alignNum++;
                        }
                    }

                    if (length == split.length) {
                        fullMatch = true;
                        break;
                    }
                    nowMappedLength += length;
                    if (nowMappedLength >= split.length) {
                        break; // no more to align
                    }
                } while (length > 0);
                if (fullMatch) break; // no need to continue, already found a full match
            }
            if (fullMatch) break; // no need to continue, already found a full match
        }

        return 0;
    }

    void GenomeIndex::write(const std::string &outDir) const {
        if (!std::filesystem::exists(outDir)) std::filesystem::create_directories(outDir);
        std::ofstream logFile(outDir + "/indexInf.txt");
        std::ofstream genomeFile(outDir + "/genome.bin", std::ios::binary | std::ios::out);
        std::ofstream saFile(outDir + "/SA.bin", std::ios::binary | std::ios::out);
        std::ofstream saIndexFile(outDir + "/SAIndex.bin", std::ios::binary | std::ios::out);
        //write genome sequence
        logFile << genomeLength << "\n";
        logFile << genome->sequence_.length() << "\n";
        logFile << genome->originalGenomeLength << "\n";
        logFile << genome->sjdbNum << "\n";
        genomeFile.write(genome->sequence_.data(), genome->sequence_.length() * sizeof(char));
        if (genome->sjdbNum > 0) {
            genomeFile.write(reinterpret_cast<const char *>(genome->sjdb.data()),
                             genome->sjdbNum * sizeof(GenomeSjdbRecord));
            genomeFile.write(reinterpret_cast<const char *>(genome->sjDonorStart_.data()),
                             genome->sjdbNum * sizeof(int32_t));
            genomeFile.write(reinterpret_cast<const char *>(genome->sjAcceptorStart_.data()),
                             genome->sjdbNum * sizeof(int32_t));
        }
        genomeFile.close();
        logFile << genome->chromosomes_.size() << '\n';
        for (auto &chromosome: genome->chromosomes_) {
            logFile << chromosome.name << " " << chromosome.start << ' ' << chromosome.length << "\n";
        }
        //write suffix array
        logFile << suffixArray.length() << " " << suffixArray.fullLength_ << "\n";
        suffixArray.sa_.writeToFile(saFile, logFile);

        //write suffix array index
        logFile << config.kMerSize << " " << config.twoDirections << "\n";

        logFile << patternMerMap_.size() << "\n";
        saIndexFile.write(reinterpret_cast<const char *>(patternMerMap_.data()),
                          patternMerMap_.size() * sizeof(MatchResultStored));
        logFile << longestCommonPrefix_.size() << "\n";
        saIndexFile.write(reinterpret_cast<const char *>(longestCommonPrefix_.data()),
                          longestCommonPrefix_.size() * sizeof(uint16_t));
        logFile << extendIndexHash.size() << "\n";
        saIndexFile.write(reinterpret_cast<const char *>(extendIndexHash.data()),
                          extendIndexHash.size() * sizeof(uint32_t));
        logFile.close();
        saFile.close();
        saIndexFile.close();
    }

    void GenomeIndex::load(const std::string &inDir) {
        genome = std::make_unique<Genome>(Genome());

        std::ifstream logFile(inDir + "/indexInf.txt");
        std::ifstream genomeFile(inDir + "/genome.bin", std::ios::binary | std::ios::in);
        std::ifstream saFile(inDir + "/SA.bin", std::ios::binary | std::ios::in);
        std::ifstream saIndexFile(inDir + "/SAIndex.bin", std::ios::binary | std::ios::in);
        
        if (!logFile.is_open() || !genomeFile.is_open() || !saFile.is_open() || !saIndexFile.is_open()) {
            throw std::runtime_error("Failed to open index files");
        }
        
        logFile >> genomeLength;
        int64_t genomeSeqLength;
        logFile >> genomeSeqLength;
        logFile >> genome->originalGenomeLength;
        logFile >> genome->sjdbNum;

        genome->sequence_.resize(genomeSeqLength);
        genomeFile.read(genome->sequence_.data(), genomeSeqLength * sizeof(char));
        if (genome->sjdbNum > 0) {
            genome->sjdb.resize(genome->sjdbNum);
            genome->sjDonorStart_.resize(genome->sjdbNum);
            genome->sjAcceptorStart_.resize(genome->sjdbNum);
            genomeFile.read(reinterpret_cast<char *>(genome->sjdb.data()),
                            genome->sjdbNum * sizeof(GenomeSjdbRecord));
            genomeFile.read(reinterpret_cast<char *>(genome->sjDonorStart_.data()),
                            genome->sjdbNum * sizeof(int32_t));
            genomeFile.read(reinterpret_cast<char *>(genome->sjAcceptorStart_.data()),
                            genome->sjdbNum * sizeof(int32_t));
        }
        genomeFile.close();
        PLOG_INFO << "Finished loading binary genome file (length: " << genomeLength << ")";

        genome->chromosomes_.clear();
        std::string chrName;
        int64_t chrStart, chrLength, chrNum;
        logFile >> chrNum;
        for (int64_t i = 0; i < chrNum; ++i) {
            logFile >> chrName >> chrStart >> chrLength;
            genome->chromosomes_.emplace_back(Genome::Chromosome{chrName, chrStart, chrLength});
            genome->chromosomeStartMap_[chrStart] = i;
        }
        int64_t saLength, saFullLength;
        logFile >> saLength >> suffixArray.fullLength_;
        suffixArray.sa_.loadFromFile(saFile, logFile);
        PLOG_INFO << "Finished loading suffix array (size: " << saLength << ")";

        logFile >> config.kMerSize >> config.twoDirections;
        setConfig(GenomeIndexConfig{config.kMerSize, config.twoDirections});

        int64_t patternMerMapSize;
        logFile >> patternMerMapSize;
        patternMerMap_.resize(patternMerMapSize);
        saIndexFile.read(reinterpret_cast<char *>(patternMerMap_.data()),
                         patternMerMapSize * sizeof(MatchResultStored));
        int64_t lcpSize;
        logFile >> lcpSize;
        longestCommonPrefix_.resize(lcpSize);
        saIndexFile.read(reinterpret_cast<char *>(longestCommonPrefix_.data()),
                         lcpSize * sizeof(uint16_t));
        int64_t extendIndexHashSize;
        logFile >> extendIndexHashSize;
        extendIndexHash.resize(extendIndexHashSize, 0);
        saIndexFile.read(reinterpret_cast<char *>(extendIndexHash.data()),
                         extendIndexHashSize * sizeof(uint32_t));
        PLOG_INFO << "Finished loading suffix array indices";

        logFile.close();
        saFile.close();
        saIndexFile.close();
    }

    void GenomeIndex::buildKMerMapSingle(int64_t h){
        int64_t nowLeftSAIndex = patternMerMap_[h].leftSAIndex;
        int64_t num = patternMerMap_[h].upperRange;
        int indNum = config.extendAlternativeByte / 4;
        if (num > 16) {
            // 4 Bytes per index
            for (int j = 0; j < indNum; ++j) {
                int length = 16;
                uint64_t index = suffixArray[nowLeftSAIndex +
                                                     j * num * 4 / config.extendAlternativeByte];
                uint64_t hash = 0;

                for (int k = 0; k < length; ++k) {
                    int c = charToIndex(genome->sequence_[index + MER_LENGTH + k]);
                    if (c != -1) {
                        hash <<= 2;
                        hash |= c;
                    } else {
                        // set it to xxTT
                        hash <<= (2 * (length - k));
                        hash |= (1 << (2 * (length - k))) - 1;
                        break;
                    }
                }
                extendIndexHash[h * indNum + j] = (uint32_t) hash;
            }
        }
    }

    void GenomeIndex::buildKMerMap() {
        extendIndexHash.resize(MER_NUM * config.extendAlternativeByte / 4, 0);

        for (size_t i = 0; i < MER_NUM - 1; ++i) {
            buildKMerMapSingle(i);
        }
    }

    MatchStats GenomeIndex::find(const std::string_view &pattern) const {

        if (pattern.length() < MER_LENGTH) {
            //only occurs on the last few bps of read
            int l = pattern.length();
            int64_t hash = encodeKMer(pattern, l);
            if (hash == -1) return {0, 0, 0, 0};


            int64_t leftHash = hash << ((MER_LENGTH - l) * 2);
            int64_t rightHash = leftHash | ((1 << (2 * (MER_LENGTH - l))) - 1);
            size_t leftSAIndex = patternMerMap_[leftHash].leftSAIndex;
            size_t rightSAIndex = patternMerMap_[rightHash].leftSAIndex + patternMerMap_[rightHash].upperRange;

            return {l, rightSAIndex - leftSAIndex + 1, leftSAIndex, rightSAIndex};
        }

        // now, find the max matching prefix through the k-mer map

        int64_t hash = encodeKMer(pattern.substr(0, MER_LENGTH), MER_LENGTH);
        int l = patternMerMap_[hash].length;
        if (l < MER_LENGTH || pattern.length() == MER_LENGTH) {
            int64_t leftHash = hash >> ((MER_LENGTH - l) * 2);
            leftHash = leftHash << ((MER_LENGTH - l) * 2);
            int64_t rightHash = leftHash | ((1 << (2 * (MER_LENGTH - l))) - 1);
            size_t leftSAIndex = patternMerMap_[leftHash].leftSAIndex;
            size_t rightSAIndex = patternMerMap_[rightHash].leftSAIndex + patternMerMap_[rightHash].upperRange;
            if (rightSAIndex - leftSAIndex + 1 > 10000) {
                return {0, 0, 0, 0};
            } else {
                int64_t nowRightSAIndex = rightSAIndex;
                while (longestCommonPrefix_[nowRightSAIndex + 1] >= l) {
                    ++nowRightSAIndex;
                    if (nowRightSAIndex >= suffixArray.length()) break;
                }
                rightSAIndex = nowRightSAIndex;
            }
            return {l, rightSAIndex - leftSAIndex + 1, leftSAIndex, rightSAIndex};
        } else {
            size_t lBound = patternMerMap_[hash].leftSAIndex;
            size_t rBound = patternMerMap_[hash].leftSAIndex + patternMerMap_[hash].upperRange + 1;
            int64_t num = patternMerMap_[hash].upperRange;

            if (num <= 16) {
                return getMatchStatsLCP(pattern, lBound, rBound, l);
            }
            const int64_t extendIndexHashBase = hash * config.extendAlternativeByte / 4;
            int64_t remainingLength = pattern.length() - MER_LENGTH;


            if (remainingLength < 16) {
                uint32_t h = encodeKMer(pattern.substr(l, remainingLength), remainingLength);
                uint32_t hLeft = h << ((16 - remainingLength) * 2);
                uint32_t hRight = (h + 1) << ((16 - remainingLength) * 2);
                size_t left = 0, right = config.extendAlternativeByte / 4;
                size_t mid;
                while (right - left > 1) {
                    mid = (left + right) / 2;
                    if (extendIndexHash[extendIndexHashBase + mid] < hLeft) left = mid;
                    else if (extendIndexHash[extendIndexHashBase + mid] > hRight) right = mid;
                    else {
                        // true left bound in [left,mid] ,true right bound in [mid,right]
                        size_t boundLeft = mid;
                        while (boundLeft - left > 1) {
                            size_t midLeft = (left + boundLeft) / 2;
                            if (extendIndexHash[extendIndexHashBase + midLeft] < hLeft) left = midLeft;
                            else if (extendIndexHash[extendIndexHashBase + midLeft] > hLeft) boundLeft = midLeft;
                            else {
                                // exact match
                                left = midLeft;
                                while (left > 0 && extendIndexHash[extendIndexHashBase + left] == hLeft)
                                    --left;
                                break;
                            }
                        }
                        size_t boundRight = mid;
                        while (right - boundRight > 1) {
                            size_t midRight = (boundRight + right) / 2;
                            if (extendIndexHash[extendIndexHashBase + midRight] < hRight) boundRight = midRight;
                            else if (extendIndexHash[extendIndexHashBase + midRight] > hRight) right = midRight;
                            else {
                                // exact match
                                right = midRight;
                                while (right < config.extendAlternativeByte / 4 &&
                                       extendIndexHash[extendIndexHashBase + right] == hRight)
                                    ++right;
                                break;
                            }
                        }
                        // now we get the range
                        break;
                    }
                }

                // reaching here means that [hLeft, hRight] is in [extendIndexHash[extendIndexHashBase + left], extendIndexHash[extendIndexHashBase + right]]

                rBound = lBound + right * num * 4 / config.extendAlternativeByte;
                if (right != config.extendAlternativeByte / 4) ++rBound;
                lBound += left * num * 4 / config.extendAlternativeByte;

                return getMatchStatsLCP(pattern, lBound, rBound, l);
            }
            if (num > 16) {
                uint32_t h = encodeKMer(pattern.substr(l, 16), 16);
                size_t left = 0, right = config.extendAlternativeByte / 4;
                while (right - left > 1) {
                    size_t mid = (left + right) / 2;
                    if (extendIndexHash[extendIndexHashBase + mid] < h) left = mid;
                    else if (extendIndexHash[extendIndexHashBase + mid] > h) right = mid;
                    else {
                        // exact match
                        left = right = mid;
                        while (left > 0 && extendIndexHash[extendIndexHashBase + left] == h)
                            --left;
                        while (right < config.extendAlternativeByte / 4 &&
                               extendIndexHash[extendIndexHashBase + right] == h)
                            ++right;
                        break;
                    }
                }
                rBound = lBound + right * num * 4 / config.extendAlternativeByte;
                if (right != config.extendAlternativeByte / 4) ++rBound;
                lBound += left * num * 4 / config.extendAlternativeByte;

            } else {
                // no need to search
            }

            auto matchStats = getMatchStatsLCP(pattern, lBound, rBound, l);

            return matchStats;
        }
    }
} // namespace rna
