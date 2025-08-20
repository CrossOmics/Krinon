#include <cassert>
#include "GenomeIndexPrefix.h"
#include <chrono>
#include <fstream>
#include <filesystem>

int LCPSearchCount = 0;
namespace rna {
    GenomeIndexPrefix::GenomeIndexPrefix(Genome &g) : genome(std::make_unique<Genome>(g)) {
        setConfig(GenomeIndexPrefixConfig());
    }


    int64_t GenomeIndexPrefix::encodeKMer(const std::string_view &pattern, const int mer_length) {
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

    void GenomeIndexPrefix::setConfig(const GenomeIndexPrefixConfig &cfg) {
        config = cfg;
        MER_LENGTH = config.kMerSize;
        MER_NUM = (1 << (MER_LENGTH * 2)) + 1;
    }

    void GenomeIndexPrefix::build() {
        if (config.twoDirections) {
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
            genome->sequence_ = genome->sequence_ + revSeq + "#########################";
            suffixArray.build(genome->sequence_);
        } else {
            genomeLength = genome->sequence_.length();
            suffixArray.build(genome->sequence_);
        }
        std::cout << "finished building suffix array\n";
        patternMerMap_.resize(MER_NUM);
        longestCommonPrefix_.resize(suffixArray.length(), 0);
        buildLCP();
        std::cout << "finished building LCP\n";
        buildKMerIndex();
        fillKMerIndex();
        std::cout << "finished building kMerIndex\n";
        appearance_flag.clear();
        buildAlternativeKMerMap();
        std::cout << "finished building index\n";
    }

    void GenomeIndexPrefix::buildLCP() {
        PackedArray &lcp_ = suffixArray.lcp_;
        longestCommonPrefix_.resize(suffixArray.length(), 0);
        for (int64_t i = 0; i < suffixArray.length(); ++i) {
            longestCommonPrefix_[i] = lcp_.get(i);
        }
        lcp_.clear();// no further use
    }

    inline void GenomeIndexPrefix::setFlag(int32_t hash, int32_t length) {
        if (length == 1) {
            appearance_flag[length][0] |= (1 << hash);
        } else {
            int32_t pos = hash >> 3;
            int32_t bit = hash & 0x7;
            appearance_flag[length][pos] |= (1 << bit);
        }
    }

    inline bool GenomeIndexPrefix::getFlag(int32_t hash, int32_t length) {
        if (length == 1) {
            return (appearance_flag[length][0] >> hash) & 1;
        } else {
            int32_t pos = hash >> 3;
            int32_t bit = hash & 0x7;
            return (appearance_flag[length][pos] >> bit) & 1;
        }
    }

    void GenomeIndexPrefix::scanGenome() {
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

    void GenomeIndexPrefix::buildKMerIndex() {
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

    void GenomeIndexPrefix::fillKMerIndex() {


        const std::string &text_ = genome->sequence_;
        const PackedArray &sa_ = suffixArray.sa_;
        for (size_t i = 0; i < sa_.length(); ++i) {
            int64_t hash = encodeKMer(text_.substr(sa_.get(i), MER_LENGTH), MER_LENGTH);
            if (hash == -1) continue;
            if (patternMerMap_[hash].leftSAIndex == -1 || patternMerMap_[hash].leftSAIndex > i)
                patternMerMap_[hash].leftSAIndex = i;

        }

        int64_t last = sa_.length();
        patternMerMap_[MER_NUM - 1].leftSAIndex = last;
        for (int64_t i = MER_NUM - 2; i >= 0; --i) {
            if (patternMerMap_[i].length < MER_LENGTH) {
                patternMerMap_[i].leftSAIndex = last;
            } else {
                last = patternMerMap_[i].leftSAIndex;
                assert(last != -1);
            }
        }
    }

    void GenomeIndexPrefix::buildExtendedKMerMap() {
        for (size_t i = 0; i < MER_NUM - 1; ++i) {
            if (patternMerMap_[i + 1].leftSAIndex - patternMerMap_[i].leftSAIndex < config.minExtendRep) {
                continue; // No need to form an extended k-mer because the number of repetitions is too small
            }
            int64_t nextLayerIndex = buildExtendedKMerMapSingleIndex(
                    MER_LENGTH,
                    patternMerMap_[i].leftSAIndex,
                    patternMerMap_[i + 1].leftSAIndex);
            extendedKMerIndex_[i] = nextLayerIndex;
        }
    }

    //Note: [leftSAIndex, rightSAIndex)
    //length: matched length on previous layer
    int64_t
    GenomeIndexPrefix::buildExtendedKMerMapSingleIndex(int32_t length, int64_t leftSAIndex, int64_t rightSAIndex,
                                                       int depth) {
        const std::string &text_ = genome->sequence_;
        const PackedArray &sa_ = suffixArray.sa_;
        if (rightSAIndex - leftSAIndex < config.minExtendRep || depth >= config.maxLayer) {
            return -1; // No need to form an extended k-mer because the number of repetitions is too small
        }
        int64_t newIndex = extendedKMerNum;
        extendKMerMap newMap;
        newMap.extendedKMerMap.resize(config.extendIndexSize);
        newMap.extendedKMerIndex.resize(config.extendIndexSize, -1);
        std::vector<std::vector<char>> localAppearanceFlag;
        localAppearanceFlag.resize(config.extendLength);
        for (int i = 1; i < config.extendLength; ++i) {
            localAppearanceFlag[i].resize(1 << (std::max(2 * i - 3, 1)), 0);
        }
        for (int64_t i = leftSAIndex; i < rightSAIndex; ++i) {
            int64_t hash = encodeKMer(text_.substr(sa_.get(i) + length, config.extendLength), config.extendLength);
            if (hash == -1) continue; // Invalid k-mer, skip
            if (newMap.extendedKMerMap[hash].leftSAIndex == -1 || newMap.extendedKMerMap[hash].leftSAIndex > i)
                newMap.extendedKMerMap[hash].leftSAIndex = i;
            for (int j = config.extendLength - 1; j > 0; --j) {
                hash >>= 2; // remove the last two bits
                localAppearanceFlag[j][hash >> 3] |= (1 << (hash & 0x7));
            }

        }


        int64_t last = rightSAIndex;
        newMap.extendedKMerMap[config.extendIndexSize - 1].leftSAIndex = last;
        for (int64_t i = config.extendIndexSize - 2; i >= 0; --i) {
            if (newMap.extendedKMerMap[i].leftSAIndex == -1) {
                newMap.extendedKMerMap[i].leftSAIndex = last;
                bool exist = false;
                int64_t l = config.extendLength - 1;
                int64_t hash = i;
                while (l > 0) {
                    hash = i >> (2 * l);
                    exist = (localAppearanceFlag[l][hash >> 3] >> (hash & 0x7)) & 1;
                    if (exist) {
                        newMap.extendedKMerMap[i].length = l;
                        break;
                    }
                    --l;
                }

            } else {
                last = newMap.extendedKMerMap[i].leftSAIndex;
                newMap.extendedKMerMap[i].length = config.extendLength;
            }
        }
        //store the new extended k-mer map
        //extendedKMerMap_.emplace_back(newMap);
        for (int64_t i = 0; i < config.extendIndexSize; ++i) {
            extendedKMerMap_.push_back(newMap.extendedKMerMap[i]);
            nextLayerIndex_.push_back(-1);
        }
        ++extendedKMerNum;
        //recursively build extended k-mer map for each extended k-mer index with enough repetitions
        for (int64_t i = 0; i < config.extendIndexSize - 1; ++i) {
            if (newMap.extendedKMerMap[i + 1].leftSAIndex - newMap.extendedKMerMap[i].leftSAIndex <
                config.minExtendRep) {
                continue; // No need to form an extended k-mer because the number of repetitions is too small
            }
            int64_t nextLayerIndex = buildExtendedKMerMapSingleIndex(
                    length + config.extendLength,
                    newMap.extendedKMerMap[i].leftSAIndex,
                    newMap.extendedKMerMap[i + 1].leftSAIndex, depth + 1);
            nextLayerIndex_[newIndex * config.extendIndexSize + i] = nextLayerIndex;
        }


        return newIndex;
    }



    MatchStats GenomeIndexPrefix::find(const std::string_view &pattern) const {
        // pattern MUST include ONLY ACGT , no N at all

        // if the MMP of the pattern is shorter than k, still return the max length and the range
        // for example, pattern ACGTTGGA, can only match 4, then return 4,leftSAIndex[ACGTAAA], leftSAIndex[ACTAAAA]-1,
        if (pattern.length() < MER_LENGTH) {
            //only occurs on the last few bps of read, just ignore them
            return {0, 0, 0, 0};
        }
        // now, find the max matching prefix through the k-mer map

        int64_t hash = encodeKMer(pattern.substr(0, MER_LENGTH), MER_LENGTH);
        int l = patternMerMap_[hash].length;
        if (l < MER_LENGTH || pattern.length() == MER_LENGTH) {
            int64_t leftHash = hash >> ((MER_LENGTH - l) * 2);
            leftHash = leftHash << ((MER_LENGTH - l) * 2);
            int64_t rightHash = (leftHash | ((1 << (2 * (MER_LENGTH - l))) - 1)) + 1;
            size_t leftSAIndex = patternMerMap_[leftHash].leftSAIndex;
            size_t rightSAIndex = patternMerMap_[rightHash].leftSAIndex - 1;
            return {(size_t) l, rightSAIndex - leftSAIndex + 1, leftSAIndex, rightSAIndex};
        } else {
            int64_t nowNodeIndex = extendedKMerIndex_[hash];
            size_t lBound = patternMerMap_[hash].leftSAIndex;
            size_t rBound = patternMerMap_[hash + 1].leftSAIndex;
            int64_t remainingLength = pattern.length() - MER_LENGTH;
            while (nowNodeIndex != -1) {
                if (remainingLength > config.extendLength) {
                    // move to the next extended k-mer
                    //const extendKMerMap &nowMap = extendedKMerMap_[nowNodeIndex];

                    int64_t extendHash = encodeKMer(pattern.substr(l, config.extendLength), config.extendLength);
                    int extendLength = extendedKMerMap_[nowNodeIndex * config.extendIndexSize + extendHash].length;
                    l += extendLength;
                    if (extendLength < config.extendLength) {
                        int64_t leftHash = (extendHash >> (2 * (config.extendLength - extendLength)))
                                << (2 * (config.extendLength - extendLength));
                        int64_t rightHash = (leftHash | ((1 << (2 * (config.extendLength - extendLength))) - 1)) + 1;
                        size_t leftSAIndex = extendedKMerMap_[nowNodeIndex * config.extendIndexSize +
                                                              leftHash].leftSAIndex;
                        size_t rightSAIndex =
                                extendedKMerMap_[nowNodeIndex * config.extendIndexSize + rightHash].leftSAIndex - 1;
                        return {(size_t) l, rightSAIndex - leftSAIndex + 1, leftSAIndex, rightSAIndex};
                    }
                    lBound = extendedKMerMap_[nowNodeIndex * config.extendIndexSize + extendHash].leftSAIndex;
                    rBound = extendedKMerMap_[nowNodeIndex * config.extendIndexSize + extendHash + 1].leftSAIndex;
                    nowNodeIndex = nextLayerIndex_[nowNodeIndex * config.extendIndexSize + extendHash];
                    remainingLength -= config.extendLength;
                } else {
                    //return the MMP in the current node
                    size_t leftSAIndex = lBound, rightSAIndex = rBound - 1;
                    int64_t length = remainingLength;

                    if (length == 0) return {(size_t) l, rightSAIndex - leftSAIndex + 1, leftSAIndex, rightSAIndex};
                    //int64_t extendHash = getHashFromTotalHash(totalHash, config.extendLength, l, totalHashLength);
                    int64_t extendHash = encodeKMer(pattern.substr(l, length), length);
                    do {
                        int64_t leftHash = extendHash << ((config.extendLength - length) * 2);
                        int64_t rightHash = (extendHash + 1) << ((config.extendLength - length) * 2);
                        leftSAIndex = extendedKMerMap_[nowNodeIndex * config.extendIndexSize + leftHash].leftSAIndex;
                        rightSAIndex =
                                extendedKMerMap_[nowNodeIndex * config.extendIndexSize + rightHash].leftSAIndex - 1;
                        if (rightSAIndex - leftSAIndex > 0) {
                            l += length;
                            break;
                        } else {
                            --length;
                            extendHash >>= 2; // remove the last two bits
                        }
                    } while (length > 0);
                    return {(size_t) l, rightSAIndex - leftSAIndex + 1, leftSAIndex, rightSAIndex};
                }

            }


            auto matchStats = getMatchStatsLCP(pattern, lBound, rBound, l);

            return matchStats;
        }

    }

    inline std::pair<size_t, bool>
    GenomeIndexPrefix::matchGenomeSeq(size_t matchedLength, size_t pos, const std::string_view &pattern) const {
        size_t l;
        bool greater = false;
        for (l = matchedLength; l < pattern.length(); ++l) {
            if (pattern[l] != genome->sequence_[pos + l]) {
                greater = pattern[l] < genome->sequence_[pos + l];
                break;
            }
        }
        return {l, greater};
    }

    // Note: [lBound, rBound)
    inline int64_t
    GenomeIndexPrefix::findSameMatchingLength(size_t targetLength, int64_t lBound, int64_t rBound,
                                              size_t matchedLength, const std::string_view &pattern) const {
        if (lBound >= rBound) return 0;
        if (rBound - lBound == 1) {
            size_t base = suffixArray.sa_.get(lBound);
            auto [length, greater] = matchGenomeSeq(matchedLength, base, pattern);
            if (length == targetLength) {
                return 1;
            } else {
                return 0;
            }
        } else {
            int64_t mid = (lBound + rBound) / 2;
            int64_t midPos = suffixArray.sa_.get(mid);
            auto [midLength, midGreater] = matchGenomeSeq(matchedLength, midPos, pattern);
            if (midGreater) {
                if (midLength < targetLength) {
                    return findSameMatchingLength(targetLength, lBound, mid, matchedLength, pattern);
                } else {
                    return mid + 1 - lBound +
                           findSameMatchingLength(targetLength, mid + 1, rBound, matchedLength, pattern);
                }
            } else {
                if (midLength < targetLength) {
                    return findSameMatchingLength(targetLength, mid + 1, rBound, matchedLength, pattern);
                } else {
                    return rBound - mid +
                           findSameMatchingLength(targetLength, lBound, mid, matchedLength, pattern);
                }
            }

        }
    }

    MatchStats GenomeIndexPrefix::getMatchStatsLCP(const std::string_view &pattern, size_t rangeLeft, size_t rangeRight,
                                                   size_t matchedLength) const {

        //first ,get the longest match
        if (rangeLeft >= rangeRight) {
            return {0, 0, 0, 0};
        }
        rangeRight--;
        if (rangeLeft == rangeRight) {
            auto [longestLength, greater] = matchGenomeSeq(matchedLength, suffixArray.sa_.get(rangeLeft), pattern);
            return {longestLength, 1, rangeLeft, rangeLeft};
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
                return {midLength, rightPos - leftPos + 1, leftPos, rightPos};
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
                if (cnt > 10000) return {longestLength,0,0,0};
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
                if (cnt > 10000) return {longestLength,0,0,0};
                rightCommonPrefix = longestCommonPrefix_[rightPos + 1];
            }
        } else rightPos = rangeRight - 1;


        return {longestLength, rightPos - leftPos + 1, leftPos, rightPos};
    }


    std::vector<Align> GenomeIndexPrefix::alignRead(const std::string_view &readSeq, const Split &split) const {
        // align the read to the genome, return the alignments
        std::vector<Align> aligns;
        aligns.reserve(readSeq.length() / 10);
        int64_t length;
        int iStart = 1 + split.length / 50;
        bool fullMatch = false;
        for (int i = 0; i < iStart; ++i) {
            int64_t nowMappedLength = i * 50;
            do {
                if (split.length - nowMappedLength < 5) break;
                auto matchStats = findAlternative(readSeq.substr(split.readStart + nowMappedLength, split.length - nowMappedLength));
                length = matchStats.maxLength;
                if (length > 0) {
                    aligns.emplace_back(
                            Align{
                                    .readStart = split.readStart + nowMappedLength,
                                    .length = matchStats.maxLength,
                                    .leftSAIndex = matchStats.LongestMatchLeftPosInSA,
                                    .rightSAIndex = matchStats.LongestMatchRightPosInSA,
                                    .rep = matchStats.count,
                                    .direction = split.direction,
                                    .isAnchor = false
                            }
                    );
                }
                if (length == split.length){
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


        return aligns;
    }

    void GenomeIndexPrefix::write(const std::string &outDir) const {
        if (!std::filesystem::exists(outDir)) std::filesystem::create_directories(outDir);
        std::ofstream logFile(outDir + "/indexInf.txt");
        std::ofstream genomeFile(outDir + "/genome.bin", std::ios::binary | std::ios::out);
        std::ofstream saFile(outDir + "/SA.bin", std::ios::binary | std::ios::out);
        std::ofstream saIndexFile(outDir + "/SAIndex.bin", std::ios::binary | std::ios::out);
        //write genome sequence
        logFile << genomeLength << "\n";
        logFile << genome->sequence_.length() << "\n";
        genomeFile.write(genome->sequence_.data(), genome->sequence_.length() * sizeof(char));
        genomeFile.close();
        logFile << genome->chromosomes_.size() << '\n';
        for (auto &chromosome: genome->chromosomes_) {
            logFile << chromosome.name << " " << chromosome.start << ' ' << chromosome.length << "\n";
        }
        //write suffix array
        logFile << suffixArray.length() << " " << suffixArray.fullLength_ << "\n";
        suffixArray.sa_.writeToFile(saFile, logFile);

        //write suffix array index
        logFile << config.kMerSize << " " << config.extendLength << " "
                << config.maxLayer << " " << config.minExtendRep << " " << config.twoDirections << "\n";

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

    void GenomeIndexPrefix::load(const std::string &inDir) {
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
        genome->sequence_.resize(genomeSeqLength);
        genomeFile.read(genome->sequence_.data(), genomeSeqLength * sizeof(char));
        genomeFile.close();
        genome->chromosomes_.clear();
        std::string chrName;
        int64_t chrStart, chrLength, chrNum;
        logFile >> chrNum;
        for (int64_t i = 0; i < chrNum; ++i) {
            logFile >> chrName >> chrStart >> chrLength;
            genome->chromosomes_.emplace_back(Genome::Chromosome{chrName, chrStart, chrLength});
            genome->chromosomeMap_[chrStart] = i;
        }
        int64_t saLength, saFullLength;
        logFile >> saLength >> suffixArray.fullLength_;
        suffixArray.sa_.loadFromFile(saFile, logFile);


        logFile >> config.kMerSize >> config.extendLength
                >> config.maxLayer >> config.minExtendRep >> config.twoDirections;
        config.extendIndexSize = (1 << (2 * config.extendLength)) + 1;
        setConfig(GenomeIndexPrefixConfig{config.kMerSize, config.extendLength, config.extendIndexSize,
                                          config.minExtendRep, config.maxLayer,
                                          config.twoDirections});

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

        logFile.close();
        saFile.close();
        saIndexFile.close();
    }


    void GenomeIndexPrefix::buildAlternativeKMerMap() {
        extendIndexHash.resize(MER_NUM * config.extendAlternativeByte/4 , 0);
        int64_t nowLeftSAIndex = patternMerMap_[0].leftSAIndex;
        int64_t nextLeftSAIndex;
        for (size_t i = 0; i < MER_NUM - 1; ++i) {
            nextLeftSAIndex = patternMerMap_[i + 1].leftSAIndex;
            int64_t num = nextLeftSAIndex - nowLeftSAIndex;
            if (num > 16) {
                // 4 Bytes per index
                int indNum = config.extendAlternativeByte / 4;
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
                    extendIndexHash[i * indNum + j] = (uint32_t)hash;
                }

            } else {

            }
            nowLeftSAIndex = nextLeftSAIndex;
        }
    }

    MatchStats GenomeIndexPrefix::findAlternative(const std::string_view &pattern) const {
        if (pattern.length() < MER_LENGTH) {
            //only occurs on the last few bps of read, just ignore them
            return {0, 0, 0, 0};
        }
        // now, find the max matching prefix through the k-mer map

        int64_t hash = encodeKMer(pattern.substr(0, MER_LENGTH), MER_LENGTH);
        int l = patternMerMap_[hash].length;
        if (l < MER_LENGTH || pattern.length() == MER_LENGTH) {
            int64_t leftHash = hash >> ((MER_LENGTH - l) * 2);
            leftHash = leftHash << ((MER_LENGTH - l) * 2);
            int64_t rightHash = (leftHash | ((1 << (2 * (MER_LENGTH - l))) - 1)) + 1;
            size_t leftSAIndex = patternMerMap_[leftHash].leftSAIndex;
            size_t rightSAIndex = patternMerMap_[rightHash].leftSAIndex - 1;
            return {(size_t) l, rightSAIndex - leftSAIndex + 1, leftSAIndex, rightSAIndex};
        } else {
            size_t lBound = patternMerMap_[hash].leftSAIndex;
            size_t rBound = patternMerMap_[hash + 1].leftSAIndex;
            size_t num = rBound - lBound;
            const int64_t extendIndexHashBase = hash * config.extendAlternativeByte/4;
            int64_t remainingLength = pattern.length() - MER_LENGTH;
            //todo this may  better
            if (remainingLength < 16) {
                uint32_t h = encodeKMer(pattern.substr(l, remainingLength), remainingLength);
                uint32_t hLeft = h <<  ((16 - remainingLength)*2);
                uint32_t hRight = (h+1) << ((16 - remainingLength)*2);
                size_t left = 0, right = config.extendAlternativeByte/4;
                size_t mid;
                while (right - left > 1){
                    mid = (left + right) / 2;
                    if (extendIndexHash[extendIndexHashBase + mid] < hLeft) left = mid;
                    else if (extendIndexHash[extendIndexHashBase + mid] > hRight) right = mid;
                    else {
                        // true left bound in [left,mid] ,true right bound in [mid,right]
                        size_t boundLeft = mid;
                        while (boundLeft - left > 1){
                            size_t midLeft = (left + boundLeft) / 2;
                            if (extendIndexHash[extendIndexHashBase + midLeft] < hLeft) left = midLeft;
                            else if (extendIndexHash[extendIndexHashBase + midLeft] > hLeft) boundLeft = midLeft;
                            else {
                                // exact match
                                left =  midLeft;
                                while (left > 0 && extendIndexHash[extendIndexHashBase + left] == hLeft)
                                    --left;
                                break;
                            }
                        }
                        size_t boundRight = mid;
                        while (right - boundRight > 1){
                            size_t midRight = (boundRight + right) / 2;
                            if (extendIndexHash[extendIndexHashBase + midRight] < hRight) boundRight = midRight;
                            else if (extendIndexHash[extendIndexHashBase + midRight] > hRight) right = midRight;
                            else {
                                // exact match
                                right = midRight;
                                while (right < config.extendAlternativeByte/4 &&
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
                if (right != config.extendAlternativeByte/4) ++rBound;
                lBound += left * num * 4 / config.extendAlternativeByte;

                return getMatchStatsLCP(pattern, lBound, rBound, l);
            }
            if (num > 16) {
                uint32_t h = encodeKMer(pattern.substr(l, 16), 16);
                size_t left = 0, right = config.extendAlternativeByte/4;
                while (right - left > 1) {
                    size_t mid = (left + right) / 2;
                    if (extendIndexHash[extendIndexHashBase + mid] < h) left = mid;
                    else if (extendIndexHash[extendIndexHashBase + mid] > h) right = mid;
                    else {
                        // exact match
                        left = right = mid;
                        while (left > 0 && extendIndexHash[extendIndexHashBase + left] == h)
                            --left;
                        while (right < config.extendAlternativeByte/4 &&
                               extendIndexHash[extendIndexHashBase + right] == h)
                            ++right;
                        break;
                    }
                }
                rBound = lBound + right * num * 4 / config.extendAlternativeByte;
                if (right != config.extendAlternativeByte/4) ++rBound;
                lBound += left * num * 4 / config.extendAlternativeByte;

            } else {
                // no need to search
            }


            auto matchStats = getMatchStatsLCP(pattern, lBound, rBound, l);

            return matchStats;
        }
    }


}