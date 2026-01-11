#include "GenomeIndex.h"
#include "../log/ErrorRecord.h"
#include "../utilsRefactored/seqFunctions.h"


namespace RefactorProcessing {
    inline bool GenomeIndex::insertAlign(std::vector<Align> &results, const Align &a) const {


        if (results.size() >= maxAlignNum) {
            rna::WarningRecord().reportWarning("Exceeding max alignments per read, some alignments are ignored");
            return false;
        }

        int posToInsert = -1;
        int64_t readStart = a.readPos;
        int64_t length = a.length;

        // find the position to insert
        for (int i = 0; i < results.size(); ++i) {
            if (results[i].readPos < readStart) continue;
            if (results[i].readPos == readStart) {
                if (results[i].length < length) continue;
                if (results[i].length == length) {
                    if (results[i].direction != a.direction) continue;
                    else {
                        // same alignment already exists
                        return false;
                    }
                }
            }
            posToInsert = i;
            break;
        }
        results.push_back(Align());
        for (int i = results.size() - 1; i > posToInsert; --i) {
            results[i] = results[i - 1];
        }
        results[posToInsert] = a;


        return true;
    }

    void GenomeIndex::find(const Split pattern, std::vector<Align> &results) const {
        // align the read to the genome
        // std::<vector> results must be reserved


        int64_t splitLength = pattern.length;
        int64_t length;
        int64_t iStart = 1 + (splitLength / 50); // todo replace 50 by parameter
        int64_t lStart = splitLength / iStart;

        bool fullMatch = false;
        std::string_view splits[2] = {pattern.forward, pattern.reverse};
        for (int dir = 0; dir <= 1; ++dir) {
            std::string_view split = splits[dir];
            int64_t splitStart = dir == 0 ? pattern.splitStart : (pattern.readLength - pattern.splitStart - pattern.length);
            for (int i = 0; i < iStart; ++i) {
                int64_t nowMappedLength = i * lStart;

                // search MMP aligns

                length = 0;
                int64_t matchedLength = nowMappedLength;
                do {

                    if (splitLength - nowMappedLength <= 5) {
                        // avoid too short matches
                        // todo replace 5 by parameter
                        length = 0;
                        break;
                    }
                    Align matchedAlign = findMMP(split.substr(matchedLength));
                    matchedAlign.readPos = nowMappedLength + splitStart;
                    length = matchedAlign.length;
                    if (length == 0) break;
                    nowMappedLength += length;
                    insertAlign(results, matchedAlign);
                    if (length == splitLength) fullMatch = true;
                } while (length > 0);
                if (fullMatch) break;
            }
            if (fullMatch) break;
        }
    }

    inline int64_t extendHashLeft(int64_t hash, int extendLength) {
        return hash << (2*extendLength);
    }

    inline int64_t extendHashRight(int64_t hash, int extendLength) {
        return (hash << (2* extendLength) ) | ((1ULL << (2* extendLength)) - 1);
    }

    inline int64_t matchedHashLeft(int64_t hash, int remainLength) {
        // xxxAAAA
        return (hash >> (2* remainLength)) << (2* remainLength);
    }

    inline int64_t matchedHashRight(int64_t hash, int remainLength) {
        // xxxTTTT
        return hash | ((1ULL << (2* remainLength)) - 1);
    }
    Align GenomeIndex::findMMP(const std::string_view& seq) const {
        Align result;
        if (seq.length() < kMerSize_){
            //the last few bps of read
            int l = seq.length();
            int64_t hash = encodeKmer(seq,l);
            if (hash == -1) return Align();
            int64_t leftHash = extendHashLeft(hash, kMerSize_ - l);
            int64_t rightHash = extendHashRight(hash, kMerSize_ - l);
            size_t leftSAIndex = patternMerMap_.get(leftHash, patternMerMap_.INDEX_LEFT_SA_INDEX);
            size_t rightSAIndex = patternMerMap_.get(rightHash, patternMerMap_.INDEX_LEFT_SA_INDEX) +
                                 patternMerMap_.get(rightHash, patternMerMap_.INDEX_UPPER_RANGE);
            result.length = l;
            result.leftSAIndex = leftSAIndex;
            result.rightSAIndex = rightSAIndex;
            result.rep = rightSAIndex - leftSAIndex + 1;

            return result;
        }

        //find the MMP

        int64_t hash = encodeKmer(seq.substr(0,kMerSize_),kMerSize_);
        int l = patternMerMap_.get(hash, patternMerMap_.INDEX_LENGTH);
        if (l < kMerSize_ || seq.length() == kMerSize_){
            int64_t leftHash = matchedHashLeft(hash, kMerSize_ - l);
            int64_t rightHash = matchedHashRight(hash, kMerSize_ - l);
            size_t leftSAIndex = patternMerMap_.get(leftHash, patternMerMap_.INDEX_LEFT_SA_INDEX);
            size_t rightSAIndex = patternMerMap_.get(rightHash, patternMerMap_.INDEX_LEFT_SA_INDEX) +
                                 patternMerMap_.get(rightHash  , patternMerMap_.INDEX_UPPER_RANGE);
            if (rightSAIndex - leftSAIndex + 1 > 10000){
                return Align();
            } else {
                int64_t nowRightSAIndex = rightSAIndex;
                while (longestCommonPrefix_[nowRightSAIndex + 1] >= l) {
                    ++nowRightSAIndex;
                    if (nowRightSAIndex >= suffixArray_.length_) break;
                }
                rightSAIndex = nowRightSAIndex;
                // no need to search left bound
                result.length = l;
                result.leftSAIndex = leftSAIndex;
                result.rightSAIndex = rightSAIndex;
                result.rep = rightSAIndex - leftSAIndex + 1;
                return result;
            }
        }

        // common situation
        size_t lBound = patternMerMap_.get(hash,patternMerMap_.INDEX_LEFT_SA_INDEX);
        size_t num = patternMerMap_.get(hash,patternMerMap_.INDEX_UPPER_RANGE) + 1;
        size_t rBound = lBound + num;

        //no need to search in extend hash table







    }

    void GenomeIndex::findMMP_GetRange(const std::string &seq, int64_t rangeLeft, int64_t rangeRight,
                                       size_t matchedLength) const {

    }

    inline std::pair<size_t, bool>
    GenomeIndex::matchGenomeSeq(const std::string &pattern, size_t matchedLength, size_t pos) const {

    }
}
