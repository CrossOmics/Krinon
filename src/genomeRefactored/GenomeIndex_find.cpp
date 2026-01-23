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
        if (num <= extendHashTableNum_) return findMMP_GetRange(seq,lBound,rBound,l);
        // start pos of current hash's table
        const int64_t extendIndexHashBase = hash * extendHashTableNum_;
        int64_t remainLength = seq.length() - kMerSize_;
        int left = 0,right = extendHashTableNum_;
        if (remainLength < 16) {
            int64_t extendHash = encodeKmer(seq.substr(l,remainLength),remainLength);
            int64_t lHash = extendHashLeft(extendHash,16 - remainLength);
            int64_t rHash = extendHashRight(extendHash,16 - remainLength);

            //find the range between lHash and rHash


            int mid;

            while (right - left > 1){
                mid = (left + right) / 2;
                if (extendedIndexHash_[extendIndexHashBase + mid] < lHash) left = mid;
                else if (extendedIndexHash_[extendIndexHashBase + mid] > rHash) right = mid;
                else {
                    int lInd = mid;
                    // binary search
                    while (lInd - left > 1){
                        int midLeft = (lInd + left) / 2;
                        if (extendedIndexHash_[extendIndexHashBase + midLeft] < lHash) left = midLeft;
                        else if (extendedIndexHash_[extendIndexHashBase +midLeft] > lHash) lInd = midLeft;
                        else {
                            left = midLeft;
                            while (left > 0 && extendedIndexHash_[extendIndexHashBase + left] == lHash) --left;
                            break;
                        }
                    }

                    int rInd = mid;
                    while (right - rInd > 1){
                        int midRight = (rInd + right) / 2;
                        if (extendedIndexHash_[extendIndexHashBase +midRight] <rHash) rInd = mid;
                        else if (extendedIndexHash_[extendIndexHashBase +midRight] > rHash) right = mid;
                        else {
                            right = midRight;
                            while (right < extendHashTableNum_ && extendedIndexHash_[extendIndexHashBase +right] == rHash) ++right;
                            break;
                        }
                    }
                    break;
                }
            }

        }else{
            int64_t extendHash = encodeKmer(seq.substr(l,16),16);
            while (right - left > 1) {
                int mid = (left + right) /2;
                if (extendedIndexHash_[extendIndexHashBase + mid] < extendHash) left = mid;
                else if (extendedIndexHash_[extendIndexHashBase +mid] >extendHash) right = mid;
                else {
                    left = right = mid;
                    while (left > 0 && extendedIndexHash_[extendIndexHashBase + left] == extendHash) --left;
                    while (right < extendHashTableNum_ && extendedIndexHash_[extendIndexHashBase + right] == extendHash) ++right;
                    break;
                }
            }
        }

        rBound = lBound + (right * num + extendHashTableNum_ - 1 )/ extendHashTableNum_; // ceil
        lBound += (left * num) /extendHashTableNum_;

        return findMMP_GetRange(seq,lBound,rBound,l);





    }

    Align GenomeIndex::findMMP_GetRange(const std::string_view &seq, int64_t rangeLeft, int64_t rangeRight,
                                       size_t matchedLength) const {
        Align result;
        if (rangeLeft >= rangeRight) {
            result.length = 0;
            result.rep = 0;
            result.leftSAIndex = 0;
            result.rightSAIndex = 0;
            return result;
        }
        rangeRight--;
        if (rangeLeft == rangeRight) {
            auto [longestLength, greater] = matchGenomeSeq(seq,matchedLength,suffixArray_[rangeLeft]);
            result.length = longestLength;
            result.rep = 1;
            result.leftSAIndex = rangeLeft;
            result.rightSAIndex = rangeLeft;
            return result;
        }
        size_t leftMaxLength = matchedLength;
        size_t rightMaxLength = matchedLength;
        size_t cnt = 0;
        while (rangeRight - rangeLeft > 1) {
            size_t mid = (rangeLeft + rangeRight) / 2;
            size_t midPos = suffixArray_[mid];
            auto [midLength, midGreater] = matchGenomeSeq(seq,matchedLength,midPos);
            if (midLength == seq.length()) {
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
                result.length = midLength;
                result.rep = rightPos - leftPos + 1;
                result.leftSAIndex = leftPos;
                result.rightSAIndex = rightPos;
                return result;
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
        auto [lLength, leftGreater] = matchGenomeSeq(seq,leftMaxLength,suffixArray_[rangeLeft]);
        auto [rLength, rightGreater] = matchGenomeSeq(seq,rightMaxLength,suffixArray_[rangeRight]);
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
                if (cnt > 10000) {
                    result.length = longestLength;
                    result.rep = 0;
                    result.leftSAIndex = 0;
                    result.rightSAIndex = 0;
                    return result;

                }
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
                if (cnt > 10000) {
                    result.length = longestLength;
                    result.rep = 0;
                    result.leftSAIndex = 0;
                    result.rightSAIndex = 0;
                    return result;

                };
                rightCommonPrefix = longestCommonPrefix_[rightPos + 1];
            }
        } else rightPos = rangeRight - 1;


        result.length = longestLength;
        result.rep = rightPos - leftPos + 1;
        result.leftSAIndex = leftPos;
        result.rightSAIndex = rightPos;
        return result;
    }

    inline std::pair<size_t, bool>
    GenomeIndex::matchGenomeSeq(const std::string_view &pattern, size_t matchedLength, size_t pos) const {
        size_t l;
        bool greater = false;
        for (l = matchedLength; l < pattern.length(); ++l) {
            if (pattern[l] != genome_.sequence_[pos + l]) {
                if (charToIndex(genome_.sequence_[pos + l]) < 0)
                    greater = true;// 'N' and '#' are regarded as the largest
                else greater = pattern[l] < genome_.sequence_[pos + l];
                break;
            }
        }
        return {l, greater};
    }
}
