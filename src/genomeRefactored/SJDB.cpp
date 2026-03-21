#include <fstream>
#include <cstring>
#include "SJDB.h"
#include "GenomeIndex.h"
#include "omp.h"
#include <string_view>

#define GTAG 1
#define CTAC 2
#define GCAG 3
#define CTGC 4
#define ATAC 5
#define GTAT 6
#define NON_CANONICAL 0
namespace RefactorProcessing {
    void SJDB::setParam(const Parameters &P) {
        sjdbGTFfeatureExon = P.sjdbGTFfeatureExon;
        sjdbGTFTagExonParentTranscriptId = P.sjdbGTFTagExonParentTranscriptId;
        sjdbGTFChrPrefix = P.sjdbGTFChrPrefix;
        sjdbGTFTagExonParentGene = P.sjdbGTFTagExonParentGene;
        sjdbGTFTagExonParentGeneName = P.sjdbGTFTagExonParentGeneName;
        sjdbGTFTagExonParentGeneType = P.sjdbGTFTagExonParentGeneType;
        sjdbOverhang = P.sjdbOverhang;
        sjdbLength = P.sjdbLength;
        limitSjdbInsertN = P.limitSjdbInsertN;

        logFile_ = std::ofstream(P.outLogFile, std::ios::app);
    }

    void SJDB::loadGTF(const std::string &gtfFile, const RefactorProcessing::Genome &genome) {
        std::ifstream sjdbStreamIn(gtfFile);
        char sjdbReadBuffer[65536];
        sjdbStreamIn.rdbuf()->pubsetbuf(sjdbReadBuffer, sizeof sjdbReadBuffer);
        std::map<std::string, int64_t> transcriptIDNumber, geneIDNumber;
        if (sjdbStreamIn.fail()) {
            //todo report error
            std::cerr << "Error: cannot open GTF file " << gtfFile << "\n";
            exit(1);
        }

        if (chromosomeNameToIndex_.empty()) {
            for (size_t i = 0; i < genome.chromosomes_.size(); ++i) {
                chromosomeNameToIndex_[genome.chromosomes_[i].name] = i;
            }
        }

        int exonN = 0;
        while (sjdbStreamIn.good()) {//count the number of exons

            std::string oneLine, chr1, ddd2, featureType;
            getline(sjdbStreamIn, oneLine);
            std::istringstream oneLineStream(oneLine);

            oneLineStream >> chr1 >> ddd2 >> featureType;
            if (chr1.substr(0, 1) != "#" && featureType == sjdbGTFfeatureExon) {
                exonN++;
            };
        };
        if (exonN == 0) return;

        sjdbStreamIn.clear();
        sjdbStreamIn.seekg(0, std::ios::beg);

        exonLoci_.resize(exonN);
        exonN = 0;
        while (sjdbStreamIn.good()) {//read the exons
            std::string oneLine, chr1, ddd2, featureType;
            getline(sjdbStreamIn, oneLine);
            std::istringstream oneLineStream(oneLine);

            oneLineStream >> chr1 >> ddd2 >> featureType;
            if (chr1.substr(0, 1) != "#" && featureType == sjdbGTFfeatureExon) {
                //exon line
                if (!sjdbGTFChrPrefix.empty()) chr1 = sjdbGTFChrPrefix + chr1;

                if (chromosomeNameToIndex_.find(chr1) == chromosomeNameToIndex_.end()) {
                    std::cerr << "Warning: chromosome " << chr1 << " not found in genome, skipping this exon\n";
                    continue;
                }

                size_t exonStart, exonEnd;
                char str1;
                oneLineStream >> exonStart >> exonEnd >> ddd2 >> str1 >> ddd2;
                if (exonEnd > genome.chromosomes_[chromosomeNameToIndex_[chr1]].length) {
                    std::cerr << "Warning: exon end " << exonEnd << " exceeds chromosome " << chr1
                              << " length, skipping this exon\n";
                    continue;
                }
                std::string nowLine;
                getline(oneLineStream, nowLine);
                std::replace(nowLine.begin(), nowLine.end(), ';', ' ');
                std::replace(nowLine.begin(), nowLine.end(), '=', ' ');
                std::replace(nowLine.begin(), nowLine.end(), '\t', ' ');
                std::replace(nowLine.begin(), nowLine.end(), '\"', ' ');

                //trID, gID, gName, gBiotype
                std::vector<std::vector<std::string>> exAttrNames(
                        {{sjdbGTFTagExonParentTranscriptId}, {sjdbGTFTagExonParentGene},
                         sjdbGTFTagExonParentGeneName, sjdbGTFTagExonParentGeneType});

                std::vector<std::string> exAttr; //trID, gID, gName, gBiotype
                exAttr.resize(exAttrNames.size());

                for (size_t ii = 0; ii < exAttrNames.size(); ii++) {
                    for (auto &attr1: exAttrNames[ii]) {//scan through possible names
                        size_t pos1 = nowLine.find(" " + attr1 + " "); //attribute name is separated by spaces
                        if (pos1 != std::string::npos)
                            pos1 = nowLine.find_first_not_of(' ', pos1 + attr1.size() + 1);
                        if (pos1 != std::string::npos) {
                            exAttr[ii] = nowLine.substr(pos1, nowLine.find_first_of(' ', pos1) - pos1);
                        }
                    }
                }

                if (exAttr[0].empty()) {//no transcript ID
                    logFile_ << "WARNING: while processing pGe.sjdbGTFfile=" << gtfFile
                             << ": no transcript_id for line:\n";
                    logFile_ << oneLine << "\n";
                    logFile_.flush();
                    exAttr[0] = "tr_" + chr1 + "_" + std::to_string(exonStart) + "_" + std::to_string(exonEnd) + "_" +
                                std::to_string(exonN); //unique name for the transcript
                }

                if (exAttr[1].empty()) {//no gene ID
                    logFile_ << "WARNING: while processing pGe.sjdbGTFfile=" << gtfFile << ": no gene_id for line:\n";
                    logFile_ << oneLine << "\n";
                    logFile_.flush();
                    exAttr[1] = "MissingGeneID";
                }

                if (exAttr[2].empty()) {//no gene name
                    exAttr[2] = exAttr[1];
                }

                if (exAttr[3].empty()) {//no gene type
                    exAttr[3] = "MissingGeneType";
                }

                transcriptIDNumber.insert(std::pair<std::string, size_t>(exAttr[0], transcriptIDNumber.size()));
                if (transcriptID_.size() < transcriptIDNumber.size()) {
                    transcriptID_.push_back(exAttr[0]);
                    if (str1 == '+') transcriptStrand_.push_back(1);
                    else if (str1 == '-') transcriptStrand_.push_back(2);
                    else transcriptStrand_.push_back(0);
                }

                geneIDNumber.insert(std::pair<std::string, size_t>(exAttr[1], geneIDNumber.size()));
                if (geneID_.size() < geneIDNumber.size()) {//new gene is added
                    geneID_.push_back(exAttr[1]);
                    geneAttr_.emplace_back(exAttr[2], exAttr[3]);
                }

                exonLoci_[exonN] = {transcriptIDNumber[exAttr[0]], static_cast<int64_t>(exonStart +
                                                                                        genome.chromosomes_[chromosomeNameToIndex_[chr1]].start -
                                                                                        1),
                                    static_cast<int64_t>(exonEnd +
                                                         genome.chromosomes_[chromosomeNameToIndex_[chr1]].start -
                                                         1), geneIDNumber[exAttr[1]]};
                exonN++;

            }
        }

        if (exonN == 0) std::cerr << "Warning: no exons found in GTF file " << gtfFile << "\n";
        exonLoci_.resize(exonN);
    }

    int SJDB::fillSjdbLoci(const std::string &dirOut, Genome &genome) {
        // from exonLoci_ to sjdbLoci_
        // return the number of junctions
        std::sort(exonLoci_.begin(), exonLoci_.end());
        size_t exonNum = exonLoci_.size();
        if (exonNum == 0) return 0;
        //todo output exon-gene info

        std::vector<exonTrLoci> extrLoci;
        extrLoci.resize(exonNum);

        size_t trExon = 0;
        for (size_t i = 0; i <= exonNum; ++i) {
            if (i == exonNum || exonLoci_[i].trID != exonLoci_[trExon].trID) {
                //process exons from trExon to i-1
                for (size_t j = trExon; j < i; ++j) {
                    extrLoci[j].trEnd = exonLoci_[i - 1].end;
                }
                if (i == exonNum) break;
                trExon = i;
            }
            extrLoci[i].trStart = exonLoci_[trExon].start;
            extrLoci[i].trID = exonLoci_[i].trID;
            extrLoci[i].exonStart = exonLoci_[i].start;
            extrLoci[i].exonEnd = exonLoci_[i].end;
            extrLoci[i].geneId = exonLoci_[i].geneId;
        }
        std::sort(extrLoci.begin(), extrLoci.end());

        //todo output exon transcript info

        std::vector<sjStride> sjLoci;
        sjLoci.reserve(exonNum);
        int64_t trID = exonLoci_[0].trID;

        for (size_t i = 1; i < exonNum; ++i) {
            if (trID == exonLoci_[i].trID) {
                size_t chr = genome.getPosChrIndex(exonLoci_[i].start);
                if (exonLoci_[i].start <= exonLoci_[i - 1].end + 1) {
                    //touching or overlapping exons
                    continue;
                }
                sjLoci.emplace_back(exonLoci_[i - 1].end + 1, exonLoci_[i].start - 1, transcriptStrand_[trID],
                                    exonLoci_[i].geneId);
            } else {
                trID = exonLoci_[i].trID;
            }
        }

        std::sort(sjLoci.begin(), sjLoci.end());
        char strandChar[3] = {'.', '+', '-'};
        size_t sjNum = sjLoci.size();
        size_t prevSjNum = sjdbLoci_.size();
        sjdbLoci_.reserve(sjNum + prevSjNum);
        for (size_t i = 0; i < sjNum; ++i) {
            if (i == 0 || sjLoci[i] != sjLoci[i - 1]) {
                size_t chr = genome.getPosChrIndex(sjLoci[i].sjStart);
                sjdbLoci_.emplace_back(genome.chromosomes_[chr].name,
                                       sjLoci[i].sjStart - genome.chromosomes_[chr].start + 1,
                                       sjLoci[i].sjEnd - genome.chromosomes_[chr].start + 1,
                                       strandChar[sjLoci[i].sjStrand], 20);
                sjdbLoci_.back().gene.insert(sjLoci[i].geneId);
            } else {
                sjdbLoci_.back().gene.insert(sjLoci[i].geneId);
            }
        }

        // todo output sjdbList

        return (int) (sjdbLoci_.size() - prevSjNum);
    }


    inline char complement(char base) {
        switch (base) {
            case 'A':
                return 'T';
            case 'T':
                return 'A';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            default:
                return base; // For 'N' or any other character
        }
    }


    struct insertRecord {
        int64_t pos;
        size_t sjPos;
    };

    class insertRecordComparator {
    private:
        const char *sjdbSeq_;
    public:
        explicit insertRecordComparator(const char *sjdbSeq) : sjdbSeq_(sjdbSeq) {}

        bool operator()(const insertRecord &a, const insertRecord &b) const {
            if (a.pos != b.pos) return a.pos < b.pos;
            const char *seqA = sjdbSeq_ + a.sjPos;
            const char *seqB = sjdbSeq_ + b.sjPos;
            int i = 0;
            while (true) {
                if (seqA[i] != seqB[i]) {
                    if (seqA[i] == 'N' || seqA[i] == '#') return false;
                    if (seqB[i] == 'N' || seqB[i] == '#') return true;
                    return seqA[i] < seqB[i];
                }
                if (seqA[i] == '#') return a.sjPos < b.sjPos;
                i++;
            }
        }
    };


    void GenomeIndex::modify(SJDB &sjdb) {
        genome_.modifyGenome(sjdb);
        //change SA and index
        size_t sjSeqLength = sjdb.sjdbLength * genome_.sjdbNum_;
        genome_.sjdbSeqLength_ = sjSeqLength;
        //notice that there will be 2*SJDB_PADDING_LENGTH '#'s between the two strands
        for (size_t i = 0; i < sjSeqLength; ++i) {
            sjdb.sjdbSeq_[2 * sjSeqLength - 1 - i] = complement(sjdb.sjdbSeq_[i]);
        }
        sjdb.sjdbSeq_ += std::string(SJDB_PADDING_LENGTH,'#'); //padding to avoid overflow

        genome_.sequence_ += sjdb.sjdbSeq_;

        //modify SA, find insert positions i in SA (i-1, new , i)
        std::vector<insertRecord> insertPos;
        insertPos.resize(2*sjSeqLength);

        std::string_view sv(sjdb.sjdbSeq_);
        omp_set_num_threads(8);
#pragma omp parallel for
        for (size_t i = 0; i < 2 * genome_.sjdbNum_; ++i) {
            for (size_t sjStart = 0; sjStart < sjdb.sjdbLength; ++sjStart) {
                size_t ind = i * sjdb.sjdbLength + sjStart;
                if (sv[i * sjdb.sjdbLength + sjStart] == '#' ||
                    sv[i * sjdb.sjdbLength + sjStart] == 'N') {
                    insertPos[ind].pos = -1;
                } else {

                    insertPos[ind].pos = sjFindInsertPosition(sv.substr(i * sjdb.sjdbLength + sjStart));
                    insertPos[ind].sjPos = i * sjdb.sjdbLength + sjStart;
                }

            }

        }
        int trueIndNum = 0;
        for (size_t i = 0; i < 2 * sjSeqLength; ++i) {
            /**
             * TODO: Again, comparison between signed and unsigned ...
             *       This one in particular looks unintended ...
             */
            if (insertPos[i].pos != -1){
                insertPos[trueIndNum] = insertPos[i];
                trueIndNum++;
            }
        }

        std::sort(insertPos.begin(),insertPos.begin() + trueIndNum,insertRecordComparator(sjdb.sjdbSeq_.c_str()));

        /**
         * TODO: What's going on here??
         *       Are we really sure that we want to assign a negative value to
         *       an unsigned integer?
         */
         /**
          * fixed by changing the type of pos to int64_t
          * mark the end of valid insert positions with -999
          */
        insertPos[trueIndNum].pos = -999;

        int64_t nowInsertSjIndex = 0;
        int64_t nowSAIndex = 0;

        SuffixArray sa;
        sa.buildFromOtherInit(suffixArray_, suffixArray_.length_ + trueIndNum);

        for (int64_t i = 0; i < (int64_t) suffixArray_.length_; ++i) {
            while (i == insertPos[nowInsertSjIndex].pos) {
                size_t sjPos = insertPos[nowInsertSjIndex].sjPos;

                sjPos += genome_.sjdbStart_;

                sa.suffixArray_.setValue(nowSAIndex,sjPos);
                nowSAIndex++;
                nowInsertSjIndex++;
            }
            size_t nowSAValue = suffixArray_[i];
            sa.suffixArray_.setValue(nowSAIndex, nowSAValue);
            nowSAIndex++;
        }

        for (; nowInsertSjIndex < trueIndNum; ++nowInsertSjIndex) {
            size_t sjPos = insertPos[nowInsertSjIndex].sjPos;

            sjPos += genome_.sjdbStart_;

            sa.suffixArray_.setValue(nowSAIndex, sjPos);
            nowSAIndex++;
        }

        suffixArray_ = std::move(sa);
        sa.suffixArray_.setDeleted();//prevent double free

        // recalculate LCP
        longestCommonPrefix_.resize(suffixArray_.length_, 0);
        PackedArray rk;
        rk.initialize(genome_.sequence_.length(), suffixArray_.wordBits_);
        for (size_t i = 0; i < genome_.sequence_.length(); ++i) rk.setValue(i, 0);
        for (size_t i = 0; i < suffixArray_.length_; ++i) rk.setValue(suffixArray_[i], i);
        size_t k = 0;
        for (size_t i = 0; i < genome_.sequence_.length(); ++i) {
            size_t j = rk[i];
            if (j == 0) {
                k = 0;
                continue;
            }
            if (k > 0) k--;
            while (genome_.sequence_[i + k] == genome_.sequence_[suffixArray_[j - 1] + k]) k++;
            if (k < std::numeric_limits<uint8_t>::max()) longestCommonPrefix_[j] = k;
            else longestCommonPrefix_[j] = std::numeric_limits<uint8_t>::max();
        }

        //update index
        std::vector<sjHash> sjHashRecord;
        sjHashRecord.resize(sjSeqLength * 2, {-1, 0});

        for (size_t i = 0; i < sjSeqLength * 2; ++i) {
            int32_t hash = 0;
            bool foundInvalid = false;
            for (unsigned int j = 0;j < kMerSize_; ++j) {
                if (i+j >= sjSeqLength * 2){
                    foundInvalid = true;
                    if (j == 0) sjHashRecord[i] = {-1, 0};
                    else {
                        hash = ((hash + 1) << (2 * (kMerSize_ - j))) - 1;
                        sjHashRecord[i] = {hash,  j};
                    }
                    break;
                }
                char c = sjdb.sjdbSeq_[i + j];
                if (c == '#' || c == 'N') {
                    foundInvalid = true;
                    if (j == 0) sjHashRecord[i] = {-1, 0};
                    else {
                        hash = ((hash + 1) << (2 * (kMerSize_ - j))) - 1;
                        sjHashRecord[i] = {hash,  j};
                    }
                    break;
                }
                hash = (hash << 2) | charToIndex(c);
            }
            if(!foundInvalid) sjHashRecord[i] = {hash, kMerSize_};
        }

        std::sort(sjHashRecord.begin(),sjHashRecord.end());

        int64_t nowShift = 0;
        int32_t prevHashInsert = -1;
        int32_t nowHashInsert = 0;

        for (size_t i = 0; i < sjSeqLength * 2; ++i) {
            if (sjHashRecord[i].hash == -1) continue;
            nowHashInsert = sjHashRecord[i].hash;
            int32_t length = sjHashRecord[i].length;
            if (length == 0) continue; // should not happen
            if (nowHashInsert != prevHashInsert) {
                if (prevHashInsert != -1) {
                    for (int32_t j = prevHashInsert + 1; j <= nowHashInsert; ++j) {
                        int64_t nowLeftSAIndex = patternMerMap_.get(j,patternMerMap_.INDEX_LEFT_SA_INDEX);
                        patternMerMap_.set(j,patternMerMap_.INDEX_LEFT_SA_INDEX,nowLeftSAIndex + nowShift);
                    }
                    //rebuild secondary index for prevHashInsert
                    buildExtendedIndexHashSingle(prevHashInsert);

                }
                prevHashInsert = nowHashInsert;
            }
            ++nowShift;
            if (length == (int) kMerSize_) {
                int64_t nowUpperRange = patternMerMap_.get(nowHashInsert,patternMerMap_.INDEX_UPPER_RANGE);
                patternMerMap_.set(nowHashInsert,patternMerMap_.INDEX_UPPER_RANGE,nowUpperRange + 1);
            }
            /**
             * TODO: The casting here makes me worried. Are you sure you want
             *       `patternMerMap_.get` to return signed valued?
             *       If so, why cast it to 32 bits here?
             */
             /**
              * `length` is at most `kMerSize_`. Therefore it only needs 4 bits to store,
              *  so it should be safe to cast it to 32 bits.
              */
            // int32_t originalLength = patternMerMap_.get(nowHashInsert,patternMerMap_.INDEX_LENGTH);
            int64_t originalLength = patternMerMap_.get(nowHashInsert,patternMerMap_.INDEX_LENGTH);
            if (originalLength < (int64_t) length) {
                patternMerMap_.set(nowHashInsert,patternMerMap_.INDEX_LENGTH,length);
            }
        }

        if (prevHashInsert != -1) {
            /**
             * TODO: This comparison also looks unintended.
             *       Comparing a 32 bit int with a 64 bit unsigned is fishy.
             *       I will cast `kMerNum_` to signed value and abort if it overflowed.
             */
             /**
              * kMerNum_ is at most 4^kMerSize_, which is at most 2^30
              * I will change the type to 32 bits to avoid confusion
              */
            uint64_t max_int64 = static_cast<uint64_t>(std::numeric_limits<int64_t>::max());
            if (kMerNum_ > max_int64) {
                std::cerr << "FATAL: kMerNum_ overflow!";
                std::abort();
            }

            for (int32_t j = prevHashInsert + 1; (int64_t) j < (int64_t) kMerNum_; ++j) {
                /**
                 * TODO: Since `patternMerMap_.get` expects unsigned values, are we completely sure that
                 *       `j` cannot be negative?
                 *       Again, I will put a runtime check for now and abort if it is not.
                 */

                /**
                 * yes, j should never be negative
                 */
                if (j < 0) {
                    std::cerr << "FATAL: index overflow!";
                    std::abort();
                }
                int64_t nowLeftSAIndex = patternMerMap_.get(j, patternMerMap_.INDEX_LEFT_SA_INDEX);
                patternMerMap_.set(j,patternMerMap_.INDEX_LEFT_SA_INDEX,nowLeftSAIndex + nowShift);
            }
            // rebuild secondary index for prevHashInsert
            buildExtendedIndexHashSingle(prevHashInsert);
        }



    }

    void Genome::modifyGenome(SJDB &sjdb) {
        if (sjdb.sjdbLoci_.empty()) return;
        size_t sjNum = sjdb.sjdbLoci_.size();
        sjdbStart_ = sequence_.length();
        sjdb.sjdbSeq_.resize(2 * sjNum * sjdb.sjdbLength);
        auto *sjdbStart = new size_t[sjNum];
        auto *sjdbEnd = new size_t[sjNum];
        auto *sjdbMotif = new uint8_t[sjNum];
        auto *sjdbShiftLeft = new uint8_t[sjNum];
        auto *sjdbShiftRight = new uint8_t[sjNum];

        std::string prevChr;
        size_t iChr = 0;
        for (size_t i = 0; i < sjNum; ++i) {
            if (prevChr != sjdb.sjdbLoci_[i].chr) {
                for (iChr = 0; iChr < chromosomes_.size(); ++iChr) {
                    if (chromosomes_[iChr].name == sjdb.sjdbLoci_[i].chr) break;
                }
                if (iChr >= chromosomes_.size()) {
                    // todo report error
                    std::cerr << "Warning: chromosome " << sjdb.sjdbLoci_[i].chr
                              << " not found in genome, skipping junction\n";
                    continue;
                }
                prevChr = sjdb.sjdbLoci_[i].chr;
            }

            sjdbStart[i] = sjdb.sjdbLoci_[i].start + chromosomes_[iChr].start - 1;
            sjdbEnd[i] = sjdb.sjdbLoci_[i].end + chromosomes_[iChr].start - 1;

            //judge the junction motif
            const char left1 = sequence_[sjdbStart[i]];
            const char left2 = sequence_[sjdbStart[i] + 1];
            const char right1 = sequence_[sjdbEnd[i] - 1];
            const char right2 = sequence_[sjdbEnd[i]];
            if (left1 == 'G' && left2 == 'T' && right1 == 'A' && right2 == 'G') sjdbMotif[i] = GTAG;
            else if (left1 == 'C' && left2 == 'T' && right1 == 'A' && right2 == 'C') sjdbMotif[i] = CTAC;
            else if (left1 == 'G' && left2 == 'C' && right1 == 'A' && right2 == 'G') sjdbMotif[i] = GCAG;
            else if (left1 == 'C' && left2 == 'T' && right1 == 'G' && right2 == 'C') sjdbMotif[i] = CTGC;
            else if (left1 == 'A' && left2 == 'T' && right1 == 'A' && right2 == 'C') sjdbMotif[i] = ATAC;
            else if (left1 == 'G' && left2 == 'T' && right1 == 'A' && right2 == 'T') sjdbMotif[i] = GTAT;
            else sjdbMotif[i] = NON_CANONICAL;

            //calc repeat length
            size_t jjL = 0;
            size_t jjR = 0;

            while (jjL < sjdbStart[i] - 1
                   && sequence_[sjdbStart[i] - 1 - jjL] == sequence_[sjdbEnd[i] - jjL]
                   && sequence_[sjdbStart[i] - 1 - jjL] != '#'
                   && sequence_[sjdbStart[i] - 1 - jjL] != 'N'
                   && jjL < 255) {
                jjL++;
            }

            while (sjdbEnd[i] + jjR < sequence_.size()
                   && sequence_[sjdbStart[i] + jjR] == sequence_[sjdbEnd[i] + 1 + jjR]
                   && sequence_[sjdbStart[i] + jjR] != '#'
                   && sequence_[sjdbStart[i] + jjR] != 'N'
                   && jjR < 255) {
                jjR++;
            }

            sjdbShiftLeft[i] = (uint8_t) jjL;
            sjdbShiftRight[i] = (uint8_t) jjR;

            //todo report too long repeats

            sjdbStart[i] -= jjL;
            sjdbEnd[i] -= jjL;
        }

        std::vector<sjdbSortRecord> sjdbSortRecords;
        sjdbSortRecords.resize(sjNum);
        for (size_t i = 0; i < sjNum; ++i) {
            size_t shift = 0;
            switch (sjdb.sjdbLoci_[i].strand) {
                case '+':
                    shift = 0;
                    break;
                case '-':
                    shift = sequence_.length();
                    break;
                default:
                    shift = sequence_.length() * 2;
            }
            sjdbSortRecords[i] = {sjdbStart[i] + shift, sjdbEnd[i] + shift, i};
        }
        std::sort(sjdbSortRecords.begin(), sjdbSortRecords.end());

        std::vector<size_t> sjdbNewIndex;
        sjdbNewIndex.reserve(sjNum);

        for (size_t i = 0; i < sjNum; ++i) {
            size_t prevIndex;
            size_t nowIndex = sjdbSortRecords[i].originalIndex;
            size_t nowSjNum = sjdbNewIndex.size();
            if (nowSjNum > 0) {
                prevIndex = sjdbNewIndex.back();
            }

            if (nowSjNum == 0 || sjdbStart[nowIndex] != sjdbStart[prevIndex]
                || sjdbEnd[nowIndex] != sjdbEnd[prevIndex]) {
                sjdbNewIndex.push_back(nowIndex);
            } else if (sjdb.sjdbLoci_[nowIndex].priority < sjdb.sjdbLoci_[prevIndex].priority) {
                continue;
            } else if (sjdb.sjdbLoci_[nowIndex].priority > sjdb.sjdbLoci_[prevIndex].priority) {
                sjdbNewIndex.back() = nowIndex;
            } else if ((sjdbMotif[nowIndex] != NON_CANONICAL && sjdbMotif[prevIndex] == NON_CANONICAL) ||
                       ((sjdbMotif[nowIndex] != NON_CANONICAL) == (sjdbMotif[prevIndex] != NON_CANONICAL) &&
                        sjdbShiftLeft[nowIndex] < sjdbShiftLeft[prevIndex])) {
                sjdbNewIndex.back() = nowIndex;
            }
        }

        sjdbSortRecords.reserve(sjdbNewIndex.size());
        for (size_t i = 0; i < sjdbNewIndex.size(); ++i) {
            size_t shift = sjdbMotif[sjdbNewIndex[i]] == NON_CANONICAL ? 0 : sjdbShiftLeft[sjdbNewIndex[i]];
            sjdbSortRecords[i] = {
                    sjdbStart[sjdbNewIndex[i]] + shift,
                    sjdbEnd[sjdbNewIndex[i]] + shift,
                    sjdbNewIndex[i]
            };
        }

        std::sort(sjdbSortRecords.begin(), sjdbSortRecords.end());

        sjDataBase_.reserve(sjdbNewIndex.size());

        for (size_t i = 0; i < sjdbNewIndex.size(); ++i) {
            size_t nowIndex = sjdbSortRecords[i].originalIndex;
            bool replace = false;

            size_t nowSjNum = sjDataBase_.size();
            if (nowSjNum > 0 && sjDataBase_.back().start == sjdbStart[nowIndex] &&
                    sjDataBase_.back().end == sjdbEnd[nowIndex]) {
                size_t prevIndex = sjdbSortRecords[i - 1].originalIndex;
                if (sjdb.sjdbLoci_[nowIndex].priority < sjdb.sjdbLoci_[prevIndex].priority) {
                    continue;
                } else if (sjdb.sjdbLoci_[nowIndex].priority > sjdb.sjdbLoci_[prevIndex].priority) {
                    replace = true;
                } else if (sjDataBase_.back().strand > 0 && sjdb.sjdbLoci_[nowIndex].strand == '.') {
                    continue;
                } else if (sjDataBase_.back().strand == 0 && sjdb.sjdbLoci_[nowIndex].strand != '.') {
                    replace = true;
                } else if (sjDataBase_.back().strand > 0 && sjdbMotif[nowIndex] == NON_CANONICAL) {
                    sjDataBase_.back().strand = 0;
                    continue;
                } else if ((sjDataBase_.back().motif != NON_CANONICAL && sjdbMotif[nowIndex] == NON_CANONICAL) ||
                           (sjDataBase_.back().motif % 2 == (2 - sjDataBase_.back().strand))) {
                    continue;
                } else {
                    replace = true;
                }
            }
            uint8_t strand = 0;
            switch (sjdb.sjdbLoci_[i].strand) {
                case '+':
                    strand = 1;
                    break;
                case '-':
                    strand = 2;
                    break;
                default:
                    strand = 0;
            }
            if (replace) {
                if (strand == 0) {
                    if (sjDataBase_.back().motif != NON_CANONICAL) {
                        strand = 2 - (sjDataBase_.back().motif % 2);
                    }
                }

                sjDataBase_.back() = {
                        .start = sjdbSortRecords[i].start,
                        .end = sjdbSortRecords[i].end,
                        .motif = sjdbMotif[nowIndex],
                        .shiftLeft = sjdbShiftLeft[nowIndex],
                        .shiftRight = sjdbShiftRight[nowIndex],
                        .strand = strand
                };
            } else {
                sjDataBase_.push_back({
                                              .start = sjdbSortRecords[i].start,
                                              .end = sjdbSortRecords[i].end,
                                              .motif = sjdbMotif[nowIndex],
                                              .shiftLeft = sjdbShiftLeft[nowIndex],
                                              .shiftRight = sjdbShiftRight[nowIndex],
                                              .strand = strand
                                      });
            }
        }
        sjdbNum_ = sjDataBase_.size();
        sjDonorStart_.resize(sjdbNum_);
        sjAcceptorStart_.resize(sjdbNum_);

        size_t sjGStart = 0;
        for (size_t i = 0; i < sjdbNum_; ++i) {
            sjDonorStart_[i] = sjDataBase_[i].start - sjdb.sjdbOverhang;
            sjAcceptorStart_[i] = sjDataBase_[i].end + 1;

            if (sjDataBase_[i].motif == NON_CANONICAL) {
                sjDonorStart_[i] += sjDataBase_[i].shiftLeft;
                sjAcceptorStart_[i] += sjDataBase_[i].shiftLeft;
            }

            memcpy((void *) (sjdb.sjdbSeq_.c_str() + sjGStart), sequence_.c_str() + sjDonorStart_[i],
                   sjdb.sjdbOverhang);
            memcpy((void *) (sjdb.sjdbSeq_.c_str() + sjGStart + sjdb.sjdbOverhang),
                   sequence_.c_str() + sjAcceptorStart_[i], sjdb.sjdbOverhang);
            memset((void *) (sjdb.sjdbSeq_.c_str() + sjGStart + sjdb.sjdbOverhang * 2),'#',SJDB_PADDING_LENGTH);
            sjGStart += sjdb.sjdbLength;
        }
        sjdb.sjdbSeq_.resize(2*sjdbNum_*sjdb.sjdbLength);



        delete[] sjdbStart;
        delete[] sjdbEnd;
        delete[] sjdbMotif;
        delete[] sjdbShiftLeft;
        delete[] sjdbShiftRight;

    }
}