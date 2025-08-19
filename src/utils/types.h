#ifndef RNAALIGNMENT_TYPES_H
#define RNAALIGNMENT_TYPES_H
#include <cstdint>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <iostream>

namespace rna{
    using GenomePos = int64_t;
    using ReadPos = int64_t;
    using Score = int64_t;
    using Length = int64_t;
    struct Read {
        std::string sequence[2]; // 0 for forward, 1 for complementary
        std::string quality;
        std::string name;
        int64_t length;
    };

    struct MatchResultStored {
        uint8_t length;
        int64_t leftSAIndex;
        MatchResultStored() : length(0), leftSAIndex(-1) {};
    };

    struct MatchStats {
        int64_t maxLength;
        size_t count;
        size_t LongestMatchLeftPosInSA;
        size_t LongestMatchRightPosInSA;

        MatchStats(size_t len = 0, size_t cnt = 0,size_t l=0,size_t r = 0) : maxLength(len), count(cnt),LongestMatchLeftPosInSA(l),LongestMatchRightPosInSA(r) {}
    };

    struct Split{
        ReadPos readStart{0};
        size_t length{0};
        int direction{0};
    };

    struct Align {
        ReadPos readStart{0};
        int64_t length{0};
        size_t leftSAIndex{0};
        size_t rightSAIndex{0};
        size_t rep{0};
        int direction{0};
        bool isAnchor{false};

        bool operator<(const Align& other) const {
            return readStart < other.readStart;
        }
    };

    struct WindowAlign{
        ReadPos readStart{0};
        GenomePos genomeStart{0};
        int64_t length{0};
        int64_t score{0};
        bool isAnchor{false};
        bool operator<(const WindowAlign& other) const {
            if(genomeStart != other.genomeStart) return genomeStart < other.genomeStart;
            return length < other.length;
        }
    };

    struct Window {
        int chrIndex{0};
        int direction{0};
        int numAnchors{0};
        int64_t startBin{0};
        int64_t endBin{0};
        std::vector<WindowAlign> aligns;
        int minLengthWhenFull{0};
    };


    struct SpliceJunction {
        int64_t type{0};
        int64_t length{0};
        bool isAnnotated{false};
    };

    struct SJDBOutput{
        std::string chr;
        int64_t left;
        int64_t right;
        bool annotated;
        int64_t type;
        int strand{0}; // 0 undefined, 1 forward, 2 reverse
    };


    struct Exon {
        GenomePos start{0};
        int64_t length{0};
        ReadPos readStart{0};
        Score score{0};
    };


    struct Transcript {
        std::string readName;
        std::string chr;
        int strand{0};
        int64_t matched{0};
        int64_t unmatched{0};
        int64_t nIns{0};
        int64_t nDel{0};
        std::vector<Exon> exons;
        std::vector<SpliceJunction> sj;
        std::vector<WindowAlign> aligns;//Actually no need to store this, for debugging only
        int64_t readStart{0};
        int64_t genomeStart{0};
        int64_t posInChr{0};
        int64_t score;
        std::string getCIGAR() const {
            std::string cigar;
            if (readStart > 0) {
                cigar += std::to_string(readStart) + "S"; // soft clipping at the start
            }
            for (int64_t i = 0; i < exons.size(); ++i) {
                cigar += std::to_string(exons[i].length) + "M"; // exon
                if (i< exons.size() - 1) {
                    if(sj[i].type >= 0) {
                        cigar += std::to_string(sj[i].length) + "N"; // splice junction
                    } else if (sj[i].type == -1){
                        cigar += std::to_string(sj[i].length) + "D"; // deletion
                    } else if (sj[i].type == -2){
                        cigar += std::to_string(-sj[i].length) + "I"; // insertion
                    }
                }
            }

            return cigar;
        }
    };

    struct SAMEntry {
        std::string readName;
        int flag{0};
        std::string chr{"*"};
        int64_t pos{0};// 1-based
        int64_t mapq{0};
        std::string CIGAR;
        std::string RNEXT{"*"};
        int64_t PNEXT{0};
        int64_t TLEN{0};
        std::string seq;
        std::string qual;
        SAMEntry(Transcript &t,Read &read){
            readName = read.name;
            flag = t.strand == 0 ? 0 : 16; // 0 for forward, 16 for reverse
            chr = t.chr;
            pos = t.posInChr + 1; // SAM is 1-based
            mapq = 255; // default value
            CIGAR = t.getCIGAR();
            RNEXT = "*";
            PNEXT = 0;
            TLEN = 0; // not used
            seq = read.sequence[0];
            qual = "*"; // quality not available in this context
        }
        friend std::ostream& operator<<(std::ostream& os, const SAMEntry& entry) {
            os << entry.readName << "\t"
               << entry.flag << "\t"
               << entry.chr << "\t"
               << entry.pos << "\t"
               << entry.mapq << "\t"
               << entry.CIGAR << "\t"
               << entry.RNEXT << "\t"
               << entry.PNEXT << "\t"
               << entry.TLEN << "\t"
               << entry.seq << "\t"
               << entry.qual;
            return os;
        }

    };


    using ReadPtr = std::shared_ptr<Read>;
    using WindowPtr = std::shared_ptr<Window>;
    using TranscriptPtr = std::shared_ptr<Transcript>;

}
#endif //RNAALIGNMENT_TYPES_H
