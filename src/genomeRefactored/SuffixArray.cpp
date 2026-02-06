#include <fstream>
#include "../utilsRefactored/defines.h"
#include "SuffixArray.h"
#include "../utilsRefactored/seqFunctions.h"

namespace RefactorProcessing{
    void SuffixArray::build(const std::string &seq) {
        buildSAIS(seq);
    }

    void SuffixArray::setParam(const Parameters &P){
        reservedLength = P.limitSjdbInsertN;
    }

    PackedArray& SuffixArray::getSuffixArray(){
        return suffixArray_;
    }

    int SuffixArray::writeToFile(const std::string &fileName) const {
        std::ofstream ofs(fileName, std::ios::binary);
        if (!ofs.is_open()) return -1;

        // magic + version
        const char magic[4] = {'S','U','F','X'}; // SUFX marker
        ofs.write(magic, 4);
        uint32_t version = 1;
        ofs.write(reinterpret_cast<const char*>(&version), sizeof(version));

        // write buildMethod
        if (!writeString(ofs, buildMethod)) return -2;

        // write numeric members
        ofs.write(reinterpret_cast<const char*>(&reservedLength), sizeof(reservedLength));
        ofs.write(reinterpret_cast<const char*>(&length_), sizeof(length_));
        ofs.write(reinterpret_cast<const char*>(&wordBits_), sizeof(wordBits_));

        if (!ofs) return -3;
        ofs.close();

        // write packed array to a separate file for speed/modularity
        std::string packedFile = fileName + ".packed";

        int ret = suffixArray_.writeToFile(packedFile);
        if (ret != 0) return -4;

        return 0;
    }

    int SuffixArray::loadFromFile(const std::string &fileName) {
        std::ifstream ifs(fileName, std::ios::binary);
        if (!ifs.is_open()) return -1;

        char magic[4];
        ifs.read(magic, 4);
        if (!ifs) return -2;
        if (!(magic[0]=='S' && magic[1]=='U' && magic[2]=='F' && magic[3]=='X')) return -3;

        uint32_t version = 0;
        ifs.read(reinterpret_cast<char*>(&version), sizeof(version));
        if (!ifs) return -4;

        if (!readString(ifs, buildMethod)) return -5;

        ifs.read(reinterpret_cast<char*>(&reservedLength), sizeof(reservedLength));
        ifs.read(reinterpret_cast<char*>(&length_), sizeof(length_));
        ifs.read(reinterpret_cast<char*>(&wordBits_), sizeof(wordBits_));

        if (!ifs) return -6;
        ifs.close();

        std::string packedFile = fileName + ".packed";
        int ret = suffixArray_.loadFromFile(packedFile);
        if (ret != 0) return -7;

        return 0;
    }


}
