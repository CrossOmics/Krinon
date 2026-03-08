#include <fstream>
#include "../utilsRefactored/defines.h"
#include "SuffixArray.h"
#include "../utilsRefactored/seqFunctions.h"



namespace RefactorProcessing{
    void SuffixArray::build(const std::string &seq) {
        if (buildMethod == "SAIS") buildSAIS(seq);
        else if (buildMethod == "CAPS") buildCAPS(seq);
        else if (buildMethod == "Traditional") buildTraditional(seq);
        /*for (int i = 0; i<10;++i){
            std::cout << suffixArray_[i] << " ";
        }*/
        /*SuffixArray debug;
        debug.buildSAIS(seq);
        std::cout << "valid length:" << length_ <<"\n";
        for (int i = 0; i < seq.length();++i){
            if (suffixArray_[i] != debug[i]){
                std::cout << "bug at "<< i <<"\n";
                std::cout << suffixArray_[i] << " " << debug[i] << "\n";
                break;
            }

        }*/
        std::cout << "\n";
    }

    void SuffixArray::buildFromOtherInit(const SuffixArray &other, int newLength) {
        buildMethod = other.buildMethod;
        length_ = newLength;
        wordBits_ = other.wordBits_;
        reservedLength = other.reservedLength;
        suffixArray_.initialize(newLength, wordBits_, reservedLength);
    }

    void SuffixArray::setParam(const Parameters &P){
        reservedLength = 0; //todo
        threadNum = P.workingThreads;
        buildMethod = P.buildMethod;
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




        // write numeric members
        ofs.write(reinterpret_cast<const char*>(&reservedLength), sizeof(reservedLength));
        ofs.write(reinterpret_cast<const char*>(&length_), sizeof(length_));
        ofs.write(reinterpret_cast<const char*>(&wordBits_), sizeof(wordBits_));

        if (!ofs) return -3;
        ofs.close();

        // write packed array to a separate file for speed/modularity
        std::string packedFile = fileName;

        int ret = suffixArray_.writeToFile(packedFile+".array");
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



        ifs.read(reinterpret_cast<char*>(&reservedLength), sizeof(reservedLength));
        ifs.read(reinterpret_cast<char*>(&length_), sizeof(length_));
        ifs.read(reinterpret_cast<char*>(&wordBits_), sizeof(wordBits_));

        if (!ifs) return -6;
        ifs.close();

        std::string packedFile = fileName;
        int ret = suffixArray_.loadFromFile(packedFile+".array");
        if (ret != 0) return -7;

        return 0;
    }


}
