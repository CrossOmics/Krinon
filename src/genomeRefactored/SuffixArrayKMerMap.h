#ifndef RNAALIGNREFACTORED_SUFFIXARRAYKMERMAP_H
#define RNAALIGNREFACTORED_SUFFIXARRAYKMERMAP_H
#include <cstdint>
#include <fstream>
namespace RefactorProcessing {
    class SuffixArrayKMerMap {
        // corresponding to SAi in STAR
        // compressed
        // a class to store k-mer data for suffix array
    private:
        uint64_t* data_;
        int ELEMENT_BIT_LENGTHS[3];
        int BIT_OFFSETS[3];
        bool dataLengthLessThanLL_;
        int64_t length_;
    public:
        SuffixArrayKMerMap(int INDEX_LENGTH_BITS = 4, int INDEX_LEFT_SA_INDEX_BITS = 36, int INDEX_UPPER_RANGE_BITS = 24);
        ~SuffixArrayKMerMap();
        static constexpr int INDEX_LENGTH = 0;
        static constexpr int INDEX_LEFT_SA_INDEX = 1;
        static constexpr int INDEX_UPPER_RANGE = 2;
        uint64_t EMPTY_UPPER_RANGE;
        uint64_t EMPTY_SA_INDEX;

        void init(uint64_t length);

        int64_t size() const;

        void set(uint64_t index, int elementIndex, int64_t value);

        int64_t get(uint64_t index, int elementIndex) const;

        void writeToFile(std::ostream &out,std::ofstream &logOut) const;

        void loadFromFile(std::ifstream &in, std::ifstream &logIn);



    };
}
#endif //RNAALIGNREFACTORED_SUFFIXARRAYKMERMAP_H
