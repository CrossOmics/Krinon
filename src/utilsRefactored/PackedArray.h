#ifndef RNAALIGNREFACTORED_PACKEDARRAYR_H
#define RNAALIGNREFACTORED_PACKEDARRAYR_H
//todo

#include <string>
namespace RefactorProcessing {
    class PackedArray{
    private:
        // configs
        int wordLengthBits_; // number of bits for each element
        int wordCompLength_;
        uint64_t bitMask_;
        int64_t wordNum_; // number of elements stored
        int64_t arrayLength_; // length of uint64_t data array
        int64_t reservedLength_; // length of reserved num in the last uint64_t element


        //data
        uint64_t* data;
        bool allocated;
        // note: notice the management of this pointer!

    public:
        PackedArray();
        ~PackedArray();

        void initialize(int64_t wordNum,int wordLengthBits,int reservedLength=0);
        void setValue(int64_t index, uint64_t value);
        uint64_t getValue(int64_t index) const;
        uint64_t operator[](int64_t index) const;
        // build from an uint64_t array, the array should be pre-allocated. Should reserve space for reservedLength.
        void buildFromReservedArray(uint64_t* array,uint32_t wordLengthBits, uint64_t length,uint64_t reservedLength);
        bool setReservedLength(int64_t newReservedLength); // return false if newReservedLength is larger than current reservedLength

        void setLength(int64_t length);
        // input/output
        int writeToFile(const std::string& fileName) const;
        int loadFromFile(const std::string& fileName, int reservedLength=0);

    };
}
#endif //RNAALIGNREFACTORED_PACKEDARRAY_H
