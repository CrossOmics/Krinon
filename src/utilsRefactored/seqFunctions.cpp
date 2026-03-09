#include <cstdint>
#include "defines.h"
#include "seqFunctions.h"

namespace RefactorProcessing {

    // encode a k-mer into an integer hash, return -1 if the seq length is less than kMerSize, or -i-2 if the i-th character is invalid
    int64_t encodeKmer(const std::string_view &seq, int kMerSize) {
        // TODO: `kMerSize` should always be unsigned!
        if (seq.length() < (size_t) kMerSize) {
            return -1;
        }
        int64_t hash = 0;
        for (int i = 0; i < kMerSize; ++i) {
            int32_t idx = charToIndex(seq[i]);
            if (idx < 0) return -i - 2;
            hash = (hash << 2) | idx;
        }
        return hash;
    }



    bool writeString(std::ofstream &ofs, const std::string &s) {
        uint32_t n = static_cast<uint32_t>(s.size());
        ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));
        if (n) ofs.write(s.data(), n);
        return bool(ofs);
    }

     bool readString(std::ifstream &ifs, std::string &s) {
        uint32_t n;
        ifs.read(reinterpret_cast<char*>(&n), sizeof(n));
        if (!ifs) return false;
        s.resize(n);
        if (n) ifs.read(&s[0], n);
        return bool(ifs);
    }


}
