#include <string>

/**
 * Error codes for handling genome index files or data structures.
 */

 enum class GenomeLoadingError {
    // NOTE: Ideally, you do not even have to bother giving specific numbers to these.
    OK = 0,
    FILE_ERROR = -1,
    MISSING_HEADER = -2,
    INVALID_HEADER = -3,
    MISSING_VERSION = -4,
    MISSING_ADDITIONAL_INDEX_TYPE = -5,
    MISSING_LCP_SIZE = -6,
    MISSING_LCP = -7,
    MISSING_EXTENDED_HASH_SIZE = -8,
    MISSING_EXTENDED_HASH = -9,
    // etc ...
};

/**
 * TODO: If you feel like it, add descriptions and macros for reporting what each error
 *       means and what to do.
 */

constexpr std::string_view getErrorDescription(GenomeLoadingError code) {
    switch (code) {
        case GenomeLoadingError::OK:
            return "";
        case GenomeLoadingError::FILE_ERROR:
            return "Could not open requested file ...";
        // etc.
    }
    return "Unknown error code.";
}

#define REPORT_GENOME_LOADING_ERROR(err) std::cerr << "[ERROR] " << getErrorDescription(err) << std::endl;