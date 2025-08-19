#ifndef RNAALIGNMENT_EXCEPTIONS_H
#define RNAALIGNMENT_EXCEPTIONS_H
#include <string>
#include <stdexcept>
namespace rna {
    enum class ExceptionType {
        NONE = 0,
        FILE_NOT_FOUND,
        INVALID_INPUT,
        FORMAT_ERROR,
        OUT_OF_RANGE,
        NOT_IMPLEMENTED,
    };

    class RNAException : public std::exception {
    public:
        RNAException(ExceptionType code, const std::string& message)
                : code_(code), message_(message) {}

        const char* what() const noexcept override {
            return message_.c_str();
        }

        ExceptionType code() const { return code_; }

    private:
        ExceptionType code_;
        std::string message_;
    };
}
#endif //RNAALIGNMENT_EXCEPTIONS_H
