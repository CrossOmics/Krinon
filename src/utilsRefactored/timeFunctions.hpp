#ifndef RNAALIGNREFACTORED_TIMEFUNCTIONS_HPP
#define RNAALIGNREFACTORED_TIMEFUNCTIONS_HPP
#include <iostream>
#include <ctime>
#include <cstring>
#include <string>

std::string getCurrentTimeString() {

    std::time_t rawtime;
    std::time(&rawtime);

    std::tm* timeinfo = std::localtime(&rawtime);

    char buffer[80];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", timeinfo);

    return std::string(buffer);
}
#endif //RNAALIGNREFACTORED_TIMEFUNCTIONS_HPP
