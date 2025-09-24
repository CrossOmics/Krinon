#ifndef RNAALIGNREFACTORED_UTILS_H
#define RNAALIGNREFACTORED_UTILS_H
std::string getTime(){

    std::time_t now = std::time(nullptr);

    std::tm* local_time = std::localtime(&now);

    char buffer[80];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", local_time);

    return std::string(buffer);
}
#endif //RNAALIGNREFACTORED_UTILS_H
