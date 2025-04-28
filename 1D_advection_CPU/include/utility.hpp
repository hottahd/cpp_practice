#pragma once
#include <string>
#include <sstream>
#include <iomanip>

namespace util {
    inline std::string zfill(int num, int width) {
        std::ostringstream oss;
        oss << std::setw(width) << std::setfill('0') << num;
        return oss.str();
    }
};