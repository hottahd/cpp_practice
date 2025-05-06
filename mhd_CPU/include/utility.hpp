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

    template <typename T>
    inline T pow2(T x) {
        return x * x;
    }

    template <typename T>
    inline T pow3(T x) {
        return x * x * x;
    }
    template <typename T>
    inline T pow4(T x) {
        return x * x * x * x;
    }
};

