#pragma once
#include "config.hpp"

template <typename Real>
struct EOS {
    Real gm;

    EOS(const Config& config) : 
        gm(config.json_obj.at("eos").at("gm").get<Real>()) {}
};