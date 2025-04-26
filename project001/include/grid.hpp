#pragma once
#include <vector>
#include <cassert>
#include "config.hpp"

template <typename Real>
struct Grid {
    int i_size, i_total, margin;
    Real xmin, xmax;
    Real dx;
    std::vector<Real> x;
        
    Grid(int i_size_, int margin_, Real xmin_, Real xmax_);
    void save(const Config& config) const;
};