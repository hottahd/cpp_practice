#pragma once
#include <cassert>
#include <vector>
#include "grid.hpp"
#include "time.hpp"

template <typename Real>
struct Quantity {
    std::vector<Real> q0, q1, q2;
    
    Quantity(const Grid<Real>& grid){
        q0.resize(grid.i_total);
        q1.resize(grid.i_total);
        q2.resize(grid.i_total);
    }

    void save(const Config& config, const Time<Real>& time) const;
    void initial_condition(const Grid<Real>& grid);
};