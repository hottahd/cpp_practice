#pragma once
#include <cassert>
#include <vector>
#include "grid.hpp"
#include "time.hpp"

template <typename Real>
struct Quantity {
    std::vector<Real> q0, q1, q2;
    Real* q0_dev;
    Real* q1_dev;
    Real* q2_dev;
    Real* qq_dev;
    Real* dqq_dev;
    
    Quantity(const Grid<Real>& grid){
        q0.resize(grid.i_total);
        q1.resize(grid.i_total);
        q2.resize(grid.i_total);
        q0_dev = nullptr;
        q1_dev = nullptr;
        q2_dev = nullptr;
        qq_dev = nullptr;
        dqq_dev = nullptr;
    }
    
    void save(const Config& config, const Time<Real>& time) const;
    void initial_condition(const Grid<Real>& grid);
    void allocate_device(const Grid<Real>& grid);
    void upload_q0();
    void download_q0();
    void free_device();
};