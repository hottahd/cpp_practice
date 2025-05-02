#pragma once

#include "time.hpp"
#include "grid.hpp"
#include "mhd.hpp"
#include <vector>

template <typename Real>
struct Model {
    Config config;
    Time<Real> time;
    Grid<Real> grid;
    MHD<Real> mhd;

    Model(Config config_, Time<Real> time_, Grid<Real> grid_, MHD<Real> mhd_)
        : config(config_), time(time_), grid(grid_), mhd(mhd_) {}
    // void cfl_condition();
    // void update();
    // void io_step();
    // void sc2ssprk();
    // static void sc2(std::vector<Real>&  qq, std::vector<Real>& dqq, Grid<Real>& grid, Real& dt, Real& vc);
    // static void bc(std::vector<Real>& qq, Grid<Real>& grid);
};
