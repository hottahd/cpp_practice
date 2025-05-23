#pragma once

#include "time.hpp"
#include "grid.hpp"
#include "quantity.hpp"
#include <vector>

template <typename Real>
struct Advection {
    Real vc; // advection velocity
    Config config;
    Time<Real> time;
    Grid<Real> grid;
    Quantity<Real> quantity;

    Advection(Real vc_, Config config_, Time<Real> time_, Grid<Real> grid_, Quantity<Real> quantity_)
        : vc(vc_), config(config_), time(time_), grid(grid_), quantity(quantity_) {
    }
    void cfl_condition();
    void update();
    void io_step();
    void sc2ssprk();
    static void sc2(std::vector<Real>&  qq, std::vector<Real>& dqq, Grid<Real>& grid, Real& dt, Real& vc);
    static void bc(std::vector<Real>& qq, Grid<Real>& grid);
};
