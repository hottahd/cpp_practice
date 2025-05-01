#pragma once

#include <vector>

#include "time.hpp"
#include "grid.cuh"
#include "quantity.cuh"

template <typename Real>
struct Advection {
    Real vc; // advection velocity
    Config config;
    Time<Real> time;
    Grid<Real> grid;
    Grid_Device<Real> grid_device;
    Quantity<Real> quantity;

    Advection(Real vc_, Config config_, Time<Real> time_, Grid<Real> grid_, Quantity<Real> quantity_)
        : vc(vc_), config(config_), time(time_), grid(grid_), quantity(quantity_) {
            grid_device = grid.to_device();
    }
    void cfl_condition();
    void run();
    void update();
    void io_step();
    void sc2ssprk();
    static void bc(std::vector<Real>& qq, Grid<Real>& grid);
};
