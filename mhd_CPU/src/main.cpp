#include <iostream>
#include <string>
#include <cstdlib>

#include "mhd.hpp"
#include "config.hpp"
#include "grid.hpp"
#include "time.hpp"
#include "types.hpp"

int main() {
    std::string save_dir = "../data";
    Config config(save_dir);
    config.create_save_directory();
    config.save();

    int i_size = 128, j_size = 1, k_size = 1;
    int margin = 2;
    Real xmin = 0.0, xmax = 1.0;
    Real ymin = 0.0, ymax = 1.0;
    Real zmin = 0.0, zmax = 1.0;
    Grid<Real> grid(i_size, j_size, k_size, margin, xmin, xmax, ymin, ymax, zmin, zmax);
    grid.save(config);

    Real tend = 1.0, dt_output = 0.1;
    Time<Real> time(tend, dt_output);

    MHD<Real> mhd(grid);
    mhd.initial_condition(grid);
    mhd.save(config, time);

    // Real vc = 1.0; // advection velocity
    // Advection<Real> advection(vc, config, time, grid, quantity);

    // advection.io_step();

    // while (advection.time.time < advection.time.tend) {
    //     advection.update();
    //     advection.time.update();

    //     advection.io_step();
    // };
        
    return 0;
}