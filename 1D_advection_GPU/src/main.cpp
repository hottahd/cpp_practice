#include <iostream>
#include <string>
#include "config.hpp"
#include "grid.hpp"
#include "time.hpp"
#include "quantity.hpp"
#include "types.hpp"
#include "advection.cup"

int main() {
    std::string save_dir = "../data";
    Config config(save_dir);
    config.create_save_directory();
    config.save();

    int i_size = 128, margin = 2;
    Real xmin = 0.0, xmax = 1.0;
    Grid<Real> grid(i_size, margin, xmin, xmax);
    grid.save(config);

    Real tend = 1.0, dt_output = 0.1;
    Time<Real> time(tend, dt_output);
    
    Quantity<Real> quantity(grid);
    quantity.initial_condition(grid);
    quantity.save(config, time);

    Real vc = 1.0; // advection velocity
    Advection<Real> advection(vc, config, time, grid, quantity);

    advection.io_step();

    while (advection.time.time < advection.time.tend) {
        advection.update();
        advection.time.update();

        advection.io_step();
    };
        
    return 0;
}