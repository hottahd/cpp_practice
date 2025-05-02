#include <iostream>
#include <string>
#include <cstdlib>

#include "types.hpp"
#include "config.hpp"
#include "time.hpp"
#include "grid.hpp"
#include "mhd.hpp"
#include "model.hpp"

int main() {
    Config config("config/config.json");
    config.create_save_directory();
    config.save();

    Grid<Real> grid = Grid<Real>::from_config(config.config_json);
    grid.save(config);

    Real tend = 1.0, dt_output = 0.1;
    Time<Real> time = Time<Real>::from_config(config.config_json);

    MHD<Real> mhd(grid);
    mhd.initial_condition(grid);
    mhd.save(config, time);

    Model<Real> model(config, time, grid, mhd);

    // advection.io_step();

    // while (advection.time.time < advection.time.tend) {
    //     advection.update();
    //     advection.time.update();

    //     advection.io_step();
    // };
        
    return 0;
}