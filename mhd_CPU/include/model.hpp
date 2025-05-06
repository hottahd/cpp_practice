#pragma once

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>  // std::setprecision, std::setw

#include "time.hpp"
#include "grid.hpp"
#include "eos.hpp"
#include "mhd.hpp"

template <typename Real>
struct Model {
    Config config;
    Time<Real> time;
    Grid<Real> grid;
    EOS<Real> eos;
    MHD<Real> mhd;

    Model(Config config_, Time<Real> time_, Grid<Real> grid_, EOS<Real> eos_, MHD<Real> mhd_)
        : config(config_), time(time_), grid(grid_), eos(eos_), mhd(mhd_) {}

    void save_metadata() {
        this->config.create_save_directory();
        this->config.save();
        this->grid.save(this->config);
    }

    void save_state() {
        this->time.save(this->config);
        this->mhd.save(this->config, this->time);
    }

    void save_if_needed() {
        // std::cout << "### aaa ### " <<  this->time.time << " " << this->time.dt_output << " " << this->time.time - this->time.dt << " " << this->time.dt_output << std::endl;
        // std::cout << "### bbb ###" << std::floor( this->time.time                 /this->time.dt_output) << " " <<std::floor((this->time.time - this->time.dt)/this->time.dt_output) << std::endl;
        if(std::floor( this->time.time                 /this->time.dt_output) != 
           std::floor((this->time.time - this->time.dt)/this->time.dt_output)) {
          this->save_state();

          std::cout << std::fixed << std::setprecision(2)
          << "time = " << std::setw(6) << this->time.time
          << ";  n_step = " << std::setw(8) << this->time.n_step
          << ";  n_output = " << std::setw(8) << this->time.n_output
          << std::endl;

          time.n_output++;
        }
    }
    static Model from_config_file(const std::string&load_filepath) {
        Config config = Config(load_filepath);
        Time<Real> time = Time<Real>::from_config(config.json_obj);
        Grid<Real> grid = Grid<Real>::from_config(config.json_obj);
        EOS<Real> eos = EOS<Real>(config);
        MHD<Real> mhd = MHD<Real>(grid);
    
        return Model(config, time, grid, eos, mhd);
    
    }
};

