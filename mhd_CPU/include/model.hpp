#pragma once

#include <string>
#include <vector>

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

    static Model from_config_file(const std::string&load_filepath) {
        Config config = Config(load_filepath);
        Time<Real> time = Time<Real>::from_config(config.json_obj);
        Grid<Real> grid = Grid<Real>::from_config(config.json_obj);
        EOS<Real> eos = EOS<Real>(config);
        MHD<Real> mhd = MHD<Real>(grid);
    
        return Model(config, time, grid, eos, mhd);
    
    }
};

