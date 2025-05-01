#pragma once

#include <vector>
#include <fstream>
#include <cassert>

#include <nlohmann/json.hpp>

#include "config.hpp"

using json = nlohmann::json;

template <typename Real>
struct GridBase {
    int i_size, i_total, margin;
    Real xmin, xmax;
    Real dx;
};

template <typename Real>
struct Grid_Device : public GridBase<Real> {
    Real* x;
};

template <typename Real>
struct Grid : public GridBase<Real> {
    std::vector<Real> x;
        
    Grid(int i_size_, int margin_, Real xmin_, Real xmax_){
        this->i_size = i_size_;
        this->margin = margin_;
        this->xmin = xmin_;
        this->xmax = xmax_;
        this->i_total = this->i_size + 2*this->margin;
        this->dx = (this->xmax - this->xmin) / this->i_size;
  
        assert(this->i_size > 0);
        assert(this->margin >= 0);
        assert(this->xmax > this->xmin);
  
        x.resize(this->i_total);
        x[this->margin] = this->xmin + 0.5*this->dx;
        for (int i = this->margin + 1; i < this->i_total; ++i) {
          x[i] = x[i - 1] + this->dx;
        }
        for (int i = this->margin - 1; i >= 0; --i) {
          x[i] = x[i + 1] - this->dx;
        }
    };    

    void save(const Config& config) const {
        json j;
        j["i_size"] = this->i_size;
        j["margin"] = this->margin;
        j["xmin"] = this->xmin;
        j["xmax"] = this->xmax;
        j["dx"] = this->dx;
    
        std::ofstream ofs_json(config.save_dir + "/grid.json");
        assert(ofs_json.is_open());
        ofs_json << j.dump(4);
    
        std::ofstream ofs_bin(config.save_dir + "/grid.bin", std::ios::binary);
        assert(ofs_bin.is_open());
        ofs_bin.write(reinterpret_cast<const char*>(x.data()), sizeof(Real)*x.size());
    }

    Grid_Device<Real> to_device() const {
        Grid_Device<Real> grid_device;

        // Copy the base class data
        static_cast<GridBase<Real>&>(grid_device) = *this;
    
        cudaMalloc(&grid_device.x, sizeof(Real)*this->i_total);
        cudaMemcpy(grid_device.x, x.data(), sizeof(Real)*this->i_total, cudaMemcpyHostToDevice);
    
        return grid_device;
    };
        
};


