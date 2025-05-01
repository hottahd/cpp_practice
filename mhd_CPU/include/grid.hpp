#pragma once
#include <vector>
#include <cassert>

#include <nlohmann/json.hpp>

#include "config.hpp"

using json = nlohmann::json;

template <typename Real>
struct Grid {
    int i_size, j_size, k_size;
    int i_total, j_total, k_total;
    int is, js, ks;
    int margin;
    int x_margin, y_margin, z_margin;
    Real xmin, xmax, ymin, ymax, zmin, zmax;
    std::vector<Real> x, y, z, dx, dy, dz;
        
    Grid(int i_size_, int j_size_, int k_size_, int margin_,
        Real xmin_, Real xmax_, Real ymin_, Real ymax_, Real zmin_, Real zmax_)
        : i_size(i_size_), j_size(j_size_), k_size(k_size_),
        margin(margin_), 
        xmin(xmin_), xmax(xmax_),
        ymin(ymin_), ymax(ymax_),
        zmin(zmin_), zmax(zmax_)
        {
            assert(i_size > 0);
            assert(j_size > 0);
            assert(k_size > 0);
            assert(margin >= 0);
            assert(xmax > xmin);
            assert(ymax > ymin);
            assert(zmax > zmin);

            is = i_size > 1 ? 1 : 0;
            js = j_size > 1 ? 1 : 0;
            ks = k_size > 1 ? 1 : 0;

            x_margin = margin*is;
            y_margin = margin*js;
            z_margin = margin*ks;

            i_total = i_size + 2*margin;
            j_total = j_size + 2*margin;
            k_total = k_size + 2*margin;
    
    
            auto set_x = [](std::vector<Real>& x, std::vector <Real>& dx, int i_size, int i_total, int x_margin, Real xmin, Real xmax) {
                x.resize(i_total);
                dx.resize(i_total);
    
                double dx0 = (xmax - xmin) / i_size;
                for (int i = 0; i < i_total; ++i) {
                    // dx at i + 1/2
                    dx[i] = dx0;
                }
    
                x[x_margin] = xmin + 0.5*dx[x_margin];
                for (int i = x_margin + 1; i < i_total; ++i) {
                    x[i] = x[i - 1] + dx[i - 1];
                }
                for (int i = x_margin - 1; i >= 0; --i) {
                    x[i] = x[i + 1] - dx[i];
                }
            };
            set_x(x, dx, i_size, i_total, x_margin, xmin, xmax);
            set_x(y, dy, j_size, j_total, y_margin, ymin, ymax);
            set_x(z, dz, k_size, k_total, z_margin, zmin, zmax);
        }

    void save(const Config& config) const {
        json j;
        j["i_size"] = i_size;
        j["j_size"] = j_size;
        j["k_size"] = k_size;    
        j["margin"] = margin;
        j["xmin"] = xmin;
        j["xmax"] = xmax;
        j["ymin"] = ymin;
        j["ymax"] = ymax;
        j["zmin"] = zmin;
        j["zmax"] = zmax;
    
        std::ofstream ofs_json(config.save_dir + "/grid.json");
        assert(ofs_json.is_open());
        ofs_json << j.dump(4);
    
        std::ofstream ofs_bin(config.save_dir + "/grid.bin", std::ios::binary);
        assert(ofs_bin.is_open());

        auto write_array = [&ofs_bin](const std::vector<Real>& x) {
            ofs_bin.write(reinterpret_cast<const char*>(x.data()), sizeof(Real)*x.size());
        };
        write_array(x);
        write_array(y);
        write_array(z);
    }  
};