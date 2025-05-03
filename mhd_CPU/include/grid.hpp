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
    std::vector<Real> x, y, z, dx, dy, dz, dxi, dyi, dzi;
        
    Grid(int i_size_, int j_size_, int k_size_, int margin_, Real xmin_, Real xmax_, Real ymin_, Real ymax_, Real zmin_, Real zmax_) :
        i_size(i_size_), j_size(j_size_), k_size(k_size_), margin(margin_),
        xmin(xmin_), xmax(xmax_), ymin(ymin_), ymax(ymax_), zmin(zmin_), zmax(zmax_) {

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

            i_total = i_size + 2*x_margin;
            j_total = j_size + 2*y_margin;
            k_total = k_size + 2*z_margin;
    
    
            auto set_coordinate = [](std::vector<Real>& x, std::vector<Real>& dx, std::vector<Real>& dxi, int i_size, int i_total, int x_margin, Real xmin, Real xmax) {
                x.resize(i_total);
                dx.resize(i_total);
                dxi.resize(i_total);
    
                double dx0 = (xmax - xmin) / i_size;
                for (int i = 0; i < i_total; ++i) {
                    // dx at i + 1/2
                    dx[i] = dx0;
                    dxi[i] = 1.0 / dx[i];
                }
    
                x[x_margin] = xmin + 0.5*dx[x_margin];
                for (int i = x_margin + 1; i < i_total; ++i) {
                    x[i] = x[i - 1] + dx[i - 1];
                }
                for (int i = x_margin - 1; i >= 0; --i) {
                    x[i] = x[i + 1] - dx[i];
                }
            };
            set_coordinate(x, dx, dxi, i_size, i_total, x_margin, xmin, xmax);
            set_coordinate(y, dy, dyi, j_size, j_total, y_margin, ymin, ymax);
            set_coordinate(z, dz, dzi, k_size, k_total, z_margin, zmin, zmax);
        }

    void save(const Config& config) const {
        std::ofstream ofs_bin(config.save_dir + "/grid.bin", std::ios::binary);
        assert(ofs_bin.is_open());

        auto write_array = [&ofs_bin](const std::vector<Real>& x) {
            ofs_bin.write(reinterpret_cast<const char*>(x.data()), sizeof(Real)*x.size());
        };
        write_array(x);
        write_array(y);
        write_array(z);
    }  


    static Grid from_config(const json& json_obj) {
        return Grid(
            json_obj.at("grid").at("i_size").get<int>(), 
            json_obj.at("grid").at("j_size").get<int>(),
            json_obj.at("grid").at("k_size").get<int>(),
            json_obj.at("grid").at("margin").get<int>(),
            json_obj.at("grid").at("xmin").get<Real>(),
            json_obj.at("grid").at("xmax").get<Real>(),
            json_obj.at("grid").at("ymin").get<Real>(),
            json_obj.at("grid").at("ymax").get<Real>(),
            json_obj.at("grid").at("zmin").get<Real>(),
            json_obj.at("grid").at("zmax").get<Real>()
        );
    }
};