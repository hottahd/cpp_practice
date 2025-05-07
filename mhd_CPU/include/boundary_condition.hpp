#pragma once

#include <string>

#include "model.hpp"

template <typename Real>
struct BoundaryCondition {
    Config& config;
    Grid<Real>& grid;
    EOS<Real>& eos;
    MHD<Real>& mhd;

    BoundaryCondition(Model<Real>& model)
        : config(model.config),
          grid(model.grid),
          eos(model.eos),
          mhd(model.mhd) {}

    void bnd_symmetric(Array3D<Real>& arr, Array3D<Real>* fac, Real sign, std::string direction, std::string bnd_type){
        int i0_, i1_, j0_, j1_, k0_, k1_;
        bnd_range(i0_, i1_, j0_, j1_, k0_, k1_, direction);
        for (int i = i0_; i < i1_; ++i) {
            for (int j = j0_; j < j1_; ++j) {
                for (int k = k0_; k < k1_; ++k) {
                    int i_ghst = i, i_trgt = i;
                    int j_ghst = j, j_trgt = j;
                    int k_ghst = k, k_trgt = k;
                    if (direction == "x") {
                        bnd_symmetric_index(i, grid.i_total, grid.i_margin, i_ghst, i_trgt, bnd_type);
                    } else if (direction == "y") {
                        bnd_symmetric_index(j, grid.j_total, grid.j_margin, j_ghst, j_trgt, bnd_type);
                    } else if (direction == "z") {
                        bnd_symmetric_index(k, grid.k_total, grid.k_margin, k_ghst, k_trgt, bnd_type);
                    }
                    if (fac) {
                        arr(i_ghst, j_ghst, k_ghst) = sign * (*fac)(i_trgt, j_trgt, k_trgt) * arr(i_trgt, j_trgt, k_trgt);
                    } else {
                        arr(i_ghst, j_ghst, k_ghst) = sign * arr(i_trgt, j_trgt, k_trgt);
                    }
                }
            }
        }
    };

    void bnd_periodic(Array3D<Real>& arr, std::string direction, std::string bnd_type) {
        int i0_, i1_, j0_, j1_, k0_, k1_;
        bnd_range(i0_, i1_, j0_, j1_, k0_, k1_, direction);
        for (int i = i0_; i < i1_; ++i) {
            for (int j = j0_; j < j1_; ++j) {
                for (int k = k0_; k < k1_; ++k) {
                    int i_ghst = i, i_trgt = i;
                    int j_ghst = j, j_trgt = j;
                    int k_ghst = k, k_trgt = k;
                    if (direction == "x") {
                        bnd_periodic_index(i, grid.i_total, grid.i_margin, i_ghst, i_trgt, bnd_type);
                    } else if (direction == "y") {
                        bnd_periodic_index(j, grid.j_total, grid.j_margin, j_ghst, j_trgt, bnd_type);
                    } else if (direction == "z") {
                        bnd_periodic_index(k, grid.k_total, grid.k_margin, k_ghst, k_trgt, bnd_type);
                    }
                    arr(i_ghst, j_ghst, k_ghst) = arr(i_trgt, j_trgt, k_trgt);
                }
            }
        }
    }

    inline void bnd_range(int& i0_, int& i1_, int& j0_, int& j1_, int& k0_, int& k1_, std::string direction) {
        i0_ = 0; i1_ = grid.i_total;
        j0_ = 0; j1_ = grid.j_total;
        k0_ = 0; k1_ = grid.k_total;
        if (direction == "x") { i1_ = grid.i_margin; }
        else if (direction == "y") { j1_ = grid.j_margin; }
        else if (direction == "z") { k1_ = grid.k_margin; }
    }

    inline void bnd_symmetric_index(int i, int i_total, int i_margin, int& i_ghst, int& i_trgt, std::string bnd_type) {
        if (bnd_type == "inner") {
            i_ghst = i;
            i_trgt = 2*i_margin - i - 1;
        } else if (bnd_type == "outer") {
            i_ghst = i_total - i_margin + i;
            i_trgt = i_total - i_margin - i - 1;
        }
    }

    inline void bnd_periodic_index(int i, int i_total, int i_margin, int& i_ghst, int& i_trgt, std::string bnd_type) {
        if (bnd_type == "inner") {
            i_ghst = i;
            i_trgt = i_total - 2*i_margin + i;
        } else if (bnd_type == "outer") {
            i_ghst = i_total - i_margin + i;
            i_trgt = i_margin + i;
        }
    }

    void apply(MHDCore<Real>& qq) {
        // Apply boundary conditions to the MHD core variables
        std::vector<std::string> directions = {"x", "y", "z"};
        std::vector<std::string> bnd_types  = {"inner", "outer"};
        
        for (const auto& dir : directions) {
            for (const auto& btype : bnd_types) {
                bnd_symmetric(qq.ro, nullptr, 1.0, dir, btype);
                bnd_symmetric(qq.vx, nullptr, 1.0, dir, btype);
                bnd_symmetric(qq.vy, nullptr, 1.0, dir, btype);
                bnd_symmetric(qq.vz, nullptr, 1.0, dir, btype);
                bnd_symmetric(qq.bx, nullptr, 1.0, dir, btype);
                bnd_symmetric(qq.by, nullptr, 1.0, dir, btype);
                bnd_symmetric(qq.bz, nullptr, 1.0, dir, btype);
                bnd_symmetric(qq.ei, nullptr, 1.0, dir, btype);
                bnd_symmetric(qq.ph, nullptr, 1.0, dir, btype);

                // bnd_periodic(qq.ro, dir, btype);
                // bnd_periodic(qq.vx, dir, btype);
                // bnd_periodic(qq.vy, dir, btype);
                // bnd_periodic(qq.vz, dir, btype);
                // bnd_periodic(qq.bx, dir, btype);
                // bnd_periodic(qq.by, dir, btype);
                // bnd_periodic(qq.bz, dir, btype);
                // bnd_periodic(qq.ei, dir, btype);
                // bnd_periodic(qq.ph, dir, btype);
}
        }        
    }
};
