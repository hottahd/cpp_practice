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
        int i0, i1, j0, j1, k0, k1;
        bnd_range(i0, i1, j0, j1, k0, k1, direction);
        for (int i = i0; i < i1; ++i) {
            for (int j = j0; j < j1; ++j) {
                for (int k = k0; k < k1; ++k) {
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

    void bnd_range(int& i0, int& i1, int& j0, int& j1, int& k0, int& k1, std::string direction) {
        i0 = 0; i1 = grid.i_total;
        j0 = 0; j1 = grid.j_total;
        k0 = 0; k1 = grid.k_total;
        if (direction == "x") { i1 = grid.i_margin; }
        else if (direction == "y") { j1 = grid.j_margin; }
        else if (direction == "z") { k1 = grid.k_margin; }
    }

    void bnd_symmetric_index(int i, int i_total, int i_margin, int& i_ghst, int& i_trgt, std::string bnd_type) {
        if (bnd_type == "inner") {
            i_ghst = i;
            i_trgt = 2*i_margin - i - 1;
        } else if (bnd_type == "outer") {
            i_ghst = i_total - i_margin + i;
            i_trgt = i_total - i_margin - i - 1;
        }
    }
};
