#pragma once

#include <cassert>
#include <vector>

#include "array3d.hpp"
#include "grid.hpp"
#include "time.hpp"

template <typename Real>
struct MHD_Core {
    Array3D<Real> ro, vx, vy, vz, bx, by, bz, ei, ph;

    MHD_Core(int i_size, int j_size, int k_size)
        : ro(i_size, j_size, k_size),
          vx(i_size, j_size, k_size),
          vy(i_size, j_size, k_size),
          vz(i_size, j_size, k_size),
          bx(i_size, j_size, k_size),
          by(i_size, j_size, k_size),
          bz(i_size, j_size, k_size),
          ei(i_size, j_size, k_size),
          ph(i_size, j_size, k_size) {}
};

template <typename Real>
struct MHD {
    MHD_Core<Real> qq, qq_argm, qq_rslt;
    
    MHD(const Grid<Real>& grid)
        : qq     (grid.i_size, grid.j_size, grid.k_size),
          qq_argm(grid.i_size, grid.j_size, grid.k_size),
          qq_rslt(grid.i_size, grid.j_size, grid.k_size) {}

    void save(const Config& config, const Time<Real>& time) const;
    void initial_condition(const Grid<Real>& grid);
};