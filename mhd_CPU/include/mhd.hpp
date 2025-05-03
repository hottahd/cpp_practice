#pragma once

#include <cassert>
#include <vector>

#include "array3d.hpp"
#include "grid.hpp"
#include "time.hpp"
#include "utility.hpp"

template <typename Real>
struct MHDCore {
    Array3D<Real> ro, vx, vy, vz, bx, by, bz, ei, ph;

    MHDCore(int i_size, int j_size, int k_size)
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
    MHDCore<Real> qq, qq_argm, qq_rslt;
    
    MHD(const Grid<Real>& grid)
        : qq     (grid.i_total, grid.j_total, grid.k_total),
          qq_argm(grid.i_total, grid.j_total, grid.k_total),
          qq_rslt(grid.i_total, grid.j_total, grid.k_total) {}

    void save(const Config& config, const Time<Real>& time) const {
        std::ofstream ofs(config.save_dir + "/mhd."+ util::zfill(time.n_output, time.n_output_digits) +".bin", std::ios::binary);
        assert(ofs.is_open());
    
        auto write_array = [&ofs](const Array3D<Real>& arr) {
            ofs.write(reinterpret_cast<const char*>(arr.data()), sizeof(Real) * arr.size());
        };
        write_array(qq.ro);
        write_array(qq.vx);
        write_array(qq.vy);
        write_array(qq.vz);
        write_array(qq.bx);
        write_array(qq.by);
        write_array(qq.bz);
        write_array(qq.ei);
        write_array(qq.ph);
        ofs.close();
    
    };
};