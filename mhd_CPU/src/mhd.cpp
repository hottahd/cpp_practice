#include "mhd.hpp"

#include <cmath>
#include <omp.h>

#include "utility.hpp"
#include "time.hpp"

template <typename Real>
void MHD<Real>::initial_condition(const Grid<Real>& grid) {
    double pr; // pressure
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < grid.i_total; ++i) {
        for (int j = 0; j < grid.j_total; ++j) {
            for (int k = 0; k < grid.k_total; ++k) {
                qq.ro(i,j,k) = 0.0;
                qq.vx(i,j,k) = 0.0;
                qq.vy(i,j,k) = 0.0;
                qq.vz(i,j,k) = 0.0;
                qq.bx(i,j,k) = 0.0;
                qq.by(i,j,k) = 0.0;
                qq.bz(i,j,k) = 0.0;
                qq.ei(i,j,k) = 0.0;
                qq.ph(i,j,k) = 0.0;

                if (grid.x[i] < 0.5) {
                    qq.ro(i,j,k) = 1.0;
                    pr = 1.0;
                    qq.vx(i,j,k) = 0.0;
                } else {
                    qq.ro(i,j,k) = 0.125;
                    pr = 0.1;
                    qq.vx(i,j,k) = 0.0;
                }

            }
        }
    }
}

template <typename Real>
void MHD<Real>::save(const Config& config, const Time<Real>& time) const {
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
}

template struct MHD<double>;
template struct MHD<float>;