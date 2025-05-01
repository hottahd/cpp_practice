#include "mhd.hpp"
#include "utility.hpp"
#include "time.hpp"
#include <cmath>

template <typename Real>
void MHD<Real>::initial_condition(const Grid<Real>& grid) {
    Real dd = 0.1; // width of Gaussian 
    Real xm = 0.5*(grid.xmax + grid.xmin); // center of Gaussian
    
    // for (int i = 0; i < grid.i_total; ++i) {
    //     q0[i] = std::exp(-std::pow((grid.x[i] - xm) / dd, 2));
    // }
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