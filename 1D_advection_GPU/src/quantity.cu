#include "quantity.hpp"
#include "utility.hpp"
#include "time.hpp"
#include <cmath>

template<typename Real>
void Quantity<Real>::allocate_device(const Grid<Real>& grid) {
    cudaMalloc(&q0_dev,  sizeof(Real) * grid.i_total);
    cudaMalloc(&q1_dev,  sizeof(Real) * grid.i_total);
    cudaMalloc(&q2_dev,  sizeof(Real) * grid.i_total);
    cudaMalloc(&qq_dev,  sizeof(Real) * grid.i_total);
    cudaMalloc(&dqq_dev, sizeof(Real) * grid.i_total);
}

template<typename Real>
void Quantity<Real>::upload_q0() {
    cudaMemcpy(q0_dev, q0.data(), sizeof(Real) * q0.size(), cudaMemcpyHostToDevice);
}

template<typename Real>
void Quantity<Real>::download_q0() {
    cudaMemcpy(q0.data(), q0_dev, sizeof(Real) * q0.size(), cudaMemcpyDeviceToHost);
}

template<typename Real>
void Quantity<Real>::free_device() {
    cudaFree(q1_dev);
    cudaFree(q2_dev);
    cudaFree(qq_dev);
    cudaFree(dqq_dev);
}

template <typename Real>
void Quantity<Real>::initial_condition(const Grid<Real>& grid) {
    Real dd = 0.1; // width of Gaussian 
    Real xm = 0.5*(grid.xmax + grid.xmin); // center of Gaussian
    
    for (int i = 0; i < grid.i_total; ++i) {
        q0[i] = std::exp(-std::pow((grid.x[i] - xm) / dd, 2));
    }
}

template <typename Real>
void Quantity<Real>::save(const Config& config, const Time<Real>& time) const {
    std::ofstream ofs(config.save_dir + "/quantity."+ util::zfill(time.n_output, time.n_output_digits) +".bin", std::ios::binary);
    assert(ofs.is_open());
    ofs.write(reinterpret_cast<const char*>(q0.data()), sizeof(Real)*q0.size());
}

template struct Quantity<double>;
template struct Quantity<float>;