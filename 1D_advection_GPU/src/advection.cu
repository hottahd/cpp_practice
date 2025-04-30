#include "advection.cup"
#include "grid.hpp"
#include "iostream"
#include "quantity.hpp"
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cstdio>

template<typename Real>
__device__ void update_core_kernel(Real* q0, Real* q1, Real dx, Real dt, Real vc, int i) {
    q1[i] = q0[i] - vc*(q0[i + 1] - q0[i - 1]) / (2.0 * dx) * dt;
}

template<typename Real>
__global__ void update_kernel(Real* q0, Real* q1, int i_total, Real dx, Real dt, Real vc) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i>=1 && i <= i_total-2) {
        update_core_kernel(q0, q1, dx, dt, vc, i);
    }
}

template<typename Real>
__global__ void bc_kernel(Real* qq, int i_total, int margin) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < margin) {
        qq[i] = qq[i_total - 2*margin + i];
        qq[i_total - margin + i] = qq[margin + i];
    }
}

template <typename Real>
void Advection<Real>::cfl_condition() {
    Real safety = 0.5; // safety factor

    time.dt = 1.e10;
    for (int i = 0; i < grid.i_total; ++i) {
        time.dt = std::min(time.dt, safety * grid.dx / std::abs(vc));
    }
};

template <typename Real>
void Advection<Real>::update() {
    cfl_condition();
    sc2ssprk();
}

template <typename Real>
void Advection<Real>::io_step() {
    if (time.time >= time.n_output * time.dt_output) {
        quantity.save(config, time);
        std::cout
        << "time = " 
        << std::fixed << std::setprecision(2) << std::setw(8) << time.time
        << ";  n_step = " 
        << std::setw(8) << time.n_step
        << ";  n_output = " 
        << std::setw(8) << time.n_output
        << std::endl;

        time.n_output++;
    }
}


template <typename Real>
void Advection<Real>::sc2ssprk() {

    quantity.allocate_device(grid);
    quantity.upload_q0();
    
    int block_size = 256;
    int grid_size = (grid.i_total + block_size - 1) / block_size;
    update_kernel<<<grid_size, block_size>>>(quantity.q0_dev, quantity.q1_dev, grid.i_total, grid.dx, time.dt, vc);
    bc_kernel<<<grid_size, block_size>>>(quantity.q1_dev, grid.i_total, grid.margin);

    update_kernel<<<grid_size, block_size>>>(quantity.q1_dev, quantity.q2_dev, grid.i_total, grid.dx, time.dt, vc);
    bc_kernel<<<grid_size, block_size>>>(quantity.q2_dev, grid.i_total, grid.margin);
    cudaDeviceSynchronize();

    cudaMemcpy(quantity.q1.data(), quantity.q1_dev, sizeof(Real) * grid.i_total, cudaMemcpyDeviceToHost);
    cudaMemcpy(quantity.q2.data(), quantity.q2_dev, sizeof(Real) * grid.i_total, cudaMemcpyDeviceToHost);
    
    for ( int i = 0; i < grid.i_total; ++i) {
        quantity.q0[i] = 0.5*( quantity.q0[i] + quantity.q2[i] );
    }

    quantity.free_device();

}

// template<typename Real>
// void Advection<Real>::sc2(std::vector<Real>& qq, std::vector<Real>& dqq, Grid<Real>& grid, Real& dt, Real& vc) {

//     for (int i = 1; i < grid.i_total - 1; ++i) {
//         dqq[i] = -vc * (qq[i + 1] - qq[i - 1]) / (2.0 * grid.dx)*dt;
//     };
// }

template<typename Real>
void Advection<Real>::bc(std::vector<Real>& qq, Grid<Real>& grid) {
    for (int i = 0; i < grid.margin; ++i) {
        qq[i] = qq[grid.i_total - 2*grid.margin + i];
        qq[grid.i_total - grid.margin + i] = qq[grid.margin + i];
    }
}

template struct Advection<double>;
template struct Advection<float>;