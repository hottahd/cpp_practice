#include "advection.cuh"

#include <vector>
#include <cstdlib>
#include <iostream>
#include <cstdio>

#include "grid.cuh"
#include "quantity.cuh"


template<typename Real>
__global__ void update_kernel(Real* q0, Real* q1, Grid_Device<Real> grid, Real dt, Real vc) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i>=1 && i <= grid.i_total-2) {
        q1[i] = q0[i] - vc*(q0[i + 1] - q0[i - 1]) / (2.0 * grid.dx) * dt;
    }
}

template<typename Real>
__global__ void bc_kernel(Real* qq, Grid_Device<Real> grid) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < grid.margin) {
        qq[i] = qq[grid.i_total - 2*grid.margin + i];
        qq[grid.i_total - grid.margin + i] = qq[grid.margin + i];
    }
}

template<typename Real>
__global__ void final_kernel(Real* q0, Real* q2, Grid_Device<Real> grid) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < grid.i_total) {
        q0[i] = 0.5*(q0[i] + q2[i]);
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
void Advection<Real>::run() {

    cfl_condition();

    quantity.allocate_device(grid);
    quantity.upload_q0();
    
    while (time.time < time.tend) {
        update();
        time.update();
        io_step();
    };
    quantity.free_device();

}


template <typename Real>
void Advection<Real>::update() {
    sc2ssprk();
}

template <typename Real>
void Advection<Real>::io_step() {
    if (time.time >= time.n_output * time.dt_output) {
        quantity.download_q0();
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
    
    int block_size = 256;
    int grid_size = (grid.i_total + block_size - 1) / block_size;
    update_kernel<<<grid_size, block_size>>>(quantity.q0_dev, quantity.q1_dev, grid_device, time.dt, vc);
    bc_kernel<<<grid_size, block_size>>>(quantity.q1_dev, grid_device);

    update_kernel<<<grid_size, block_size>>>(quantity.q1_dev, quantity.q2_dev, grid_device, time.dt, vc);
    bc_kernel<<<grid_size, block_size>>>(quantity.q2_dev, grid_device);
    
    final_kernel<<<grid_size, block_size>>>(quantity.q0_dev, quantity.q2_dev, grid_device);

}

template<typename Real>
void Advection<Real>::bc(std::vector<Real>& qq, Grid<Real>& grid) {
    for (int i = 0; i < grid.margin; ++i) {
        qq[i] = qq[grid.i_total - 2*grid.margin + i];
        qq[grid.i_total - grid.margin + i] = qq[grid.margin + i];
    }
}

template struct Advection<double>;
template struct Advection<float>;