#include "advection.hpp"
#include "grid.hpp"
#include "iostream"
#include <vector>

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
    std::vector<Real> dqq(grid.i_total);
    sc2(quantity.q0, dqq, grid, time.dt, vc);
    for (int i = 1; i < grid.i_total - 1; ++i) {
        quantity.q1[i] = quantity.q0[i] + dqq[i];
    }
    bc(quantity.q1, grid);
    
    sc2(quantity.q1, dqq, grid, time.dt, vc);
    for (int i = 1; i < grid.i_total - 1; ++i) {
        quantity.q2[i] = 0.5*(quantity.q0[i] + quantity.q1[i] + dqq[i]);
    }
    bc(quantity.q2, grid);

    for ( int i = 0; i < grid.i_total; ++i) {
        quantity.q0[i] = quantity.q2[i];
    }
}

template<typename Real>
void Advection<Real>::sc2(std::vector<Real>& qq, std::vector<Real>& dqq, Grid<Real>& grid, Real& dt, Real& vc) {

    for (int i = 1; i < grid.i_total - 1; ++i) {
        dqq[i] = -vc * (qq[i + 1] - qq[i - 1]) / (2.0 * grid.dx)*dt;
    };
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