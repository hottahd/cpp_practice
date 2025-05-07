#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>

#include "config.hpp"
#include "model.hpp"
#include "time_integrator.hpp"
#include "utility.hpp"

using util::pow2;

template <typename Real>
struct Init {
    Real xm;
    Real rol, prl, vvl;
    Real ror, prr, vvr;
};

void initial_condition(Model<Real>& model, const Init<Real>& init);

int main() {
    Model<Real> model = Model<Real>::from_config_file("../config/config.json");
    model.save_metadata();

    Init<Real> init;
    init.rol = 1.0  ; init.prl = 1.0; init.vvl = 0.0;
    init.ror = 0.125; init.prr = 0.1; init.vvr = 0.0;

    initial_condition(model, init);
    model.save_state();

    TimeIntegrator<Real> time_integrator(model);
    time_integrator.run();
    return 0;
}

void initial_condition(Model<Real>& model, const Init<Real>& init) {
    MHDCore<Real>& qq = model.mhd.qq;
    const auto& grid = model.grid;
    const auto& eos = model.eos;

    for (int i = 0; i < grid.i_total; ++i) {
        for (int j = 0; j < grid.j_total; ++j) {
            for (int k = 0; k < grid.k_total; ++k) {

                qq.vy(i, j, k) = 0.0;
                qq.vz(i, j, k) = 0.0;
                qq.bx(i, j, k) = 0.0;
                qq.by(i, j, k) = 0.0;
                qq.bz(i, j, k) = 0.0;
                qq.ph(i, j, k) = 0.0;

                qq.vx(i, j, k) = 0.0;
                qq.ei(i, j, k) = 1.0;
                qq.ro(i, j, k) = 1.0 + 0.3*std::exp(-pow2((grid.x[i] - 0.5) / 0.1));

                if (grid.x[i] < 0.5) {
                    qq.ro(i, j, k) = init.rol;
                    qq.ei(i, j, k) = init.prl / (eos.gm - 1.0) / qq.ro(i, j, k);
                    qq.vx(i, j, k) = init.vvl;
                } else {
                    qq.ro(i, j, k) = init.ror;
                    qq.ei(i, j, k) = init.prr / (eos.gm - 1.0) / qq.ro(i, j, k);
                    qq.vx(i, j, k) = init.vvr;
                }

                

            }
        }
    }
}