#pragma once

#include <cmath>
#include <algorithm>
#include <initializer_list>

#include "model.hpp"
#include "constants.hpp"
#include "boundary_condition.hpp"
#include "artificial_viscosity.hpp"

template <typename Real>
inline Real space_centered_4th (const Array3D<Real>& qq, Real dxyzi, int i, int j, int k, int is, int js, int ks) {
    return (
        -     qq(i + 2*is, j + 2*js, k + 2*ks)
        + 8.0*qq(i +   is, j +   js, k +   ks)
        - 8.0*qq(i -   is, j -   js, k -   ks)
        +     qq(i - 2*is, j - 2*js, k - 2*ks)
    )*inv12<Real>*dxyzi;
};

template <typename Real>
inline Real space_centered_4th (const Array3D<Real>& qq1, const Array3D<Real>& qq2, Real dxyzi, int i, int j, int k, int is, int js, int ks) {
    return (
        -     qq1(i + 2*is, j + 2*js, k + 2*ks)*qq2(i + 2*is, j + 2*js, k + 2*ks)
        + 8.0*qq1(i +   is, j +   js, k +   ks)*qq2(i +   is, j +   js, k +   ks)
        - 8.0*qq1(i -   is, j -   js, k -   ks)*qq2(i -   is, j -   js, k -   ks)
        +     qq1(i - 2*is, j - 2*js, k - 2*ks)*qq2(i - 2*is, j - 2*js, k - 2*ks)
    )*inv12<Real>*dxyzi;
};

template <typename Real>
inline Real space_centered_4th (const Array3D<Real>& qq1, const Array3D<Real>& qq2, const Array3D<Real>& qq3, Real dxyzi, int i, int j, int k, int is, int js, int ks) {
    return (
        -     qq1(i + 2*is, j + 2*js, k + 2*ks)*qq2(i + 2*is, j + 2*js, k + 2*ks)*qq3(i + 2*is, j + 2*js, k + 2*ks)
        + 8.0*qq1(i +   is, j +   js, k +   ks)*qq2(i +   is, j +   js, k +   ks)*qq3(i +   is, j +   js, k +   ks)
        - 8.0*qq1(i -   is, j -   js, k -   ks)*qq2(i -   is, j -   js, k -   ks)*qq3(i -   is, j -   js, k -   ks)
        +     qq1(i - 2*is, j - 2*js, k - 2*ks)*qq2(i - 2*is, j - 2*js, k - 2*ks)*qq3(i - 2*is, j - 2*js, k - 2*ks)
    )*inv12<Real>*dxyzi;
};

template <typename Real>
struct TimeIntegrator {
    // Disallow copying and assignment since this class manages resources
    TimeIntegrator(const TimeIntegrator&) = delete;
    TimeIntegrator& operator=(const TimeIntegrator&) = delete;

    Model<Real>& model;
    Config& config;
    Time<Real>& time;
    Grid<Real>& grid;
    EOS<Real>& eos;
    MHD<Real>& mhd;

    BoundaryCondition<Real> bc;
    ArtificialViscosity<Real> artdiff;

    Array3D<Real> pr, bb, ht, vb;
    Real cfl_number;
    
    TimeIntegrator(Model<Real>& model_)
        : model(model_),
          config(model_.config),
          time(model_.time),
          grid(model_.grid),
          eos(model_.eos),
          mhd(model_.mhd),
          bc(model_),
          artdiff(model_),
          pr(grid.i_total, grid.j_total, grid.k_total),
          bb(grid.i_total, grid.j_total, grid.k_total),
          ht(grid.i_total, grid.j_total, grid.k_total),
          vb(grid.i_total, grid.j_total, grid.k_total) {
            cfl_number = config.json_obj.at("time_integrator").at("cfl_number").get<Real>();
          }
        
    // core function for MHD time integration
    void update_sc4(MHDCore<Real>& qq_orgn, MHDCore<Real>& qq_argm, MHDCore<Real>& qq_rslt, Real dt) {
        for (int i = 0; i < grid.i_total; ++i) {
            for (int j = 0; j < grid.j_total; ++j) {
                for (int k = 0; k < grid.k_total; ++k) {
                    // gas pressure
                    pr(i, j, k) = qq_argm.ro(i,j,k)*qq_argm.ei(i,j,k)*(eos.gm - 1.0);
                    // squared magnetic strength
                    bb(i, j, k) = qq_argm.bx(i,j,k)*qq_argm.bx(i,j,k) 
                                + qq_argm.by(i,j,k)*qq_argm.by(i,j,k) 
                                + qq_argm.bz(i,j,k)*qq_argm.bz(i,j,k);
                    // enthalpy + 2*magnetic energy + kinetic energy      
                    ht(i, j, k) = 
                        + qq_argm.ro(i,j,k)*qq_argm.ei(i,j,k) + pr(i,j,k)
                        + bb(i,j,k)*pii4<Real>
                        + 0.5*qq_argm.ro(i,j,k)*(
                            + qq_argm.vx(i,j,k)*qq_argm.vx(i,j,k)
                            + qq_argm.vy(i,j,k)*qq_argm.vy(i,j,k)
                            + qq_argm.vz(i,j,k)*qq_argm.vz(i,j,k)
                        );
                    // v dot b
                    vb(i, j, k) = 
                        + qq_argm.vx(i,j,k)*qq_argm.bx(i,j,k)
                        + qq_argm.vy(i,j,k)*qq_argm.by(i,j,k)
                        + qq_argm.vz(i,j,k)*qq_argm.bz(i,j,k);
                }
            }
        }


        for (int i = grid.i_margin; i < grid.i_total-grid.i_margin; ++i) {
            const Real dxi = grid.dxi[i];
            for (int j = grid.j_margin; j < grid.j_total-grid.j_margin; ++j) {
                const Real dyi = grid.dyi[j];
                for (int k = grid.k_margin; k < grid.k_total-grid.k_margin; ++k) {
                    const Real dzi = grid.dzi[k];

                    // equation of continuity
                    qq_rslt.ro(i, j, k) = qq_orgn.ro(i, j, k) + dt * (
                    -space_centered_4th(qq_argm.ro, qq_argm.vx, dxi, i, j, k, grid.is, 0, 0)
                    -space_centered_4th(qq_argm.ro, qq_argm.vx, dyi, i, j, k, 0, grid.js, 0)
                    -space_centered_4th(qq_argm.ro, qq_argm.vx, dzi, i, j, k, 0, 0, grid.ks)
                    );

                    // x equation of motion
                    qq_rslt.vx(i, j, k) = (
                        qq_orgn.ro(i, j, k) * qq_orgn.vx(i, j, k) + dt * (
                        -space_centered_4th(qq_argm.ro, qq_argm.vx, qq_argm.vx, dxi, i, j, k, grid.is, 0, 0)
                        -space_centered_4th(qq_argm.ro, qq_argm.vx, qq_argm.vy, dyi, i, j, k, 0, grid.js, 0)
                        -space_centered_4th(qq_argm.ro, qq_argm.vx, qq_argm.vz, dzi, i, j, k, 0, 0, grid.ks)
                        -space_centered_4th(pr, dxi, i, j, k, grid.is, 0, 0)
                        -space_centered_4th(bb, dxi, i, j, k, grid.is, 0, 0)*pii8<Real>
                        +space_centered_4th(qq_argm.bx, qq_argm.bx, dxi, i, j, k, grid.is, 0, 0)*pii4<Real>
                        +space_centered_4th(qq_argm.bx, qq_argm.by, dyi, i, j, k, 0, grid.js, 0)*pii4<Real>
                        +space_centered_4th(qq_argm.bx, qq_argm.bz, dzi, i, j, k, 0, 0, grid.ks)*pii4<Real>
                        )
                    )/qq_rslt.ro(i, j, k);

                    // y equation of motion
                    qq_rslt.vy(i, j, k) = (
                        qq_orgn.ro(i, j, k) * qq_orgn.vy(i, j, k) + dt * (
                        -space_centered_4th(qq_argm.ro, qq_argm.vy, qq_argm.vx, dxi, i, j, k, grid.is, 0, 0)
                        -space_centered_4th(qq_argm.ro, qq_argm.vy, qq_argm.vy, dyi, i, j, k, 0, grid.js, 0)
                        -space_centered_4th(qq_argm.ro, qq_argm.vy, qq_argm.vz, dzi, i, j, k, 0, 0, grid.ks)
                        -space_centered_4th(pr, dyi, i, j, k, 0, grid.js, 0)
                        -space_centered_4th(bb, dyi, i, j, k, 0, grid.js, 0)*pii8<Real>
                        +space_centered_4th(qq_argm.by, qq_argm.bx, dxi, i, j, k, grid.is, 0, 0)*pii4<Real>
                        +space_centered_4th(qq_argm.by, qq_argm.by, dyi, i, j, k, 0, grid.js, 0)*pii4<Real>
                        +space_centered_4th(qq_argm.by, qq_argm.bz, dzi, i, j, k, 0, 0, grid.ks)*pii4<Real>
                        )
                    )/qq_rslt.ro(i, j, k);

                    // z equation of motion
                    qq_rslt.vz(i, j, k) = (
                        qq_orgn.ro(i, j, k) * qq_orgn.vz(i, j, k) + dt * (
                        -space_centered_4th(qq_argm.ro, qq_argm.vz, qq_argm.vx, dxi, i, j, k, grid.is, 0, 0)
                        -space_centered_4th(qq_argm.ro, qq_argm.vz, qq_argm.vy, dyi, i, j, k, 0, grid.js, 0)
                        -space_centered_4th(qq_argm.ro, qq_argm.vz, qq_argm.vz, dzi, i, j, k, 0, 0, grid.ks)
                        -space_centered_4th(pr, dzi, i, j, k, 0, 0, grid.ks)
                        -space_centered_4th(bb, dzi, i, j, k, 0, 0, grid.ks)*pii8<Real>
                        +space_centered_4th(qq_argm.bz, qq_argm.bx, dxi, i, j, k, grid.is, 0, 0)*pii4<Real>
                        +space_centered_4th(qq_argm.bz, qq_argm.by, dyi, i, j, k, 0, grid.js, 0)*pii4<Real>
                        +space_centered_4th(qq_argm.bz, qq_argm.bz, dzi, i, j, k, 0, 0, grid.ks)*pii4<Real>
                        )
                    )/qq_rslt.ro(i, j, k);

                    // x magnetic induction
                    qq_rslt.bx(i, j, k) = qq_orgn.bx(i, j, k) + dt * (
                        -space_centered_4th(qq_argm.vy, qq_argm.bx, dyi, i, j, k, 0, grid.js, 0)
                        +space_centered_4th(qq_argm.vx, qq_argm.by, dyi, i, j, k, 0, grid.js, 0)
                        -space_centered_4th(qq_argm.vz, qq_argm.bx, dzi, i, j, k, 0, 0, grid.ks)
                        +space_centered_4th(qq_argm.vx, qq_argm.bz, dzi, i, j, k, 0, 0, grid.ks)
                    );

                    // y magnetic induction
                    qq_rslt.by(i, j, k) = qq_orgn.by(i, j, k) + dt * (
                        -space_centered_4th(qq_argm.vx, qq_argm.by, dxi, i, j, k, grid.is, 0, 0)
                        +space_centered_4th(qq_argm.vy, qq_argm.bx, dxi, i, j, k, grid.is, 0, 0)
                        -space_centered_4th(qq_argm.vz, qq_argm.by, dzi, i, j, k, 0, 0, grid.ks)
                        +space_centered_4th(qq_argm.vy, qq_argm.bz, dzi, i, j, k, 0, 0, grid.ks)
                    );

                    // z magnetic induction
                    qq_rslt.bz(i, j, k) = qq_orgn.bz(i, j, k) + dt * (
                        -space_centered_4th(qq_argm.vx, qq_argm.bz, dxi, i, j, k, grid.is, 0, 0)
                        +space_centered_4th(qq_argm.vz, qq_argm.bx, dxi, i, j, k, grid.is, 0, 0)
                        -space_centered_4th(qq_argm.vy, qq_argm.bz, dyi, i, j, k, 0, grid.js, 0)
                        +space_centered_4th(qq_argm.vz, qq_argm.by, dyi, i, j, k, 0, grid.js, 0)
                    );

                // Et: total energy per unit volume
                // ei: is the internal energy per unit mass 
                    const Real Et = 
                        + qq_orgn.ro(i,j,k)*qq_orgn.ei(i,j,k)
                        + 0.5*qq_orgn.ro(i, j, k)*(
                            + qq_orgn.vx(i, j, k)*qq_orgn.vx(i, j, k)
                            + qq_orgn.vy(i, j, k)*qq_orgn.vy(i, j, k)
                            + qq_orgn.vz(i, j, k)*qq_orgn.vz(i, j, k)
                        )
                        + pii8<Real>*(
                            + qq_orgn.bx(i,j,k)*qq_orgn.bx(i,j,k)
                            + qq_orgn.by(i,j,k)*qq_orgn.by(i,j,k)
                            + qq_orgn.bz(i,j,k)*qq_orgn.bz(i,j,k)
                        )
                        ;

                    qq_rslt.ei(i, j, k) = ( Et + dt * (
                        -space_centered_4th(ht, qq_argm.vx, dxi, i, j, k, grid.is, 0, 0)
                        -space_centered_4th(ht, qq_argm.vy, dyi, i, j, k, 0, grid.js, 0)
                        -space_centered_4th(ht, qq_argm.vz, dzi, i, j, k, 0, 0, grid.ks)
                        -space_centered_4th(vb, qq_argm.bx, dxi, i, j, k, grid.is, 0, 0)*pii4<Real>
                        -space_centered_4th(vb, qq_argm.by, dyi, i, j, k, 0, grid.js, 0)*pii4<Real>
                        -space_centered_4th(vb, qq_argm.bz, dzi, i, j, k, 0, 0, grid.ks)*pii4<Real>  )
                        -0.5*qq_rslt.ro(i,j,k)*(
                            + qq_rslt.vx(i,j,k)*qq_rslt.vx(i,j,k)
                            + qq_rslt.vy(i,j,k)*qq_rslt.vy(i,j,k)
                            + qq_rslt.vz(i,j,k)*qq_rslt.vz(i,j,k) )
                        - pii8<Real>*(
                            + qq_rslt.bx(i,j,k)*qq_rslt.bx(i,j,k)
                            + qq_rslt.by(i,j,k)*qq_rslt.by(i,j,k)
                            + qq_rslt.bz(i,j,k)*qq_rslt.bz(i,j,k) )
                    )/qq_rslt.ro(i, j, k);

                }
            }
        }
    }

    void runge_kutta_4step(){
        MHDCore<Real>& qq     = mhd.qq;
        MHDCore<Real>& qq_argm = mhd.qq_argm;
        MHDCore<Real>& qq_rslt = mhd.qq_rslt;

        update_sc4(qq     , qq     , qq_rslt, time.dt/4.0);
        qq_argm.copy_from(qq_rslt);
        bc.apply(qq_argm);
        update_sc4(qq     , qq_argm, qq_rslt, time.dt/3.0);
        qq_argm.copy_from(qq_argm);
        bc.apply(qq_argm);
        update_sc4(qq     , qq_argm, qq_rslt, time.dt/2.0);
        qq_argm.copy_from(qq_rslt);
        bc.apply(qq_argm);
        update_sc4(qq     , qq_argm, qq_rslt, time.dt     );
        qq.copy_from(qq_rslt);
        bc.apply(qq);
    }

    void cfl_condition() {

        this->time.dt = 1.e10;

        for (int i = grid.i_margin; i < grid.i_total-grid.i_margin; ++i) {
            for (int j = grid.j_margin; j < grid.j_total-grid.j_margin; ++j) {
                for (int k = grid.k_margin; k < grid.k_total-grid.k_margin; ++k) {
                    // cs: sound speed, vv: fluid velocity, ca: Alfvén speed
                    Real cs = std::sqrt(this->eos.gm*(this->eos.gm-1.0)*this->mhd.qq.ei(i,j,k));
                    Real vv = std::sqrt(
                        + this->mhd.qq.vx(i,j,k)*this->mhd.qq.vx(i,j,k)
                        + this->mhd.qq.vy(i,j,k)*this->mhd.qq.vy(i,j,k)
                        + this->mhd.qq.vz(i,j,k)*this->mhd.qq.vz(i,j,k)
                    );  
                    Real ca = std::sqrt( (
                        + this->mhd.qq.bx(i,j,k)*this->mhd.qq.bx(i,j,k)
                        + this->mhd.qq.by(i,j,k)*this->mhd.qq.by(i,j,k)
                        + this->mhd.qq.bz(i,j,k)*this->mhd.qq.bz(i,j,k)
                    )/this->mhd.qq.ro(i,j,k)*pii4<Real>);
                    this->time.dt = std::min(this->time.dt, 
                            cfl_number*std::min<Real>({this->grid.dx[i], this->grid.dy[j], this->grid.dz[k]})/(cs + vv));
                }
            }
        }

    }

    void run() {
        // Time integration loop
        this->time.dt = 0.1;
        model.save_if_needed();
        while (this->time.time < this->time.tend) {

            // basic MHD time integration
            cfl_condition();
            runge_kutta_4step();

            // artificial vicosities
            artdiff.characteristic_velocity_eval();
            artdiff.update(mhd.qq, mhd.qq_rslt, artdiff.cc, grid.dxi, "x");
            bc.apply(mhd.qq_rslt);
            mhd.qq.copy_from(mhd.qq_rslt);
            artdiff.update(mhd.qq, mhd.qq_rslt, artdiff.cc, grid.dyi, "y");
            bc.apply(mhd.qq_rslt);
            mhd.qq.copy_from(mhd.qq_rslt);
            artdiff.update(mhd.qq, mhd.qq_rslt, artdiff.cc, grid.dzi, "z");
            bc.apply(mhd.qq_rslt);
            mhd.qq.copy_from(mhd.qq_rslt);

            // Time is update after all procedures
            this->time.update();
            model.save_if_needed();
        }
    }
};

