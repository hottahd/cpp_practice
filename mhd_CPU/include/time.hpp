#pragma once
#include <cassert>

template <typename Real>
struct Time{
    Real time, tend;
    Real dt_output;
    Real dt;
    int n_step, n_output, n_output_digits;

    Time(Real tend_, Real dt_output_, int n_output_digits_ = 8)
        :tend(tend_), dt_output(dt_output_), n_output_digits(n_output_digits_) {
            assert(tend > 0);
            assert(dt_output > 0);

            time = 0;
            n_step = 0;
            n_output = 0;
        }

    void update();
};