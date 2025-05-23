#include "time.hpp"

template <typename Real>
void Time<Real>::update() {
    time += dt;
    n_step++;
};

template struct Time<double>;
template struct Time<float>;