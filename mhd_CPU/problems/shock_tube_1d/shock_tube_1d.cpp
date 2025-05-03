#include "initial_condition.hpp"

template <typename Real>
void ShockTube_1D_IC<Real>::operator()(MHD<Real>& mhd, const Grid<Real>& grid) const{

}
template struct ShockTube_1D_IC<double>;
template struct ShockTube_1D_IC<float>;
