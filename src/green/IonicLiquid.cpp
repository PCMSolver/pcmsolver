#include "IonicLiquid.hpp"

#include <cmath>
#include <ostream>
#include <stdexcept>

#include "Config.hpp"

#include "EigenPimpl.hpp"
#include "TaylorPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "GreensFunction.hpp"
#include "IGreensFunction.hpp"

template<typename T>
double IonicLiquid<T>::derivative(const Eigen::Vector3d & direction,
                                        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    // NORMALIZATION TEMPORARILY REMOVED /direction.norm();
    return epsilon_ * (this->derivativeProbe(direction, p1, p2));  
}

template<typename T>
T IonicLiquid<T>::operator()(T * sp, T * pp) const
{
    T distance = sqrt((sp[0] - pp[0]) * (sp[0] - pp[0]) +
                      (sp[1] - pp[1]) * (sp[1] - pp[1]) +
                      (sp[2] - pp[2]) * (sp[2] - pp[2]));
    return (exp(-kappa_ * distance) / (epsilon_ * distance));
}

template <typename T>
std::ostream & IonicLiquid<T>::printObject(std::ostream & os)
{
    os << "Green's function type: ionic liquid" << std::endl;
    os << "Permittivity         = " << epsilon_ << std::endl;
    os << "Inverse Debye length = " << kappa_;
    return os;
}

template class IonicLiquid<double>;
template class IonicLiquid<AD_directional>;
template class IonicLiquid<AD_gradient>;
template class IonicLiquid<AD_hessian>;
