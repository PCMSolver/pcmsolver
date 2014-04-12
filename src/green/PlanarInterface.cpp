#include "PlanarInterface.hpp"

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
double PlanarInterface<T>::derivative(const Eigen::Vector3d & direction,
                                        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    throw std::runtime_error("Green's function for a planar interface has not yet been implemented!");
}

template<typename T>
T PlanarInterface<T>::operator()(T * sp, T * pp) const
{
    throw std::runtime_error("Green's function for a planar interface has not yet been implemented!");
}

template <typename T>
std::ostream & PlanarInterface<T>::printObject(std::ostream & os)
{
    os << "Green's function type: planar interface" << std::endl;
    os << "Permittivity (layer 1) = " << eps1_ << std::endl;
    os << "Permittivity (layer 2) = " << eps2_ << std::endl;
    os << "Position               = " << pos_ << std::endl;
    os << "Width                  = " << width_;
    return os;
}

template class PlanarInterface<double>;
template class PlanarInterface<AD_directional>;
template class PlanarInterface<AD_gradient>;
template class PlanarInterface<AD_hessian>;
