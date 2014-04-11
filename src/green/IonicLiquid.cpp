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
T IonicLiquid<T>::evaluate(T * sp, T * pp) const
{
    T distance = sqrt((sp[0] - pp[0]) * (sp[0] - pp[0]) +
                      (sp[1] - pp[1]) * (sp[1] - pp[1]) +
                      (sp[2] - pp[2]) * (sp[2] - pp[2]));
    return (exp(-kappa_ * distance) / (epsilon_ * distance));
}

template <typename T>
void IonicLiquid<T>::operator()(Eigen::MatrixXd & S, Eigen::MatrixXd & D,
                                      const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                                      const Eigen::VectorXd & areas, const Eigen::VectorXd & radii) const
{
    throw std::runtime_error("Green's function for an ionic liquid has not yet been implemented!");
}

template <typename T>
void IonicLiquid<T>::operator()(Eigen::MatrixXd & S,
                                      const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                                      const Eigen::VectorXd & areas) const
{
    throw std::runtime_error("Green's function for an ionic liquid has not yet been implemented!");
}

template <typename T>
void IonicLiquid<T>::operator()(Eigen::MatrixXd & D,
                                      const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                                      const Eigen::VectorXd & areas, const Eigen::VectorXd & radii) const
{
    throw std::runtime_error("Green's function for an ionic liquid has not yet been implemented!");
}

template<typename T>
double IonicLiquid<T>::compDiagonalElementS(double area) const
{
    throw std::runtime_error("Green's function for an ionic liquid has not yet been implemented!");
}

template<typename T>
double IonicLiquid<T>::compDiagonalElementD(double area, double radius) const
{
    throw std::runtime_error("Green's function for an ionic liquid has not yet been implemented!");
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
