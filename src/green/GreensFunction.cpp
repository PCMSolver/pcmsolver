#include "GreensFunction.hpp"

#include <iostream>
#include <cmath>
#include <stdexcept>

#include "Config.hpp"

#include "EigenPimpl.hpp"
#include "TaylorPimpl.hpp"

#include "DerivativeTypes.hpp"

template <typename T>
void GreensFunction<T>::delta(double value)
{
    if (value <= 1.0e-10) {
        throw std::invalid_argument("Delta value must be larger than 1.0e-10");
    }
    delta_ = value;
}

template <typename T>
double GreensFunction<T>::function(const Eigen::Vector3d & source,
                                   const Eigen::Vector3d & probe) const
{
    T sp[3], pp[3], res;
    sp[0] = source(0);
    sp[1] = source(1);
    sp[2] = source(2);
    pp[0] = probe(0);
    pp[1] = probe(1);
    pp[2] = probe(2);
    res = evaluate(sp, pp);
    return res[0];
}

template <>
double GreensFunction<double>::function(const Eigen::Vector3d & source,
                                        const Eigen::Vector3d & probe) const
{
    double sp[3], pp[3], res;
    sp[0] = source(0);
    sp[1] = source(1);
    sp[2] = source(2);
    pp[0] = probe(0);
    pp[1] = probe(1);
    pp[2] = probe(2);
    res = evaluate(sp, pp);
    return res;
}

template <typename T>
double GreensFunction<T>::derivativeSource(const Eigen::Vector3d & direction,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    T t1[3], t2[3], derivative;
    //	direction.normalize();
    t1[0] = p1(0);
    t1[0][1] = direction(0);
    t1[1] = p1(1);
    t1[1][1] = direction(1);
    t1[2] = p1(2);
    t1[2][1] = direction(2);
    t2[0] = p2(0);
    t2[1] = p2(1);
    t2[2] = p2(2);
    derivative = evaluate(t1, t2);
    return derivative[1];
}

template <>
double GreensFunction<double>::derivativeSource(const Eigen::Vector3d & direction,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    Eigen::Vector3d deltaPlus  = p1 + direction * delta_ / direction.norm();
    Eigen::Vector3d deltaMinus = p1 - direction * delta_ / direction.norm();
    double funcPlus  = function(deltaPlus,  p2);
    double funcMinus = function(deltaMinus, p2);
    return (funcPlus - funcMinus)/(2.0*delta_);
}

template <typename T>
double GreensFunction<T>::derivativeProbe(const Eigen::Vector3d & direction,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    T t1[3], t2[3], derivative;
    //	direction.normalize();
    t1[0] = p1(0);
    t1[1] = p1(1);
    t1[2] = p1(2);
    t2[0] = p2(0);
    t2[0][1] = direction(0);
    t2[1] = p2(1);
    t2[1][1] = direction(1);
    t2[2] = p2(2);
    t2[2][1] = direction(2);
    derivative = evaluate(t1, t2);
    return derivative[1];
}

template <>
double GreensFunction<double>::derivativeProbe(const Eigen::Vector3d & direction,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    Eigen::Vector3d deltaPlus  = p2 + direction * delta_ / direction.norm();
    Eigen::Vector3d deltaMinus = p2 - direction * delta_ / direction.norm();
    double funcPlus  = function(p1, deltaPlus);
    double funcMinus = function(p1, deltaMinus);
    return (funcPlus - funcMinus)/(2.0*delta_);
}

template <typename T>
Eigen::Vector3d GreensFunction<T>::gradientSource(const Eigen::Vector3d & p1,
        const Eigen::Vector3d & p2) const
{
    Eigen::Vector3d g;
    gradientSource(g, p1, p2);
    return g;
}

template <typename T>
void GreensFunction<T>::gradientSource(Eigen::Vector3d & g,
                                       const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    T t1[3], t2[3], grad;
    t1[0] = p1(0);
    t1[0][1] = 1;
    t1[1] = p1(1);
    t1[1][2] = 1;
    t1[2] = p1(2);
    t1[2][3] = 1;
    t2[0] = p2(0);
    t2[1] = p2(1);
    t2[2] = p2(2);
    grad = evaluate(t1, t2);
    g << grad[1], grad[2], grad[3];
    return;
}

template <>
void GreensFunction<double>::gradientSource(Eigen::Vector3d & g,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    Eigen::Vector3d direction;
    direction << 1.0, 0.0, 0.0;
    g(0) = derivativeSource(direction, p1, p2);
    direction << 0.0, 1.0, 0.0;
    g(1) = derivativeSource(direction, p1, p2);
    direction << 0.0, 0.0, 1.0;
    g(2) = derivativeSource(direction, p1, p2);
}

template <>
void GreensFunction<AD_directional>::gradientSource(Eigen::Vector3d & g,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    Eigen::Vector3d direction;
    direction << 1.0, 0.0, 0.0;
    g(0) = derivativeSource(direction, p1, p2);
    direction << 0.0, 1.0, 0.0;
    g(1) = derivativeSource(direction, p1, p2);
    direction << 0.0, 0.0, 1.0;
    g(2) = derivativeSource(direction, p1, p2);
}

template <typename T>
Eigen::Vector3d GreensFunction<T>::gradientProbe(const Eigen::Vector3d & p1,
        const Eigen::Vector3d & p2) const
{
    Eigen::Vector3d g;
    gradientProbe(g, p1, p2);
    return g;
}

template <typename T>
void GreensFunction<T>::gradientProbe(Eigen::Vector3d & g,
                                      const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    T t1[3], t2[3], grad;
    t1[0] = p1(0);
    t1[1] = p1(1);
    t1[2] = p1(2);
    t2[0] = p2(0);
    t2[0][1] = 1;
    t2[1] = p2(1);
    t2[1][2] = 1;
    t2[2] = p2(2);
    t2[2][3] = 1;
    grad = evaluate(t1, t2);
    g << grad[1], grad[2], grad[3];
    return;
}

template <>
void GreensFunction<double>::gradientProbe(Eigen::Vector3d & g,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    Eigen::Vector3d direction;
    direction << 1.0, 0.0, 0.0;
    g(0) = derivativeProbe(direction, p1, p2);
    direction << 0.0, 1.0, 0.0;
    g(1) = derivativeProbe(direction, p1, p2);
    direction << 0.0, 0.0, 1.0;
    g(2) = derivativeProbe(direction, p1, p2);
}

template <>
void GreensFunction<AD_directional>::gradientProbe(Eigen::Vector3d & g,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    Eigen::Vector3d direction;
    direction << 1.0, 0.0, 0.0;
    g(0) = derivativeProbe(direction, p1, p2);
    direction << 0.0, 1.0, 0.0;
    g(1) = derivativeProbe(direction, p1, p2);
    direction << 0.0, 0.0, 1.0;
    g(2) = derivativeProbe(direction, p1, p2);
}

template <typename T>
std::ostream & GreensFunction<T>::printObject(std::ostream & os)
{
    os << "Green's Function" << std::endl;
    os << "Delta = " << delta_ << std::endl;
    os << "Uniform = " << uniform_;
    return os;
}

template class GreensFunction<double>;
template class GreensFunction<AD_directional>;
template class GreensFunction<AD_gradient>;
template class GreensFunction<AD_hessian>;
