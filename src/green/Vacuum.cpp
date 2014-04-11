#include "Vacuum.hpp"

#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"
#include "TaylorPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "GreensFunction.hpp"
#include "IGreensFunction.hpp"

static double factor = 1.07;

template<typename T>
double Vacuum<T>::derivative(const Eigen::Vector3d & direction,
                             const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    return this->derivativeProbe(direction, p1, p2);
    //    return direction.dot(g);  // NORMALIZTION TEMPORARY REMOVED /direction.norm();
}

template<typename T>
T Vacuum<T>::evaluate(T * sp, T * pp) const
{
    T res;
    res = 1.0/sqrt((sp[0]-pp[0])*(sp[0]-pp[0])+
                   (sp[1]-pp[1])*(sp[1]-pp[1])+
                   (sp[2]-pp[2])*(sp[2]-pp[2]));
    return res;
}

template <typename T>
void Vacuum<T>::operator()(Eigen::MatrixXd & S, Eigen::MatrixXd & D,
                           const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                           const Eigen::VectorXd & areas, const Eigen::VectorXd & radii) const
{
    int size = S.rows();

    for (int i = 0; i < size; ++i) {
        S(i, i) = factor * std::sqrt(4 * M_PI / areas(i));
        D(i, i) = -factor * std::sqrt(M_PI/ areas(i)) * (1.0 / radii(i));
        Eigen::Vector3d source = centers.col(i);
        for (int j = 0; j < size; ++j) {
            Eigen::Vector3d probe = centers.col(j);
            Eigen::Vector3d probeNormal = normals.col(j);
            probeNormal.normalize();
            if (i != j) {
                S(i, j) = this->function(source, probe);
                D(i, j) = this->derivative(probeNormal, source, probe);
            }
        }
    }
}

template <typename T>
void Vacuum<T>::operator()(Eigen::MatrixXd & S,
                           const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                           const Eigen::VectorXd & areas) const
{
    int size = S.rows();

    for (int i = 0; i < size; ++i) {
        S(i, i) = factor * std::sqrt(4 * M_PI / areas(i));
        Eigen::Vector3d source = centers.col(i);
        for (int j = 0; j < size; ++j) {
            Eigen::Vector3d probe = centers.col(j);
            Eigen::Vector3d probeNormal = normals.col(j);
            probeNormal.normalize();
            if (i != j) {
                S(i, j) = this->function(source, probe);
            }
        }
    }
}

template <typename T>
void Vacuum<T>::operator()(Eigen::MatrixXd & D,
                           const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                           const Eigen::VectorXd & areas, const Eigen::VectorXd & radii) const
{
    int size = D.rows();

    for (int i = 0; i < size; ++i) {
        D(i, i) = -factor * std::sqrt(M_PI/ areas(i)) * (1.0 / radii(i));
        Eigen::Vector3d source = centers.col(i);
        for (int j = 0; j < size; ++j) {
            Eigen::Vector3d probe = centers.col(j);
            Eigen::Vector3d probeNormal = normals.col(j);
            probeNormal.normalize();
            if (i != j) {
                D(i, j) = this->derivative(probeNormal, source, probe);
            }
        }
    }
}

template<typename T>
double Vacuum<T>::compDiagonalElementS(double area) const
{
    return (factor * sqrt(4 * M_PI / area));
}

template<typename T>
double Vacuum<T>::compDiagonalElementD(double area, double radius) const
{
    return (- factor * sqrt(M_PI / area) / radius);
}

template <typename T>
std::ostream & Vacuum<T>::printObject(std::ostream & os)
{
    os << "Vacuum" << std::endl;
    os << "Delta = " << this->delta_ << std::endl;
    os << "Uniform = " << this->uniform_ << std::endl;
    os << "epsilon = " << epsilon_;
    return os;
}

template class Vacuum<double>;
template class Vacuum<AD_directional>;
template class Vacuum<AD_gradient>;
template class Vacuum<AD_hessian>;
