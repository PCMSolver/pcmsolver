/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
 *
 *     This file is part of PCMSolver.
 *
 *     PCMSolver is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#include "UniformDielectric.hpp"

#include <cmath>
#include <ostream>

#include "Config.hpp"

#include <Eigen/Dense>
#include "taylor.hpp"

#include "DerivativeTypes.hpp"
#include "DiagonalIntegrator.hpp"
#include "GreensFunction.hpp"
#include "IGreensFunction.hpp"

template<typename T>
double UniformDielectric<T>::derivative(const Eigen::Vector3d & direction,
                                        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    return epsilon_ * (this->derivativeProbe(direction, p1,
                       p2));  // NORMALIZTION TEMPORARY REMOVED /direction.norm();
}

template<typename T>
T UniformDielectric<T>::operator()(T sp[3], T pp[3]) const
{
    T distance = sqrt((sp[0] - pp[0]) * (sp[0] - pp[0]) +
                      (sp[1] - pp[1]) * (sp[1] - pp[1]) +
                      (sp[2] - pp[2]) * (sp[2] - pp[2]));
    return 1/(epsilon_ * distance);
}

template <typename T>
double UniformDielectric<T>::diagonalS(double area) const {
        return this->diagonal_->computeS(this, area);
}

template <typename T>
double UniformDielectric<T>::diagonalD(double area, double radius) const {
        return this->diagonal_->computeD(this, area, radius);
}
/*
template <typename T>
void UniformDielectric<T>::operator()(Eigen::MatrixXd & S, Eigen::MatrixXd & D,
                                      const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                                      const Eigen::VectorXd & areas, const Eigen::VectorXd & radii) const
{
    int size = S.rows();

    for (int i = 0; i < size; ++i) {
        S(i, i) = factor * std::sqrt(4 * M_PI / areas(i)) * (1.0 / epsilon_);
        D(i, i) = -factor * std::sqrt(M_PI/ areas(i)) * (1.0 / radii(
                      i)); // Check this! There might be an epsilon_?
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
void UniformDielectric<T>::operator()(Eigen::MatrixXd & S,
                                      const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                                      const Eigen::VectorXd & areas) const
{
    int size = S.rows();

    for (int i = 0; i < size; ++i) {
        S(i, i) = factor * std::sqrt(4 * M_PI / areas(i)) * (1.0 / epsilon_);
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
void UniformDielectric<T>::operator()(Eigen::MatrixXd & D,
                                      const Eigen::MatrixXd & centers, const Eigen::MatrixXd & normals,
                                      const Eigen::VectorXd & areas, const Eigen::VectorXd & radii) const
{
    int size = D.rows();

    for (int i = 0; i < size; ++i) {
        D(i, i) = -factor * std::sqrt(M_PI/ areas(i)) * (1.0 / radii(
                      i)); // Check this! There might be an epsilon_?
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
double UniformDielectric<T>::compDiagonalElementS(double area) const
{
    return factor * sqrt(4 * M_PI / area) / epsilon_;
}

template<typename T>
double UniformDielectric<T>::compDiagonalElementD(double area, double radius) const
{
    //	double s = factor * sqrt(4 * M_PI / area);
    //	return s / (2 * radius);
    return - factor * sqrt(M_PI / area) / radius;
}
*/

template <typename T>
std::ostream & UniformDielectric<T>::printObject(std::ostream & os)
{
    os << "Green's function type: uniform dielectric" << std::endl;
    os << "Permittivity = " << epsilon_ << std::endl;
    return os;
}

template class UniformDielectric<double>;
template class UniformDielectric<AD_directional>;
template class UniformDielectric<AD_gradient>;
template class UniformDielectric<AD_hessian>;
