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
 *     PCMSolver API, see: <https://repo.ctcc.no/projects/pcmsolver>
 */
/* pcmsolver_copyright_end */

#include "Vacuum.hpp"

#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"
#include "TaylorPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "DiagonalIntegrator.hpp"
#include "GreensFunction.hpp"
#include "IGreensFunction.hpp"

template<typename T>
double Vacuum<T>::derivative(const Eigen::Vector3d & direction,
                             const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    return this->derivativeProbe(direction, p1, p2);
    //    return direction.dot(g);  // NORMALIZTION TEMPORARY REMOVED /direction.norm();
}

template<typename T>
T Vacuum<T>::operator()(T * sp, T * pp) const
{
    T res;
    res = 1.0/sqrt((sp[0]-pp[0])*(sp[0]-pp[0])+
                   (sp[1]-pp[1])*(sp[1]-pp[1])+
                   (sp[2]-pp[2])*(sp[2]-pp[2]));
    return res;
}

template <typename T>
double Vacuum<T>::diagonalS(double area) const {
        return this->diagonal_->computeS(this, area);
}

template <typename T>
double Vacuum<T>::diagonalD(double area, double radius) const {
        return this->diagonal_->computeD(this, area, radius);
}

/*
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
*/

template <typename T>
std::ostream & Vacuum<T>::printObject(std::ostream & os)
{
    os << "Green's function type: vacuum";
    return os;
}

template class Vacuum<double>;
template class Vacuum<AD_directional>;
template class Vacuum<AD_gradient>;
template class Vacuum<AD_hessian>;
