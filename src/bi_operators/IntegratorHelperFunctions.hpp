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

#ifndef INTEGRATORHELPERFUNCTIONS_HPP
#define INTEGRATORHELPERFUNCTIONS_HPP

#include <cmath>
#include <functional>
#include <vector>

#include "Config.hpp"

#include <Eigen/Dense>

#include "Element.hpp"

namespace integrator {
/*! \typedef DiagonalS
 *  \brief functor handle to the calculation of the diagonal of S
 */
typedef std::function<double(double)> DiagonalS;

/*! \typedef KernelS
 *  \brief functor handle to the kernelS method in IGreensFunction
 */
typedef std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> KernelS;

/*! \typedef DiagonalD
 *  \brief functor handle to the calculation of the diagonal of D
 */
typedef std::function<double(double, double)> DiagonalD;

/*! \typedef KernelD
 *  \brief functor handle to the kernelD method in IGreensFunction
 */
typedef std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &, const Eigen::Vector3d &)> KernelD;

/*! Returns matrix representation of the single layer operator by collocation
 *  \param[in] elements list of finite elements
 *  \param[in] diag     functor for the evaluation of the diagonal of S
 *  \param[in] kern     function for the evaluation of the off-diagonal of S
 */
inline Eigen::MatrixXd singleLayer(const std::vector<Element> & elements,
                                   const DiagonalS & diag, const KernelS & kern)
{
    size_t mat_size = elements.size();
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(mat_size, mat_size);
    for (size_t i = 0; i < mat_size; ++i) {
        // Fill diagonal
        S(i, i) = diag(elements[i].area());
        Eigen::Vector3d source = elements[i].center();
        for (size_t j = 0; j < mat_size; ++j) {
            // Fill off-diagonal
            Eigen::Vector3d probe = elements[j].center();
            if (i != j) S(i, j) = kern(source, probe);
        }
    }
    return S;
}

/*! Returns matrix representation of the double layer operator by collocation
 *  \param[in] elements list of finite elements
 *  \param[in] diag     functor for the evaluation of the diagonal of D
 *  \param[in] kern     function for the evaluation of the off-diagonal of D
 */
inline Eigen::MatrixXd doubleLayer(const std::vector<Element> & elements,
                                   const DiagonalD & diag, const KernelD & kern)
{
    size_t mat_size = elements.size();
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(mat_size, mat_size);
    for (size_t i = 0; i < mat_size; ++i) {
        // Fill diagonal
        D(i, i) = diag(elements[i].area(), elements[i].sphere().radius());
        Eigen::Vector3d source = elements[i].center();
        for (size_t j = 0; j < mat_size; ++j) {
            // Fill off-diagonal
            Eigen::Vector3d probe = elements[j].center();
            Eigen::Vector3d probeNormal = elements[j].normal();
            probeNormal.normalize();
            if (i != j) D(i, j) = kern(probeNormal, source, probe);
        }
    }
    return D;
}
} // namespace integrator

#endif // INTEGRATORHELPERFUNCTIONS_HPP
