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

#include "AnisotropicLiquid.hpp"

#include <cmath>
#include <ostream>
#include <stdexcept>

#include "Config.hpp"

#include "EigenPimpl.hpp"
#include "TaylorPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "DiagonalIntegrator.hpp"
#include "Element.hpp"
#include "GreensFunction.hpp"
#include "IGreensFunction.hpp"
#include "MathUtils.hpp"

template <typename T>
void AnisotropicLiquid<T>::build() {
	// Initialize some internals: molecule-fixed to lab-fixed frame rotation
	// matrix, permittivity tensor in molecule-fixed frame and its inverse
	// 1. construct rotation matrix from Euler angles
	eulerRotation(R_, eulerAngles_);
	// 2. Apply the rotation matrix: epsilon_ = R_^t * epsilonLab_ * R_
	epsilon_ = R_.transpose() * epsilonLab_.asDiagonal() * R_;
	// 3. Obtain epsilonInv_ = R_ * epsilonLab_^-1 * R_^t
	Eigen::Vector3d scratch; 
	scratch << (1.0/epsilonLab_(0)), (1.0/epsilonLab_(1)), (1.0/epsilonLab_(2));
	epsilonInv_ = R_ * scratch.asDiagonal() * R_.transpose();
}

template<typename T>
double AnisotropicLiquid<T>::derivative(const Eigen::Vector3d & direction,
                                  const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    // NORMALIZATION TEMPORARILY REMOVED /direction.norm();
    return 0.0; //epsilon_ * (this->derivativeProbe(direction, p1, p2));
}

template<typename T>
T AnisotropicLiquid<T>::operator()(T * sp, T * pp) const
{
    T distance = sqrt((sp[0] - pp[0]) * (sp[0] - pp[0]) +
                      (sp[1] - pp[1]) * (sp[1] - pp[1]) +
                      (sp[2] - pp[2]) * (sp[2] - pp[2]));
    return 0.0;//(exp(-kappa_ * distance) / (epsilon_ * distance));
}

template <typename T>
double AnisotropicLiquid<T>::diagonalS(const Element & e) const {
        return 0.0; //this->diagonal_->computeS(this, e);
}

template <typename T>
double AnisotropicLiquid<T>::diagonalD(const Element & e) const {
        return 0.0; //this->diagonal_->computeD(this, e);
}

template <typename T>
std::ostream & AnisotropicLiquid<T>::printObject(std::ostream & os)
{
    os << "Green's function type: anisotropic liquid" << std::endl;
    os << "Permittivity tensor diagonal (lab frame) = " << epsilon_.transpose() << std::endl;
    os << "Euler angles (molecule-to-lab frame)     = " << eulerAngles_.transpose();
    return os;
}

template class AnisotropicLiquid<double>;
template class AnisotropicLiquid<AD_directional>;
template class AnisotropicLiquid<AD_gradient>;
template class AnisotropicLiquid<AD_hessian>;
