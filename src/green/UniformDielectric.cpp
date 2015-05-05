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
#include "Element.hpp"
#include "GreensFunction.hpp"
#include "IGreensFunction.hpp"

template<typename T>
double UniformDielectric<T>::derivative(const Eigen::Vector3d & direction,
                                        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{ // NORMALIZTION TEMPORARY REMOVED /direction.norm();
    return epsilon_ * (this->derivativeProbe(direction, p1, p2)); 
}

template<typename T>
T UniformDielectric<T>::operator()(T sp[3], T pp[3]) const
{
    T distance = sqrt((sp[0] - pp[0]) * (sp[0] - pp[0]) +
                      (sp[1] - pp[1]) * (sp[1] - pp[1]) +
                      (sp[2] - pp[2]) * (sp[2] - pp[2]));
    return (1.0/(epsilon_ * distance));
}

template <typename T>
double UniformDielectric<T>::diagonalS(const Element & e) const {
        return this->diagonal_->computeS(this, e);
}

template <typename T>
double UniformDielectric<T>::diagonalD(const Element & e) const {
        return this->diagonal_->computeD(this, e);
}

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
