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

#include "IonicLiquid.hpp"

#include <cmath>
#include <ostream>
#include <stdexcept>

#include "Config.hpp"

#include <Eigen/Dense>
#include "taylor.hpp"

#include "DerivativeTypes.hpp"
#include "DiagonalIntegrator.hpp"
#include "GreensFunction.hpp"
#include "IGreensFunction.hpp"

template<typename T>
double IonicLiquid<T>::derivative(const Eigen::Vector3d & direction,
                                  const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
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
double IonicLiquid<T>::diagonalS(double area) const {
        this->diagonal_->computeS(this, area);
        return 1.0;
}

template <typename T>
double IonicLiquid<T>::diagonalD(double area, double radius) const {
        this->diagonal_->computeD(this, area, radius);
        return 1.0;
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
