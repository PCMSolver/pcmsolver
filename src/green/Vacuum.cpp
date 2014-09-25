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

#include "Vacuum.hpp"

#include <cmath>

#include "Config.hpp"

#include <Eigen/Dense>
#include "taylor.hpp"

#include "DerivativeTypes.hpp"
#include "DiagonalIntegrator.hpp"
#include "Element.hpp"
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
double Vacuum<T>::diagonalS(const Element & e) const {
        return this->diagonal_->computeS(this, e);
}

template <typename T>
double Vacuum<T>::diagonalD(const Element & e) const {
        return this->diagonal_->computeD(this, e);
}

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
