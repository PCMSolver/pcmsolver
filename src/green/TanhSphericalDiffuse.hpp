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

#ifndef TANHSPHERICALDIFFUSE_HPP
#define TANHSPHERICALDIFFUSE_HPP

#include <cmath>
#include <iosfwd>
#include <string>

#include "Config.hpp"

#include <Eigen/Dense>
#include "taylor.hpp"

#include "DerivativeTypes.hpp"
#include "DiagonalIntegratorFactory.hpp"
#include "DiagonalIntegrator.hpp"
#include "ForIdGreen.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "GreensFunctionFactory.hpp"
#include "SphericalDiffuse.hpp"
#include "TanhDiffuse.hpp"

/*! \file TanhSphericalDiffuse.hpp
 *  \typedef TanhSphericalDiffuse
 *  \brief Green's function for a diffuse interface with spherical symmetry and tanh profile
 *  \author Hui Cao, Ville Weijo, Luca Frediani and Roberto Di Remigio
 *  \date 2010-2015
 */

typedef SphericalDiffuse<TanhDiffuse> TanhSphericalDiffuse;

template <>
inline Numerical GreensFunction<Numerical, TanhDiffuse>::function(const Eigen::Vector3d & source,
                                        const Eigen::Vector3d & probe) const
{
    // To be implemented
    return 0.0;
}

template <>
inline Numerical GreensFunction<Numerical, TanhDiffuse>::derivativeSource(const Eigen::Vector3d & normal_p1,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    // To be implemented
    return 0.0;
}

template <>
inline Numerical GreensFunction<Numerical, TanhDiffuse>::derivativeProbe(const Eigen::Vector3d & normal_p2,
        const Eigen::Vector3d & p1, const Eigen::Vector3d & p2) const
{
    // To be implemented
    return 0.0;
}

namespace
{
    IGreensFunction * createTanhSphericalDiffuse(const greenData & _data)
    {
        double eL, eR, w, c; // To be read from _data
        eL = 0.0, eR = 0.0, w = 0.0, c = 0.0;
        DiagonalIntegrator * integrator =
		DiagonalIntegratorFactory::TheDiagonalIntegratorFactory().createDiagonalIntegrator(_data.integratorType);
        return new TanhSphericalDiffuse(eL, eR, w, c, integrator);
    }
    const std::string TANHSPHERICALDIFFUSE("TANHSPHERICALDIFFUSE");
    const bool registeredTanhSphericalDiffuse =
        GreensFunctionFactory::TheGreensFunctionFactory().registerGreensFunction(
            TANHSPHERICALDIFFUSE, createTanhSphericalDiffuse);
}

#endif // TANHSPHERICALDIFFUSE_HPP
