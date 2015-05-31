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

#ifndef INTEGRATORTYPES_HPP
#define INTEGRATORTYPES_HPP

#include "Config.hpp"

#include <boost/mpl/vector.hpp>

#include "DerivativeTypes.hpp"
#include "CollocationIntegrator.hpp"
#include "NumericalIntegrator.hpp"

typedef boost::mpl::vector<
    CollocationIntegrator<Numerical, Uniform>,
    CollocationIntegrator<AD_directional, Uniform>,
    CollocationIntegrator<AD_gradient, Uniform>,
    CollocationIntegrator<AD_hessian, Uniform>,
    NumericalIntegrator<Numerical, Uniform>,
    NumericalIntegrator<AD_directional, Uniform>,
    NumericalIntegrator<AD_gradient, Uniform>,
    NumericalIntegrator<AD_hessian, Uniform>
> integrator_types_uniform;

typedef boost::mpl::vector<
    CollocationIntegrator<Numerical, Yukawa>,
    CollocationIntegrator<AD_directional, Yukawa>,
    CollocationIntegrator<AD_gradient, Yukawa>,
    CollocationIntegrator<AD_hessian, Yukawa>,
    NumericalIntegrator<Numerical, Yukawa>,
    NumericalIntegrator<AD_directional, Yukawa>,
    NumericalIntegrator<AD_gradient, Yukawa>,
    NumericalIntegrator<AD_hessian, Yukawa>
> integrator_types_yukawa;

typedef boost::mpl::vector<
    CollocationIntegrator<Numerical, Anisotropic>,
    CollocationIntegrator<AD_directional, Anisotropic>,
    CollocationIntegrator<AD_gradient, Anisotropic>,
    CollocationIntegrator<AD_hessian, Anisotropic>,
    NumericalIntegrator<Numerical, Anisotropic>,
    NumericalIntegrator<AD_directional, Anisotropic>,
    NumericalIntegrator<AD_gradient, Anisotropic>,
    NumericalIntegrator<AD_hessian, Anisotropic>
> integrator_types_anisotropic;

typedef boost::mpl::vector<
    CollocationIntegrator<Numerical, TanhDiffuse>,
    CollocationIntegrator<AD_directional, TanhDiffuse>,
    CollocationIntegrator<AD_gradient, TanhDiffuse>,
    CollocationIntegrator<AD_hessian, TanhDiffuse>,
    NumericalIntegrator<Numerical, TanhDiffuse>,
    NumericalIntegrator<AD_directional, TanhDiffuse>,
    NumericalIntegrator<AD_gradient, TanhDiffuse>,
    NumericalIntegrator<AD_hessian, TanhDiffuse>
> integrator_types_tanhdiffuse;

#endif // INTEGRATORTYPES_HPP
