/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2017 Roberto Di Remigio, Luca Frediani and collaborators.
 *
 * This file is part of PCMSolver.
 *
 * PCMSolver is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PCMSolver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to the
 * PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */

#pragma once

#include "Config.hpp"

#include "taylor.hpp"

#include <boost/mpl/vector.hpp>

namespace pcm {
typedef double Stencil;
typedef taylor<double, 1, 1> AD_directional;
typedef taylor<double, 3, 1> AD_gradient;
typedef taylor<double, 3, 2> AD_hessian;

namespace green {
typedef boost::mpl::vector<Stencil, AD_directional, AD_gradient, AD_hessian>
    derivative_types;
} // namespace green
} // namespace pcm
