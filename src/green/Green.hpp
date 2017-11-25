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

#include "AnisotropicLiquid.hpp"
#include "GreenData.hpp"
#include "IGreensFunction.hpp"
#include "IonicLiquid.hpp"
#include "SphericalDiffuse.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"
#include "utils/Factory.hpp"

/*!
 * \file Solver.hpp
 * \brief Top-level include file for solvers
 * \author Roberto Di Remigio
 * \date 2016
 *
 * Includes all solver-related headers and defines the bootstrap function
 * for the Factory<ISolver, SolverData>
 */

namespace pcm {
namespace green {
namespace detail {
typedef pcm::function<IGreensFunction *(const GreenData &)> CreateGreensFunction;
} // namespace detail

inline Factory<detail::CreateGreensFunction> bootstrapFactory() {
  Factory<detail::CreateGreensFunction> factory_;

  factory_.subscribe("VACUUM", createVacuum);
  factory_.subscribe("UNIFORMDIELECTRIC", createUniformDielectric);
  factory_.subscribe("SPHERICALDIFFUSE", createSphericalDiffuse);
  factory_.subscribe("IONICLIQUID", createIonicLiquid);
  factory_.subscribe("ANISOTROPICLIQUID", createAnisotropicLiquid);

  return factory_;
}
} // namespace green
} // namespace pcm
