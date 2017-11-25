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

#include "CPCMSolver.hpp"
#include "IEFSolver.hpp"
#include "ISolver.hpp"
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
namespace solver {
namespace detail {
typedef pcm::function<ISolver *(const SolverData &)> CreateSolver;
} // namespace detail

inline Factory<detail::CreateSolver> bootstrapFactory() {
  Factory<detail::CreateSolver> factory_;

  factory_.subscribe("CPCM", createCPCMSolver);
  factory_.subscribe("IEFPCM", createIEFSolver);

  return factory_;
}
} // namespace solver
} // namespace pcm
