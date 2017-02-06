/**
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

#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "Config.hpp"

#include "ISolver.hpp"
#include "SolverData.hpp"
#include "IEFSolver.hpp"
#include "CPCMSolver.hpp"
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
inline ISolver * createCPCMSolver(const SolverData & data) {
  return new CPCMSolver(data.hermitivitize, data.correction);
}

inline ISolver * createIEFSolver(const SolverData & data) {
  return new IEFSolver(data.hermitivitize);
}

inline void bootstrapFactory() {
  const bool registeredCPCMSolver =
      Factory<ISolver, SolverData>::TheFactory().registerObject("CPCM",
                                                                createCPCMSolver);
  if (!registeredCPCMSolver)
    PCMSOLVER_ERROR("Subscription of solver for CPCM to factory failed!");

  const bool registeredIEFSolver =
      Factory<ISolver, SolverData>::TheFactory().registerObject("IEFPCM",
                                                                createIEFSolver);
  if (!registeredIEFSolver)
    PCMSOLVER_ERROR("Subscription of solver for IEF to factory failed!");
}
} // namespace solver
} // namespace pcm

#endif // SOLVER_HPP
