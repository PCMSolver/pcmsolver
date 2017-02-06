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

#ifndef CAVITY_HPP
#define CAVITY_HPP

#include "Config.hpp"

#include "ICavity.hpp"
#include "CavityData.hpp"
#include "GePolCavity.hpp"
#include "RestartCavity.hpp"
#include "utils/Factory.hpp"

/*!
 * \file Cavity.hpp
 * \brief Top-level include file for cavities
 * \author Roberto Di Remigio
 * \date 2016
 *
 * Includes all cavity-related headers and defines the bootstrap function
 * for the Factory<ICavity, CavityData>
 */

namespace pcm {
namespace cavity {
inline ICavity * createGePolCavity(const CavityData & data) {
  return new GePolCavity(
      data.molecule, data.area, data.probeRadius, data.minimalRadius);
}

inline ICavity * createRestartCavity(const CavityData & data) {
  return new RestartCavity(data.filename);
}

inline void bootstrapFactory() {
  const bool registeredGePol =
      Factory<ICavity, CavityData>::TheFactory().registerObject("GEPOL",
                                                                createGePolCavity);
  if (!registeredGePol)
    PCMSOLVER_ERROR("Subscription of GePol cavity to factory failed!");

  const bool registeredRestart =
      Factory<ICavity, CavityData>::TheFactory().registerObject("RESTART",
                                                                createRestartCavity);
  if (!registeredRestart)
    PCMSOLVER_ERROR("Subscription of restart cavity to factory failed!");
}
} // namespace cavity
} // namespace pcm

#endif // CAVITY_HPP
