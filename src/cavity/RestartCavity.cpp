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

#include "RestartCavity.hpp"

#include <iostream>
#include <string>
#include <vector>

#include "Config.hpp"

#include "CavityData.hpp"

namespace pcm {
namespace cavity {
std::ostream & RestartCavity::printCavity(std::ostream & os) {
  os << "Cavity type: Restart" << std::endl;
  os << "Number of finite elements = " << nElements_ << std::endl;
  return os;
}

ICavity * createRestartCavity(const CavityData & data) {
  return new RestartCavity(data.filename);
}
} // namespace cavity
} // namespace pcm
