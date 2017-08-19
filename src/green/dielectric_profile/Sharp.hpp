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

#pragma once

#include <iosfwd>

#include "Config.hpp"

/*! \file Sharp.hpp
 *  \struct Sharp
 *  \brief A sharp dielectric separation
 *  \author Roberto Di Remigio
 *  \date 2015
 */

namespace pcm {
namespace dielectric_profile {
struct Sharp __final {
  double epsilon;
  double epsilonSolvent;
  double radius;
  Sharp() : epsilon(1.0), epsilonSolvent(1.0), radius(1.0) {}
  Sharp(double eL, double eR, double c)
      : epsilon(eL), epsilonSolvent(eR), radius(c) {}
  friend std::ostream & operator<<(std::ostream & os, Sharp & obj) {
    os << "Sphere permittivity  = " << obj.epsilon << std::endl;
    os << "Solvent permittivity = " << obj.epsilonSolvent << std::endl;
    os << "Sphere radius        = " << obj.radius << " AU";
    return os;
  }
};
} // namespace dielectric_profile
} // namespace pcm
