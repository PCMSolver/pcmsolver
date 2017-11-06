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

#include <iostream>

#include "Config.hpp"

/*! \file Yukawa.hpp
 *  \struct Yukawa
 *  \brief describes a medium with damping, i.e. ionic liquid
 *  \author Roberto Di Remigio
 *  \date 2015
 */

namespace pcm {
namespace dielectric_profile {
struct Yukawa __final {
  double epsilon;
  double kappa;
  Yukawa() : epsilon(1.0), kappa(0.0) {}
  Yukawa(double eps, double k) : epsilon(eps), kappa(k) {}
  friend std::ostream & operator<<(std::ostream & os, Yukawa & arg) {
    os << "Permittivity         = " << arg.epsilon << std::endl;
    os << "Inverse Debye length = " << arg.kappa;
    return os;
  }
};
} // namespace dielectric_profile
} // namespace pcm
