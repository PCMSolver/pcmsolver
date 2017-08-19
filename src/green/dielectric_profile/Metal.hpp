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

#include <complex>
#include <iosfwd>

#include "Config.hpp"

/*! \file Metal.hpp
 *  \struct Metal
 *  \brief An object with complex permittivity
 *  \author Roberto Di Remigio
 *  \date 2015
 */

namespace pcm {
namespace dielectric_profile {
struct Metal __final {
  std::complex<double> epsilon;
  Metal() : epsilon(std::complex<double>(1.0, 1.0)) {}
  Metal(double eRe, double eIm) : epsilon(std::complex<double>(eRe, eIm)) {}
  Metal(const std::complex<double> & e) : epsilon(e) {}
  friend std::ostream & operator<<(std::ostream & os, Metal & arg) {
    os << "Permittivity = " << arg.epsilon;
    return os;
  }
};
} // namespace dielectric_profile
} // namespace pcm
