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

template <typename DerivativeTraits>
inline DerivativeTraits distance(DerivativeTraits u[3], DerivativeTraits v[3]) {
  return sqrt(pow(u[0] - v[0], 2) + pow(u[1] - v[1], 2) + pow(u[2] - v[2], 2));
}

template <typename DerivativeTraits>
inline DerivativeTraits norm(DerivativeTraits u[3]) {
  return sqrt(pow(u[0], 2) + pow(u[1], 2) + pow(u[2], 2));
}

template <typename DerivativeTraits>
inline DerivativeTraits dot_product(DerivativeTraits u[3], DerivativeTraits v[3]) {
  return (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]);
}
