/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2020 Roberto Di Remigio, Luca Frediani and contributors.
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

#include <array>
#include <cmath>

#include "Config.hpp"

/*! \file Symmetry.hpp
 *  \class Symmetry
 *  \brief Contains very basic info about symmetry (only Abelian groups)
 *  \author Roberto Di Remigio
 *  \date 2014
 *
 *  Just a wrapper around a vector containing the generators of the group
 */

class Symmetry {
private:
  /*!
   * Number of generators
   */
  int nrGenerators_{0};
  /*!
   * Generators
   */
  std::array<int, 3> generators_{0};
  /*!
   * Number of irreps
   */
  int nrIrrep_{1};

public:
  Symmetry() = default;
  Symmetry(int nr_gen, const std::array<int, 3> gens)
      : nrGenerators_(nr_gen), generators_(gens) {}
  Symmetry(int nr_gen, int g1, int g2, int g3)
      : nrGenerators_(nr_gen), generators_({g1, g2, g3}) {}
  int nrGenerators() const { return nrGenerators_; }
  int generators(int i) const { return generators_[i]; }
  int nrIrrep() const { return std::pow(2, nrGenerators_); }
};
