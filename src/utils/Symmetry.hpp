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

#include <algorithm>
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
  int nrGenerators_;
  /*!
   * Generators
   */
  int generators_[3];
  /*!
   * Number of irreps
   */
  int nrIrrep_;

public:
  /*! \brief Default constructor sets up C1 point group
   */
  Symmetry() : nrGenerators_(0) {
    std::fill(generators_, generators_ + 3, 0);
    nrIrrep_ = int(std::pow(2.0, nrGenerators_));
  }
  Symmetry(int nr_gen, int gen[3]) : nrGenerators_(nr_gen) {
    // Transfer the passed generators array into generators_
    std::copy(gen, gen + nrGenerators_, generators_);
    // We can now initialize the number of irreps
    nrIrrep_ = int(std::pow(2.0, nrGenerators_));
  }
  Symmetry(const Symmetry & other)
      : nrGenerators_(other.nrGenerators_), nrIrrep_(other.nrIrrep_) {
    std::copy(other.generators_, other.generators_ + nrGenerators_, generators_);
  }
  ~Symmetry() {}
  int nrGenerators() const { return nrGenerators_; }
  int generators(int i) const { return generators_[i]; }
  int nrIrrep() const { return nrIrrep_; }
};

/*! Builds Symmetry object.
 *
 * \note C1 is built as Symmetry C1 = buildGroup(0, 0, 0, 0);
 */
Symmetry buildGroup(int _nr_gen, int _gen1, int _gen2, int _gen3);
