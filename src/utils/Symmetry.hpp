/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
 *
 *     This file is part of PCMSolver.
 *
 *     PCMSolver is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#ifndef SYMMETRY_HPP
#define SYMMETRY_HPP

#include "Config.hpp"

#include <Eigen/Dense>

/*! \file Symmetry.hpp
 *  \class Symmetry
 *  \brief Contains very basic info about symmetry (only Abelian groups)
 *  \author Roberto Di Remigio
 *  \date 2014
 *
 *  Just a wrapper around a vector containing the generators of the group
 */

class Symmetry
{
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
    Symmetry() {}
    Symmetry(int nr_gen, int gen[3]) : nrGenerators_(nr_gen) {
        // Transfer the passed generators array into generators_
        std::copy(gen, gen + nrGenerators_, generators_);
        // We can now initialize the number of irreps
        nrIrrep_ = int(std::pow(2.0, nrGenerators_));
    }
    Symmetry(const Symmetry & other) :
        nrGenerators_(other.nrGenerators_), nrIrrep_(other.nrIrrep_) {
        std::copy(other.generators_, other.generators_ + nrGenerators_, generators_);
    }
    ~Symmetry() {}
    int nrGenerators() const { return nrGenerators_; }
    int generators(int i) const { return generators_[i]; }
    int nrIrrep() const { return nrIrrep_; }
};

/*! Builds Symmetry object.
 */
Symmetry buildGroup(int _nr_gen, int _gen1, int _gen2, int _gen3);

/*! Returns parity of input integer.
 *      zyx         Parity
 *   0  000    E      1.0
 *   1  001   Oyz    -1.0
 *   2  010   Oxz    -1.0
 *   3  011   C2z     1.0
 *   4  100   Oxy    -1.0
 *   5  101   C2y     1.0
 *   6  110   C2x     1.0
 *   7  111    i     -1.0
 */
double parity(int i);

#endif // SYMMETRY_HPP
