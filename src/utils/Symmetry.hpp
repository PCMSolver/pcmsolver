/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2018 Roberto Di Remigio, Luca Frediani and contributors.
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
#include <bitset>
#include <cmath>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

/*! \file Symmetry.hpp */

namespace pcm {
namespace utils {
/*! \class Symmetry
 *  \brief Contains very basic info about symmetry (only Abelian groups)
 *  \author Roberto Di Remigio
 *  \date 2014, 2018
 *  \note C1 is built as Symmetry C1(0, 0, 0, 0);
 *
 *  Indexing of symmetry operations and their mapping to a bitstring:
 *       zyx         Parity
 *    0  000    E      1.0
 *    1  001   Oyz    -1.0
 *    2  010   Oxz    -1.0
 *    3  011   C2z     1.0
 *    4  100   Oxy    -1.0
 *    5  101   C2y     1.0
 *    6  110   C2x     1.0
 *    7  111    i     -1.0
 */
class Symmetry __final {
public:
#ifdef HAS_CXX11
  Symmetry()
      : nrGenerators_(0),
        nrIrrep_(static_cast<int>(std::pow(2.0, nrGenerators_))),
        generators_({{0, 0, 0}}) {}
  Symmetry(int nr_gen, int gen1, int gen2, int gen3)
      : nrGenerators_(nr_gen),
        nrIrrep_(static_cast<int>(std::pow(2.0, nrGenerators_))),
        generators_({{gen1, gen2, gen3}}) {}
#else
  Symmetry()
      : nrGenerators_(0), nrIrrep_(static_cast<int>(std::pow(2.0, nrGenerators_))) {
    std::fill(generators_, generators_ + 3, 0);
  }
  Symmetry(int nr_gen, int gen1, int gen2, int gen3)
      : nrGenerators_(nr_gen),
        nrIrrep_(static_cast<int>(std::pow(2.0, nrGenerators_))),
  {
    generators_[0] = gen1;
    generators_[1] = gen2;
    generators_[2] = gen3;
  }
#endif
  Symmetry(int nr_gen, int gen[3])
      : nrGenerators_(nr_gen),
        nrIrrep_(static_cast<int>(std::pow(2.0, nrGenerators_))) {
    generators_[0] = gen[0];
    generators_[1] = gen[1];
    generators_[2] = gen[2];
  }

  int nrGenerators() const { return nrGenerators_; }
  int nrIrrep() const { return nrIrrep_; }
  int generators(size_t i) const { return generators_[i]; }

private:
  /*! Number of generators */
  int nrGenerators_;
  /*! Number of irreps */
  int nrIrrep_;
  /*! Generators */
  pcm::array<int, 3> generators_;
};

/*! \brief Transform matrix to block diagonal form by symmetry
 *  \param[out] matrix the matrix to be block-diagonalized
 *  \param[in]  cavitySize the size of the cavity (size of the matrix)
 *  \param[in]  ntsirr     the size of the irreducible portion of the cavity (size of
 * the blocks)
 *  \param[in]  nr_irrep   the number of irreducible representations (number of
 * blocks)
 *
 * This function implements the symmetry-blocking of the PCM matrix due to
 * point group symmetry as reported in:
 * L. Frediani, R. Cammi, C. S. Pomelli, J. Tomasi and K. Ruud, J. Comput.Chem. 25,
 * 375 (2003)
 */
void symmetryBlocking(Eigen::MatrixXd & matrix,
                      PCMSolverIndex cavitySize,
                      PCMSolverIndex ntsirr,
                      int nr_irrep);

/*! \brief Packs symmetry blocked matrix, i.e. only stores non-zero blocks
 *  \param[out] blockedMatrix the result of packing fullMatrix
 *  \param[in]  fullMatrix the matrix to be packed
 *  \param[in]  dimBlock the dimension of the square blocks
 *  \param[in]  nrBlocks the number of square blocks
 *
 *  This function packs the square block diagonal fullMatrix with nrBlocks of
 *  dimension dimBlock into a std::vector<Eigen::MatrixXd> containing nrBlocks
 *  square matrices of dimenion dimBlock.
 */
void symmetryPacking(std::vector<Eigen::MatrixXd> & blockedMatrix,
                     const Eigen::MatrixXd & fullMatrix,
                     int dimBlock,
                     int nrBlocks);

/*! \brief Calculate the parity of the bitset
 *  \param[in] bitrep a bitset
 *  \tparam nBits length of the input bitset
 *
 *  Calculate the parity of the bitset as defined by:
 *     bitrep[0] XOR bitrep[1] XOR ... XOR bitrep[nBits-1]
 */
template <size_t nBits> inline int parity(std::bitset<nBits> bitrep) {
  int parity = 0;
  for (size_t i = 0; i < bitrep.size(); ++i) {
    parity ^= bitrep[i];
  }
  return parity;
}

/*! \brief Returns parity of input integer.
 *  \param[in] i an integer, usually an index for an irrep or a symmetry operation
 *
 * The parity is defined as the result of using XOR on the bitrep
 * of the given integer. For example:
 *   2 -> 010 -> 0^1^0 = 1 -> -1.0
 *   6 -> 110 -> 1^1^0 = 0 ->  1.0
 *
 * It can also be interpreted as the action of a given operation on the
 * Cartesian axes:
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
inline double parity(unsigned int i) {
  return (parity(std::bitset<3>(i)) ? -1.0 : 1.0);
}
} // namespace utils
} // namespace pcm
