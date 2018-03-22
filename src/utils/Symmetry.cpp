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

#include "Symmetry.hpp"

#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include "MathUtils.hpp"

namespace pcm {
namespace utils {
void symmetryBlocking(Eigen::MatrixXd & matrix,
                      PCMSolverIndex cavitySize,
                      PCMSolverIndex ntsirr,
                      int nr_irrep) {
  // u is the character table for the group (t in the paper)
  Eigen::MatrixXd u = Eigen::MatrixXd::Zero(nr_irrep, nr_irrep);
  for (int i = 0; i < nr_irrep; ++i) {
    for (int j = 0; j < nr_irrep; ++j) {
      u(i, j) = parity(i & j);
    }
  }
  // Naming of indices:
  //     a, b, c, d   run over the total size of the cavity (N)
  //     i, j, k, l   run over the number of irreps (n)
  //     p, q, r, s   run over the irreducible size of the cavity (N/n)
  // Instead of forming U (T in the paper) and then perform the dense
  // multiplication, we multiply block-by-block using just the u matrix.
  //      matrix = U * matrix * Ut; U * Ut = Ut * U = id
  // First half-transformation, i.e. first_half = matrix * Ut
  Eigen::MatrixXd first_half = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
  for (int i = 0; i < nr_irrep; ++i) {
    int ioff = i * ntsirr;
    for (int k = 0; k < nr_irrep; ++k) {
      int koff = k * ntsirr;
      for (int j = 0; j < nr_irrep; ++j) {
        int joff = j * ntsirr;
        double ujk = u(j, k) / nr_irrep;
        for (int p = 0; p < ntsirr; ++p) {
          int a = ioff + p;
          for (int q = 0; q < ntsirr; ++q) {
            int b = joff + q;
            int c = koff + q;
            first_half(a, c) += matrix(a, b) * ujk;
          }
        }
      }
    }
  }
  // Second half-transformation, i.e. matrix = U * first_half
  matrix.setZero(cavitySize, cavitySize);
  for (int i = 0; i < nr_irrep; ++i) {
    int ioff = i * ntsirr;
    for (int k = 0; k < nr_irrep; ++k) {
      int koff = k * ntsirr;
      for (int j = 0; j < nr_irrep; ++j) {
        int joff = j * ntsirr;
        double uij = u(i, j);
        for (int p = 0; p < ntsirr; ++p) {
          int a = ioff + p;
          int b = joff + p;
          for (int q = 0; q < ntsirr; ++q) {
            int c = koff + q;
            matrix(a, c) += uij * first_half(b, c);
          }
        }
      }
    }
  }
  // Traverse the matrix and discard numerical zeros
  for (PCMSolverIndex a = 0; a < cavitySize; ++a) {
    for (PCMSolverIndex b = 0; b < cavitySize; ++b) {
      if (numericalZero(matrix(a, b))) {
        matrix(a, b) = 0.0;
      }
    }
  }
}

void symmetryPacking(std::vector<Eigen::MatrixXd> & blockedMatrix,
                     const Eigen::MatrixXd & fullMatrix,
                     int dimBlock,
                     int nrBlocks) {
  int j = 0;
  for (int i = 0; i < nrBlocks; ++i) {
    blockedMatrix.push_back(fullMatrix.block(j, j, dimBlock, dimBlock));
    j += dimBlock;
  }
}
} // namespace utils
} // namespace pcm
