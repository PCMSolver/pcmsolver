/**
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2016 Roberto Di Remigio, Luca Frediani and collaborators.
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

#include "BoundaryIntegralOperator.hpp"

#include "Config.hpp"

#include <Eigen/Core>

#include "cavity/Cavity.hpp"
#include "green/IGreensFunction.hpp"
#include "utils/MathUtils.hpp"

Eigen::MatrixXd BoundaryIntegralOperator::computeS(
    const Cavity & cav, const IGreensFunction & gf) const {
  Eigen::MatrixXd biop = computeS_impl(cav.elements(), gf);
  // Perform symmetry blocking
  // The total size of the cavity
  PCMSolverIndex cavitySize = cav.size();
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();
  if (cav.pointGroup().nrGenerators() != 0) {
    TIMER_ON("Symmetry blocking");
    symmetryBlocking(biop, cavitySize, dimBlock, nrBlocks);
    TIMER_OFF("Symmetry blocking");
  }
  return biop;
}

Eigen::MatrixXd BoundaryIntegralOperator::computeD(
    const Cavity & cav, const IGreensFunction & gf) const {
  Eigen::MatrixXd biop = computeD_impl(cav.elements(), gf);
  // Perform symmetry blocking
  // The total size of the cavity
  PCMSolverIndex cavitySize = cav.size();
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();
  if (cav.pointGroup().nrGenerators() != 0) {
    TIMER_ON("Symmetry blocking");
    symmetryBlocking(biop, cavitySize, dimBlock, nrBlocks);
    TIMER_OFF("Symmetry blocking");
  }
  return biop;
}
