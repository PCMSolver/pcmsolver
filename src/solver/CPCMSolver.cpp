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

#include "CPCMSolver.hpp"

#include <fstream>
#include <iostream>
#include <string>

#include "Config.hpp"

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include "SolverData.hpp"
#include "bi_operators/IBoundaryIntegralOperator.hpp"
#include "cavity/Element.hpp"
#include "cavity/ICavity.hpp"
#include "green/IGreensFunction.hpp"
#include "utils/MathUtils.hpp"

namespace pcm {
namespace solver {
void CPCMSolver::buildSystemMatrix_impl(const ICavity & cavity,
                                        const IGreensFunction & gf_i,
                                        const IGreensFunction & gf_o,
                                        const IBoundaryIntegralOperator & op) {
  if (!isotropic_)
    PCMSOLVER_ERROR("C-PCM is defined only for isotropic environments!");
  TIMER_ON("Computing S");
  double epsilon = dielectric_profile::epsilon(gf_o.permittivity());
  S_ = op.computeS(cavity, gf_i);
  S_ /= (epsilon - 1.0) / (epsilon + correction_);
  // Get in Hermitian form
  if (hermitivitize_)
    utils::hermitivitize(S_);
  TIMER_OFF("Computing S");

  // Symmetry-pack
  // The number of irreps in the group
  int nrBlocks = cavity.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cavity.irreducible_size();
  utils::symmetryPacking(blockS_, S_, dimBlock, nrBlocks);

  built_ = true;
}

Eigen::VectorXd CPCMSolver::computeCharge_impl(const Eigen::VectorXd & potential,
                                               int irrep) const {
  // The potential and charge vector are of dimension equal to the
  // full dimension of the cavity. We have to select just the part
  // relative to the irrep needed.
  int fullDim = S_.rows();
  Eigen::VectorXd charge = Eigen::VectorXd::Zero(fullDim);
  int nrBlocks = blockS_.size();
  int irrDim = fullDim / nrBlocks;
  charge.segment(irrep * irrDim, irrDim) =
      -blockS_[irrep].ldlt().solve(potential.segment(irrep * irrDim, irrDim));

  return charge;
}

std::ostream & CPCMSolver::printSolver(std::ostream & os) {
  os << "Solver Type: C-PCM" << std::endl;
  if (hermitivitize_) {
    os << "PCM matrix hermitivitized" << std::endl;
  } else {
    os << "PCM matrix NOT hermitivitized (matches old DALTON)" << std::endl;
  }
  os << "Correction = " << correction_;

  return os;
}

ISolver * createCPCMSolver(const SolverData & data) {
  return new CPCMSolver(data.hermitivitize, data.correction);
}
} // namespace solver
} // namespace pcm
