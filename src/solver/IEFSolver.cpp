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

#include "IEFSolver.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/LU>

#include "SolverData.hpp"
#include "SolverImpl.hpp"
#include "bi_operators/IBoundaryIntegralOperator.hpp"
#include "cavity/Element.hpp"
#include "cavity/ICavity.hpp"
#include "green/IGreensFunction.hpp"
#include "utils/MathUtils.hpp"

namespace pcm {
namespace solver {
void IEFSolver::buildSystemMatrix_impl(const ICavity & cavity,
                                       const IGreensFunction & gf_i,
                                       const IGreensFunction & gf_o,
                                       const IBoundaryIntegralOperator & op) {
  isotropic_ = (gf_i.uniform() && gf_o.uniform());
  isotropic_ ? buildIsotropicMatrix(cavity, gf_i, gf_o, op)
             : buildAnisotropicMatrix(cavity, gf_i, gf_o, op);
}

void IEFSolver::buildAnisotropicMatrix(const ICavity & cav,
                                       const IGreensFunction & gf_i,
                                       const IGreensFunction & gf_o,
                                       const IBoundaryIntegralOperator & op) {
  Tepsilon_ = detail::anisotropicTEpsilon(cav, gf_i, gf_o, op);
  Rinfinity_ = detail::anisotropicRinfinity(cav, gf_i, gf_o, op);

  // Pack into a block diagonal matrix
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();
  // For the moment just packs into a std::vector<Eigen::MatrixXd>
  utils::symmetryPacking(blockTepsilon_, Tepsilon_, dimBlock, nrBlocks);
  utils::symmetryPacking(blockRinfinity_, Rinfinity_, dimBlock, nrBlocks);

  built_ = true;
}

void IEFSolver::buildIsotropicMatrix(const ICavity & cav,
                                     const IGreensFunction & gf_i,
                                     const IGreensFunction & gf_o,
                                     const IBoundaryIntegralOperator & op) {
  Tepsilon_ = detail::isotropicTEpsilon(
      cav, gf_i, dielectric_profile::epsilon(gf_o.permittivity()), op);
  Rinfinity_ = detail::isotropicRinfinity(cav, gf_i, op);

  // Pack into a block diagonal matrix
  // The number of irreps in the group
  int nrBlocks = cav.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cav.irreducible_size();
  // For the moment just packs into a std::vector<Eigen::MatrixXd>
  utils::symmetryPacking(blockTepsilon_, Tepsilon_, dimBlock, nrBlocks);
  utils::symmetryPacking(blockRinfinity_, Rinfinity_, dimBlock, nrBlocks);

  built_ = true;
}

Eigen::VectorXd IEFSolver::computeCharge_impl(const Eigen::VectorXd & potential,
                                              int irrep) const {
  // The potential and charge vector are of dimension equal to the
  // full dimension of the cavity. We have to select just the part
  // relative to the irrep needed.
  int fullDim = Rinfinity_.rows();
  Eigen::VectorXd charge = Eigen::VectorXd::Zero(fullDim);
  int nrBlocks = blockRinfinity_.size();
  int irrDim = fullDim / nrBlocks;
  charge.segment(irrep * irrDim, irrDim) =
      -blockTepsilon_[irrep].partialPivLu().solve(
          blockRinfinity_[irrep] * potential.segment(irrep * irrDim, irrDim));

  // Obtain polarization weights
  if (hermitivitize_) {
    Eigen::VectorXd adj_asc = Eigen::VectorXd::Zero(fullDim);
    // Form T^\dagger * v = c
    adj_asc.segment(irrep * irrDim, irrDim) =
        blockTepsilon_[irrep].adjoint().partialPivLu().solve(
            potential.segment(irrep * irrDim, irrDim));
    // Form R^\dagger * c = q^* ("transposed" polarization charges)
    adj_asc.segment(irrep * irrDim, irrDim) =
        -blockRinfinity_[irrep].adjoint() *
        (adj_asc.segment(irrep * irrDim, irrDim).eval());
    // Get polarization weights
    charge = 0.5 * (adj_asc + charge.eval());
  }

  return charge;
}

std::ostream & IEFSolver::printSolver(std::ostream & os) {
  std::string type;
  if (isotropic_) {
    type = "IEFPCM, isotropic";
  } else {
    type = "IEFPCM, anisotropic";
  }
  os << "Solver Type: " << type << std::endl;
  if (hermitivitize_) {
    os << "PCM matrix hermitivitized";
  } else {
    os << "PCM matrix NOT hermitivitized (matches old DALTON)";
  }

  return os;
}

ISolver * createIEFSolver(const SolverData & data) {
  return new IEFSolver(data.hermitivitize);
}
} // namespace solver
} // namespace pcm
