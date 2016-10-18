/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2016 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */
/* pcmsolver_copyright_end */

#include "CPCMSolver.hpp"

#include <fstream>
#include <iostream>
#include <string>

#include "Config.hpp"

#include <Eigen/Core>

#include "bi_operators/BoundaryIntegralOperator.hpp"
#include "cavity/Cavity.hpp"
#include "cavity/Element.hpp"
#include "green/IGreensFunction.hpp"
#include "utils/MathUtils.hpp"
#include "SolverData.hpp"
#include "utils/Factory.hpp"

void CPCMSolver::buildSystemMatrix_impl(const Cavity & cavity,
                                        const IGreensFunction & gf_i,
                                        const IGreensFunction & gf_o,
                                        const BoundaryIntegralOperator & op) {
  if (!isotropic_)
    PCMSOLVER_ERROR("C-PCM is defined only for isotropic environments!",
                    BOOST_CURRENT_FUNCTION);
  TIMER_ON("Computing S");
  double epsilon = profiles::epsilon(gf_o.permittivity());
  S_ = op.computeS(cavity, gf_i);
  S_ /= (epsilon - 1.0) / (epsilon + correction_);
  // Get in Hermitian form
  if (hermitivitize_)
    hermitivitize(S_);
  TIMER_OFF("Computing S");

  // Symmetry-pack
  // The number of irreps in the group
  int nrBlocks = cavity.pointGroup().nrIrrep();
  // The size of the irreducible portion of the cavity
  int dimBlock = cavity.irreducible_size();
  symmetryPacking(blockS_, S_, dimBlock, nrBlocks);

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
      -blockS_[irrep].llt().solve(potential.segment(irrep * irrDim, irrDim));

  return charge;
}

std::ostream & CPCMSolver::printSolver(std::ostream & os) {
  os << "Solver Type: C-PCM" << std::endl;
  if (hermitivitize_) {
    os << "PCM matrix hermitivitized";
  } else {
    os << "PCM matrix NOT hermitivitized (matches old DALTON)";
  }

  return os;
}

namespace {
PCMSolver * createCPCMSolver(const solverData & data) {
  return new CPCMSolver(data.hermitivitize, data.correction);
}
const std::string CPCMSOLVER("CPCM");
const bool registeredCPCMSolver =
    Factory<PCMSolver, solverData>::TheFactory().registerObject(CPCMSOLVER,
                                                                createCPCMSolver);
}
