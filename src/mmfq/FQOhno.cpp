/**
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

#include "FQOhno.hpp"

#include <iostream>

#include "Config.hpp"

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include "utils/MMFQ.hpp"

namespace pcm {
namespace mmfq {
using utils::MMFQ;

FQOhno::FQOhno(const MMFQ & ff, bool nonpol) : mmfq_(ff), nonPolarizable_(nonpol) {
  if (!nonPolarizable_) buildSystemMatrix_impl();
}

void FQOhno::buildSystemMatrix_impl() {
  PCMSolverIndex nFragments = mmfq_.nFragments;
  PCMSolverIndex nSitesPerFragment = mmfq_.nSitesPerFragment;
  PCMSolverIndex nSites = nFragments * nSitesPerFragment;
  PCMSolverIndex fq_dim = nSites + nFragments;
  TIMER_ON("Computing Dlambda");
  Dlambda_ = Eigen::MatrixXd::Zero(fq_dim, fq_dim);
  // Fill J block
  for (PCMSolverIndex i = 0; i < nSites; ++i) {
    for (PCMSolverIndex j = 0; j < nSites; ++j) {
      if (i == j) {
        Dlambda_(i, i) = mmfq_.eta(i);
      } else {
        double eta_ij = 0.5 * (mmfq_.eta(i) + mmfq_.eta(j));
        double dist = (mmfq_.sites.col(i) - mmfq_.sites.col(j)).norm();
        Dlambda_(i, j) = eta_ij / std::sqrt(1.0 + std::pow(eta_ij * dist, 2));
      }
    }
  }
  for (PCMSolverIndex i = 0; i < nFragments; ++i) {
    for (PCMSolverIndex j = 0; j < nSitesPerFragment; ++j) {
      Dlambda_(i + nSites, j + i * nSitesPerFragment) = 1.0;
      Dlambda_(j + i * nSitesPerFragment, i + nSites) = 1.0;
    }
  }
  TIMER_OFF("Computing Dlambda");

  built_ = true;
}

Eigen::VectorXd FQOhno::computeCharge_impl(const Eigen::VectorXd & potential,
                                           bool scf) const {
  Eigen::VectorXd RHS = Eigen::VectorXd::Zero(Dlambda_.rows());
  RHS.head(potential.size()) = potential;
  // If doing SCF, we have to add electronegativities on top of the potential
  if (scf)
    RHS.head(potential.size()) += mmfq_.chi;
  return -Dlambda_.ldlt().solve(RHS).head(potential.size());
}

std::ostream & FQOhno::printSolver(std::ostream & os) {
  os << "Fluctuating charge solver type: Ohno" << std::endl;
  if (nonPolarizable_) os << "Nonpolarizable force field" << std::endl;
  os << "Number of fragments = " << mmfq_.nFragments << std::endl;
  os << "Number of sites per fragment = " << mmfq_.nSitesPerFragment << std::endl;
  os << "Number of sites = " << mmfq_.nFragments * mmfq_.nSitesPerFragment;
  return os;
}
} // namespace mmfq
} // namespace pcm
