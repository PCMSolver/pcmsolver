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

#include "ChargeDistribution.hpp"
#include <iostream>

#include "Config.hpp"

#include <Eigen/Core>

#include "Molecule.hpp"

namespace pcm {
namespace utils {
Eigen::VectorXd computeNewtonPotential(const GFValue & gf,
                                       const Eigen::Matrix3Xd & grid,
                                       const ChargeDistribution & dist) {
  Eigen::VectorXd newton = Eigen::VectorXd::Zero(grid.cols());
  for (int i = 0; i < dist.monopoles.size(); ++i) {
    for (int j = 0; j < grid.cols(); ++j) {
      newton(j) += dist.monopoles(i) * gf(grid.col(j), dist.monopolesSites.col(i));
    }
  }
  return newton;
}

Eigen::VectorXd computeDipolarPotential(const GFDerivative & gf,
                                        const Eigen::Matrix3Xd & grid,
                                        const ChargeDistribution & dist) {
  Eigen::VectorXd retval = Eigen::VectorXd::Zero(grid.cols());
  for (int i = 0; i < dist.dipoles.cols(); ++i) {
    for (int j = 0; j < grid.cols(); ++j) {
      retval(j) += gf(dist.dipoles.col(i), grid.col(j), dist.dipolesSites.col(i));
    }
  }
  return retval;
}

Eigen::VectorXd computeDipolarPotential(const Eigen::Matrix3Xd & grid,
                                        const ChargeDistribution & dist) {
  Eigen::VectorXd retval = Eigen::VectorXd::Zero(grid.cols());
  for (int i = 0; i < dist.dipoles.cols(); ++i) {
    for (int j = 0; j < grid.cols(); ++j) {
      Eigen::Vector3d distance = grid.col(j) - dist.dipolesSites.col(i);
      retval(j) +=
          (distance.dot(dist.dipoles.col(i))) / std::pow(distance.norm(), 3);
    }
  }
  return retval;
}

ChargeDistribution nuclearChargeDistribution(const Molecule & mol) {
  ChargeDistribution chg;
  chg.monopoles = mol.charges();
  chg.monopolesSites = mol.geometry();
  return chg;
}
} // namespace utils
} // namespace pcm
