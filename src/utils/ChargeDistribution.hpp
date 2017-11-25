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

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

/*! \file ChargeDistribution.hpp */

namespace pcm {
class Molecule;
} // namespace pcm

namespace pcm {
namespace utils {
/*!  \struct ChargeDistribution
 *  \brief POD representing a classical charge distribution
 *  \author Roberto Di Remigio
 *  \date 2016
 */
struct ChargeDistribution {
  /*! Monopoles */
  Eigen::VectorXd monopoles;
  /*! Monopoles sites */
  Eigen::Matrix3Xd monopolesSites;
  /*! Dipoles */
  Eigen::Matrix3Xd dipoles;
  /*! Dipoles sites */
  Eigen::Matrix3Xd dipolesSites;
};

/*! \typedef GFValue
 *  \brief functor handle to the calculation of the value of a Greens's function in a
 * point
 */
typedef pcm::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)>
    GFValue;

/*! \typedef GFDerivative
 *  \brief functor handle to the derivative of a Green's function in a point
 */
typedef pcm::function<double(const Eigen::Vector3d &,
                             const Eigen::Vector3d &,
                             const Eigen::Vector3d &)>
    GFDerivative;

/*! \brief Computes Newton potential in a set of points for given Green's function
 * and classical charge distribution
 *  \param[in] gf   the Green's function
 *  \param[in] grid where to evaluate the Newton potential
 *  \param[in] dist classical charge distribution
 *  \return the Newton potential on the grid
 */
Eigen::VectorXd computeNewtonPotential(const GFValue & gf,
                                       const Eigen::Matrix3Xd & grid,
                                       const ChargeDistribution & dist);

/*! \brief Computes dipolar potential in vacuum in a set of points and classical
 * charge distribution
 *  \param[in] grid where to evaluate the Newton potential
 *  \param[in] dist classical charge distribution
 *  \return the Newton potential on the grid
 */
Eigen::VectorXd computeDipolarPotential(const Eigen::Matrix3Xd & grid,
                                        const ChargeDistribution & dist);

/*! \brief Computes dipolar potential in a set of points for given Green's function
 * and classical charge distribution
 *  \param[in] gf   the Green's function derivative (dipolar interaction tensor)
 *  \param[in] grid where to evaluate the Newton potential
 *  \param[in] dist classical charge distribution
 *  \return the Newton potential on the grid
 */
Eigen::VectorXd computeDipolarPotential(const GFDerivative & gf,
                                        const Eigen::Matrix3Xd & grid,
                                        const ChargeDistribution & dist);

/*! \brief Return classical nuclear charge distribution from a Molecule object
 *  \param[in] mol the Molecule object
 *  \return the classical nuclear charge distribution
 */
ChargeDistribution nuclearChargeDistribution(const Molecule & mol);
} // namespace utils
} // namespace pcm
