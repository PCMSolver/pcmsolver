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

#ifndef MMFQ_HPP
#define MMFQ_HPP

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

namespace pcm {
namespace utils {
/*! \file MMFQ.hpp
 *  \struct MMFQ
 *  \brief POD representing a classical fluctuating charges MM force field
 *  \author Roberto Di Remigio
 *  \date 2017
 */
struct MMFQ {
  /*! Number of FQ fragments */
  PCMSolverIndex nFragments;
  /*! Number of FQ sites per MM fragment
   *  \note This would be three for a standard representation of water
   */
  PCMSolverIndex nSitesPerFragment;
  /*! FQ electronegativities */
  Eigen::VectorXd chi;
  /*! FQ hardnesses */
  Eigen::VectorXd eta;
  /*! FQ sites */
  Eigen::Matrix3Xd sites;
};
} // namespace utils
} // namespace pcm

#endif // MMFQ_HPP
