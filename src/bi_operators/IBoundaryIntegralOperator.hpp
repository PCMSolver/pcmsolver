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

#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

namespace pcm {
class ICavity;
namespace cavity {
class Element;
} // namespace cavity
class IGreensFunction;
} // namespace pcm

namespace pcm {
class IBoundaryIntegralOperator {
public:
  virtual ~IBoundaryIntegralOperator() {}
  Eigen::MatrixXd computeS(const ICavity & cav, const IGreensFunction & gf) const;
  Eigen::MatrixXd computeD(const ICavity & cav, const IGreensFunction & gf) const;

private:
  virtual Eigen::MatrixXd computeS_impl(const std::vector<cavity::Element> & elems,
                                        const IGreensFunction & gf) const = 0;
  virtual Eigen::MatrixXd computeD_impl(const std::vector<cavity::Element> & elems,
                                        const IGreensFunction & gf) const = 0;
};
} // namespace pcm
