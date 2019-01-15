/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2019 Roberto Di Remigio, Luca Frediani and contributors.
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

namespace pcm {
namespace utils {
namespace detail {
typedef pcm::tuple<std::vector<double>, std::vector<double> > QuadratureRule;

QuadratureRule initializeRule(int nNodes);

QuadratureRule GaussLegendre16();

QuadratureRule GaussLegendre32();

QuadratureRule GaussLegendre64();
} // namespace detail

template <int nNodes> class GaussLegendreRule {
public:
  GaussLegendreRule()
      : nPoints_(nNodes), abscissa_((nPoints_ / 2), 0), weight_((nPoints_ / 2), 0) {
    pcm::tie(abscissa_, weight_) = detail::initializeRule(nPoints_);
  }
  int nPoints() const { return nPoints_; }
  double abscissa(int i) const { return abscissa_[i]; }
  double weight(int i) const { return weight_[i]; }

protected:
  int nPoints_;
  std::vector<double> abscissa_;
  std::vector<double> weight_;
};
} // namespace utils
} // namespace pcm
