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

#include "Collocation.hpp"

#include "Config.hpp"

#include <Eigen/Core>

#include "BIOperatorData.hpp"
#include "cavity/Element.hpp"
#include "green/IGreensFunction.hpp"

namespace pcm {
using cavity::Element;

namespace bi_operators {
Collocation::Collocation() : factor_(1.07) {}

Collocation::Collocation(double fac) : factor_(fac) {}

Eigen::MatrixXd Collocation::computeS_impl(const std::vector<Element> & elems,
                                           const IGreensFunction & gf) const {
  PCMSolverIndex cavitySize = elems.size();
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
  for (PCMSolverIndex i = 0; i < cavitySize; ++i) {
    Element source = elems[i];
    S(i, i) = gf.singleLayer(source, factor_);
    for (PCMSolverIndex j = 0; j < cavitySize; ++j) {
      Element probe = elems[j];
      if (i != j)
        S(i, j) = gf.kernelS(source.center(), probe.center());
    }
  }
  return S;
}

Eigen::MatrixXd Collocation::computeD_impl(const std::vector<Element> & elems,
                                           const IGreensFunction & gf) const {
  PCMSolverIndex cavitySize = elems.size();
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
  for (PCMSolverIndex i = 0; i < cavitySize; ++i) {
    Element source = elems[i];
    D(i, i) = gf.doubleLayer(source, factor_);
    for (PCMSolverIndex j = 0; j < cavitySize; ++j) {
      Element probe = elems[j];
      if (i != j)
        D(i, j) =
            gf.kernelD(probe.normal().normalized(), source.center(), probe.center());
    }
  }
  return D;
}

IBoundaryIntegralOperator * createCollocation(const BIOperatorData & data) {
  return new Collocation(data.scaling);
}
} // namespace bi_operators
} // namespace pcm
