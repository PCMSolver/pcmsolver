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

#ifndef BOUNDARYINTEGRALOPERATOR_HPP
#define BOUNDARYINTEGRALOPERATOR_HPP

#include "Config.hpp"

#include "IBoundaryIntegralOperator.hpp"
#include "BIOperatorData.hpp"
#include "Collocation.hpp"
#include "Purisima.hpp"
#include "Numerical.hpp"
#include "utils/Factory.hpp"

/*!
 * \file BoundaryIntegralOperator.hpp
 * \brief Top-level include file for boundary integral operators
 * \author Roberto Di Remigio
 * \date 2016
 *
 * Includes all boundary integral operator-related headers and defines the bootstrap
 *function
 * for the Factory<IBoundaryIntegralOperator, BIOperatorData>
 */

namespace pcm {
namespace bi_operators {

inline IBoundaryIntegralOperator * createCollocation(const BIOperatorData & data) {
  return new Collocation(data.scaling);
}

inline IBoundaryIntegralOperator * createNumerical(
    const BIOperatorData & /* data */) {
  return new Numerical();
}

inline IBoundaryIntegralOperator * createPurisima(const BIOperatorData & data) {
  return new Purisima(data.scaling);
}

inline void bootstrapFactory() {
  const bool registeredCollocation =
      Factory<IBoundaryIntegralOperator, BIOperatorData>::TheFactory()
          .registerObject("COLLOCATION", createCollocation);
  if (!registeredCollocation)
    PCMSOLVER_ERROR("Subscription of collocation integrator to factory failed!");

  const bool registeredNumerical =
      Factory<IBoundaryIntegralOperator, BIOperatorData>::TheFactory()
          .registerObject("NUMERICAL", createNumerical);
  if (!registeredNumerical)
    PCMSOLVER_ERROR("Subscription of numerical integrator to factory failed!");

  const bool registeredPurisima =
      Factory<IBoundaryIntegralOperator, BIOperatorData>::TheFactory()
          .registerObject("PURISIMA", createPurisima);
  if (!registeredPurisima)
    PCMSOLVER_ERROR("Subscription of Purisima integrator to factory failed!");
}
} // namespace bi_operators
} // namespace pcm

#endif // BOUNDARYINTEGRALOPERATOR_HPP
