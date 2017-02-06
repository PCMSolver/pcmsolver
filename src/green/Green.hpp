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

#ifndef GREEN_HPP
#define GREEN_HPP

#include "Config.hpp"

#include "IGreensFunction.hpp"
#include "GreenData.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "AnisotropicLiquid.hpp"
#include "IonicLiquid.hpp"
#include "SphericalDiffuse.hpp"
#include "utils/ForId.hpp"
#include "utils/Factory.hpp"

/*!
 * \file Solver.hpp
 * \brief Top-level include file for solvers
 * \author Roberto Di Remigio
 * \date 2016
 *
 * Includes all solver-related headers and defines the bootstrap function
 * for the Factory<ISolver, SolverData>
 */

namespace pcm {
namespace green {
inline IGreensFunction * createVacuum(const GreenData & data) {
  detail::buildVacuum build;
  return for_id<derivative_types, IGreensFunction>(build, data, data.howDerivative);
}

inline IGreensFunction * createUniformDielectric(const GreenData & data) {
  detail::buildUniformDielectric build;
  return for_id<derivative_types, IGreensFunction>(build, data, data.howDerivative);
}

inline IGreensFunction * createAnisotropicLiquid(const GreenData & data) {
  detail::buildAnisotropicLiquid build;
  return for_id<derivative_types, IGreensFunction>(build, data, data.howDerivative);
}

inline IGreensFunction * createIonicLiquid(const GreenData & data) {
  detail::buildIonicLiquid build;
  return for_id<derivative_types, IGreensFunction>(build, data, data.howDerivative);
}

inline IGreensFunction * createSphericalDiffuse(const GreenData & data) {
  detail::buildSphericalDiffuse build;
  return for_id<dielectric_profile::onelayer_diffuse_profile_types, IGreensFunction>(
      build, data, data.howProfile);
}

inline void bootstrapFactory() {
  const bool registeredVacuum =
      Factory<IGreensFunction, GreenData>::TheFactory().registerObject("VACUUM",
                                                                       createVacuum);
  if (!registeredVacuum)
    PCMSOLVER_ERROR("Subscription of vacuum Green's function to factory failed!");

  const bool registeredUniformDielectric =
      Factory<IGreensFunction, GreenData>::TheFactory().registerObject(
          "UNIFORMDIELECTRIC", createUniformDielectric);
  if (!registeredUniformDielectric)
    PCMSOLVER_ERROR(
        "Subscription of uniform dielectric Green's function to factory failed!");

  const bool registeredSphericalDiffuse =
      Factory<IGreensFunction, GreenData>::TheFactory().registerObject(
          "SPHERICALDIFFUSE", createSphericalDiffuse);
  if (!registeredSphericalDiffuse)
    PCMSOLVER_ERROR(
        "Subscription of spherical diffuse Green's function to factory failed!");

  const bool registeredIonicLiquid =
      Factory<IGreensFunction, GreenData>::TheFactory().registerObject(
          "IONICLIQUID", createIonicLiquid);
  if (!registeredIonicLiquid)
    PCMSOLVER_ERROR(
        "Subscription of ionic liquid Green's function to factory failed!");

  const bool registeredAnisotropicLiquid =
      Factory<IGreensFunction, GreenData>::TheFactory().registerObject(
          "ANISOTROPICLIQUID", createAnisotropicLiquid);
  if (!registeredAnisotropicLiquid)
    PCMSOLVER_ERROR(
        "Subscription of anisotropic liquid Green's function to factory failed!");
}
} // namespace green
} // namespace pcm

#endif // GREEN_HPP
