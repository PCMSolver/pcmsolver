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

#include "Config.hpp"

#include "AnisotropicLiquid.hpp"
#include "DerivativeTypes.hpp"
#include "GreenData.hpp"
#include "IGreensFunction.hpp"
#include "IonicLiquid.hpp"
#include "SphericalDiffuse.hpp"
#include "SphericalSharp.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"
#include "dielectric_profile/OneLayerErf.hpp"
#include "dielectric_profile/OneLayerTanh.hpp"
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
namespace detail {
typedef pcm::function<IGreensFunction *(const GreenData &)> CreateGreensFunction;
} // namespace detail

inline Factory<detail::CreateGreensFunction> bootstrapFactory() {
  Factory<detail::CreateGreensFunction> factory_;

  factory_.subscribe("VACUUM_NUMERICAL", createVacuum<Stencil>);
  factory_.subscribe("VACUUM_DERIVATIVE", createVacuum<AD_directional>);
  factory_.subscribe("VACUUM_GRADIENT", createVacuum<AD_gradient>);
  factory_.subscribe("VACUUM_HESSIAN", createVacuum<AD_hessian>);

  factory_.subscribe("UNIFORMDIELECTRIC_NUMERICAL",
                     createUniformDielectric<Stencil>);
  factory_.subscribe("UNIFORMDIELECTRIC_DERIVATIVE",
                     createUniformDielectric<AD_directional>);
  factory_.subscribe("UNIFORMDIELECTRIC_GRADIENT",
                     createUniformDielectric<AD_gradient>);
  factory_.subscribe("UNIFORMDIELECTRIC_HESSIAN",
                     createUniformDielectric<AD_hessian>);

  factory_.subscribe("IONICLIQUID_NUMERICAL", createIonicLiquid<Stencil>);
  factory_.subscribe("IONICLIQUID_DERIVATIVE", createIonicLiquid<AD_directional>);
  factory_.subscribe("IONICLIQUID_GRADIENT", createIonicLiquid<AD_gradient>);
  factory_.subscribe("IONICLIQUID_HESSIAN", createIonicLiquid<AD_hessian>);

  factory_.subscribe("ANISOTROPICLIQUID_NUMERICAL",
                     createAnisotropicLiquid<Stencil>);
  factory_.subscribe("ANISOTROPICLIQUID_DERIVATIVE",
                     createAnisotropicLiquid<AD_directional>);
  factory_.subscribe("ANISOTROPICLIQUID_GRADIENT",
                     createAnisotropicLiquid<AD_gradient>);
  factory_.subscribe("ANISOTROPICLIQUID_HESSIAN",
                     createAnisotropicLiquid<AD_hessian>);

  factory_.subscribe("SPHERICALDIFFUSE_NUMERICAL_TANH",
                     createSphericalDiffuse<dielectric_profile::OneLayerTanh>);
  factory_.subscribe("SPHERICALDIFFUSE_NUMERICAL_ERF",
                     createSphericalDiffuse<dielectric_profile::OneLayerErf>);
  factory_.subscribe("SPHERICALDIFFUSE_NUMERICAL_LOG",
                     createSphericalDiffuse<dielectric_profile::OneLayerLog>);

  factory_.subscribe("SPHERICALSHARP_NUMERICAL", createSphericalSharp<Stencil>);
  factory_.subscribe("SPHERICALSHARP_DERIVATIVE",
                     createSphericalSharp<AD_directional>);
  factory_.subscribe("SPHERICALSHARP_GRADIENT", createSphericalSharp<AD_gradient>);
  factory_.subscribe("SPHERICALSHARP_HESSIAN", createSphericalSharp<AD_hessian>);

  return factory_;
}
} // namespace green
} // namespace pcm
