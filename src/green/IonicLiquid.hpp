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

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

/*! \file IonicLiquid.hpp */

namespace pcm {
namespace cavity {
class Element;
} // namespace cavity
} // namespace pcm

#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "dielectric_profile/Yukawa.hpp"

namespace pcm {
namespace green {
/*! \class IonicLiquid
 *  \brief Green's functions for ionic liquid, described by the linearized
 * Poisson-Boltzmann equation.
 *  \author Luca Frediani, Roberto Di Remigio
 *  \date 2013-2016
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 */
template <typename DerivativeTraits = AD_directional>
class IonicLiquid final
    : public GreensFunction<DerivativeTraits, dielectric_profile::Yukawa> {
public:
  IonicLiquid(double eps, double k);

  virtual double permittivity() const override final {
    PCMSOLVER_ERROR("permittivity() only implemented for uniform dielectrics");
  }

private:
  virtual DerivativeTraits operator()(DerivativeTraits * sp,
                                      DerivativeTraits * pp) const override;
  virtual double kernelD_impl(const Eigen::Vector3d & direction,
                              const Eigen::Vector3d & p1,
                              const Eigen::Vector3d & p2) const override;

  virtual KernelS exportKernelS_impl() const override;
  virtual KernelD exportKernelD_impl() const override;
  virtual DerivativeProbe exportDerivativeProbe_impl() const override;

  [[noreturn]] virtual double singleLayer_impl(const Element & /* e */,
                                               double /* factor */) const override;
  [[noreturn]] virtual double doubleLayer_impl(const Element & /* e */,
                                               double /* factor */) const override;

  virtual std::ostream & printObject(std::ostream & os) override;
};

template <typename DerivativeTraits>
IGreensFunction * createIonicLiquid(const GreenData & data) {
  return new IonicLiquid<DerivativeTraits>(data.epsilon, data.kappa);
}
} // namespace green
} // namespace pcm
