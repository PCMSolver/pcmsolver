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

/*! \file UniformDielectric.hpp */

namespace pcm {
namespace cavity {
class Element;
} // namespace cavity
} // namespace pcm

#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "dielectric_profile/Uniform.hpp"

namespace pcm {
namespace green {
/*! \class UniformDielectric
 *  \brief Green's function for uniform dielectric.
 *  \author Luca Frediani and Roberto Di Remigio
 *  \date 2012-2016
 *  \tparam DerivativeTraits evaluation strategy for the function and its derivatives
 */
template <typename DerivativeTraits = AD_directional>
class UniformDielectric final
    : public GreensFunction<DerivativeTraits, dielectric_profile::Uniform> {
public:
  UniformDielectric(double eps);

  virtual double permittivity() const override final {
    return this->profile_.epsilon;
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

  virtual double singleLayer_impl(const Element & e, double factor) const override;
  virtual double doubleLayer_impl(const Element & e, double factor) const override;

  virtual std::ostream & printObject(std::ostream & os) override;
};

template <typename DerivativeTraits>
IGreensFunction * createUniformDielectric(const GreenData & data) {
  return new UniformDielectric<DerivativeTraits>(data.epsilon);
}
} // namespace green
} // namespace pcm
