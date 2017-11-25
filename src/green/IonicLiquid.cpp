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

#include "IonicLiquid.hpp"

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

#include "DerivativeTypes.hpp"
#include "DerivativeUtils.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "cavity/Element.hpp"
#include "dielectric_profile/Yukawa.hpp"
#include "utils/ForId.hpp"

namespace pcm {
using cavity::Element;
using dielectric_profile::Yukawa;
namespace green {
template <typename DerivativeTraits>
IonicLiquid<DerivativeTraits>::IonicLiquid(double eps, double k)
    : GreensFunction<DerivativeTraits, Yukawa>() {
  this->profile_ = Yukawa(eps, k);
}

template <typename DerivativeTraits>
DerivativeTraits IonicLiquid<DerivativeTraits>::operator()(
    DerivativeTraits * sp,
    DerivativeTraits * pp) const {
  double eps = this->profile_.epsilon;
  double k = this->profile_.kappa;
  return (exp(-k * distance(sp, pp)) / (eps * distance(sp, pp)));
}

template <typename DerivativeTraits>
double IonicLiquid<DerivativeTraits>::kernelD_impl(
    const Eigen::Vector3d & direction,
    const Eigen::Vector3d & p1,
    const Eigen::Vector3d & p2) const {
  return this->profile_.epsilon * (this->derivativeProbe(direction, p1, p2));
}

template <typename DerivativeTraits>
KernelS IonicLiquid<DerivativeTraits>::exportKernelS_impl() const {
  return pcm::bind(&IonicLiquid<DerivativeTraits>::kernelS, *this, pcm::_1, pcm::_2);
}

template <typename DerivativeTraits>
KernelD IonicLiquid<DerivativeTraits>::exportKernelD_impl() const {
  return pcm::bind(
      &IonicLiquid<DerivativeTraits>::kernelD, *this, pcm::_1, pcm::_2, pcm::_3);
}

template <typename DerivativeTraits>
DerivativeProbe IonicLiquid<DerivativeTraits>::exportDerivativeProbe_impl() const {
  return pcm::bind(&IonicLiquid<DerivativeTraits>::derivativeProbe,
                   *this,
                   pcm::_1,
                   pcm::_2,
                   pcm::_3);
}

template <typename DerivativeTraits>
double IonicLiquid<DerivativeTraits>::singleLayer_impl(const Element & /* e */,
                                                       double /* factor */) const {
  PCMSOLVER_ERROR("Not implemented yet for IonicLiquid");
  // return 0.0;
}

template <typename DerivativeTraits>
double IonicLiquid<DerivativeTraits>::doubleLayer_impl(const Element & /* e */,
                                                       double /* factor */) const {
  PCMSOLVER_ERROR("Not implemented yet for IonicLiquid");
  // return 0.0;
}

template <typename DerivativeTraits>
std::ostream & IonicLiquid<DerivativeTraits>::printObject(std::ostream & os) {
  os << "Green's function type: ionic liquid" << std::endl;
  os << this->profile_;
  return os;
}

template class IonicLiquid<Stencil>;
template class IonicLiquid<AD_directional>;
template class IonicLiquid<AD_gradient>;
template class IonicLiquid<AD_hessian>;

IGreensFunction * createIonicLiquid(const GreenData & data) {
  detail::buildIonicLiquid build;
  return for_id<derivative_types, IGreensFunction>(build, data, data.howDerivative);
}
} // namespace green
} // namespace pcm
