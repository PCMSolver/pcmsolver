/**
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2016 Roberto Di Remigio, Luca Frediani and contributors.
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
#include "GreensFunction.hpp"
#include "cavity/Element.hpp"
#include "dielectric_profile/Yukawa.hpp"

#include "GreenData.hpp"
#include "utils/ForId.hpp"
#include "utils/Factory.hpp"

template <typename DerivativeTraits>
IonicLiquid<DerivativeTraits>::IonicLiquid(double eps, double k)
    : GreensFunction<DerivativeTraits, Yukawa>() {
  this->profile_ = Yukawa(eps, k);
}

template <typename DerivativeTraits>
DerivativeTraits IonicLiquid<DerivativeTraits>::operator()(
    DerivativeTraits * sp, DerivativeTraits * pp) const {
  double eps = this->profile_.epsilon;
  double k = this->profile_.kappa;
  return (exp(-k * distance(sp, pp)) / (eps * distance(sp, pp)));
}

template <typename DerivativeTraits>
double IonicLiquid<DerivativeTraits>::kernelD_impl(
    const Eigen::Vector3d & direction, const Eigen::Vector3d & p1,
    const Eigen::Vector3d & p2) const {
  return this->profile_.epsilon * (this->derivativeProbe(direction, p1, p2));
}

template <typename DerivativeTraits>
KernelS IonicLiquid<DerivativeTraits>::exportKernelS_impl() const {
  return pcm::bind(&IonicLiquid<DerivativeTraits>::kernelS, *this, pcm::_1, pcm::_2);
}

template <typename DerivativeTraits>
KernelD IonicLiquid<DerivativeTraits>::exportKernelD_impl() const {
  return pcm::bind(&IonicLiquid<DerivativeTraits>::kernelD, *this, pcm::_1, pcm::_2,
                   pcm::_3);
}

template <typename DerivativeTraits>
double IonicLiquid<DerivativeTraits>::singleLayer_impl(const Element & /* e */,
                                                       double /* factor */) const {
  PCMSOLVER_ERROR("Not implemented yet for IonicLiquid", BOOST_CURRENT_FUNCTION);
  return 0.0;
}

template <typename DerivativeTraits>
double IonicLiquid<DerivativeTraits>::doubleLayer_impl(const Element & /* e */,
                                                       double /* factor */) const {
  PCMSOLVER_ERROR("Not implemented yet for IonicLiquid", BOOST_CURRENT_FUNCTION);
  return 0.0;
}

template <typename DerivativeTraits>
std::ostream & IonicLiquid<DerivativeTraits>::printObject(std::ostream & os) {
  os << "Green's function type: ionic liquid" << std::endl;
  os << this->profile_;
  return os;
}

namespace {
struct buildIonicLiquid {
  template <typename T> IGreensFunction * operator()(const greenData & data) {
    return new IonicLiquid<T>(data.epsilon, data.kappa);
  }
};

IGreensFunction * createIonicLiquid(const greenData & data) {
  buildIonicLiquid build;
  return for_id<derivative_types, IGreensFunction>(build, data, data.howDerivative);
}
const std::string IONICLIQUID("IONICLIQUID");
const bool registeredIonicLiquid =
    Factory<IGreensFunction, greenData>::TheFactory().registerObject(
        IONICLIQUID, createIonicLiquid);
}

template class IonicLiquid<Numerical>;
template class IonicLiquid<AD_directional>;
template class IonicLiquid<AD_gradient>;
template class IonicLiquid<AD_hessian>;
