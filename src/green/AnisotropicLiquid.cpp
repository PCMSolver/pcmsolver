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

#include "AnisotropicLiquid.hpp"

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

#include "DerivativeTypes.hpp"
#include "DerivativeUtils.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "cavity/Element.hpp"
#include "dielectric_profile/Anisotropic.hpp"
#include "utils/ForId.hpp"

namespace pcm {
using cavity::Element;
using dielectric_profile::Anisotropic;
namespace green {
template <typename DerivativeTraits>
AnisotropicLiquid<DerivativeTraits>::AnisotropicLiquid(
    const Eigen::Vector3d & eigen_eps,
    const Eigen::Vector3d & euler_ang)
    : GreensFunction<DerivativeTraits, Anisotropic>() {
  this->profile_ = Anisotropic(eigen_eps, euler_ang);
}

template <typename DerivativeTraits>
DerivativeTraits AnisotropicLiquid<DerivativeTraits>::operator()(
    DerivativeTraits * sp,
    DerivativeTraits * pp) const {
  // The distance has to be calculated using epsilonInv_ as metric:
  DerivativeTraits scratch = 0.0;
  Eigen::Matrix3d epsilonInv_ = this->profile_.epsilonInv();
  double detEps_ = this->profile_.detEps();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      scratch += (sp[i] - pp[i]) * epsilonInv_(i, j) * (sp[j] - pp[j]);
    }
  }
  DerivativeTraits distance = sqrt(scratch);

  return (1.0 / (sqrt(detEps_) * distance));
}

template <typename DerivativeTraits>
double AnisotropicLiquid<DerivativeTraits>::kernelD_impl(
    const Eigen::Vector3d & direction,
    const Eigen::Vector3d & p1,
    const Eigen::Vector3d & p2) const {
  // Since the permittivity is a tensorial quantity,
  // the full gradient is needed to get the kernel of D and D^\dagger
  Eigen::Vector3d scratch = this->profile_.epsilon() * (this->gradientProbe(p1, p2));
  return scratch.dot(direction);
}

template <typename DerivativeTraits>
KernelS AnisotropicLiquid<DerivativeTraits>::exportKernelS_impl() const {
  return pcm::bind(
      &AnisotropicLiquid<DerivativeTraits>::kernelS, *this, pcm::_1, pcm::_2);
}

template <typename DerivativeTraits>
KernelD AnisotropicLiquid<DerivativeTraits>::exportKernelD_impl() const {
  return pcm::bind(&AnisotropicLiquid<DerivativeTraits>::kernelD,
                   *this,
                   pcm::_1,
                   pcm::_2,
                   pcm::_3);
}

template <typename DerivativeTraits>
DerivativeProbe AnisotropicLiquid<DerivativeTraits>::exportDerivativeProbe_impl()
    const {
  return pcm::bind(&AnisotropicLiquid<DerivativeTraits>::derivativeProbe,
                   *this,
                   pcm::_1,
                   pcm::_2,
                   pcm::_3);
}

template <typename DerivativeTraits>
double AnisotropicLiquid<DerivativeTraits>::singleLayer_impl(
    const Element & /* e */,
    double /* factor */) const {
  PCMSOLVER_ERROR("Not implemented yet for AnisotropicLiquid");
  // return 0.0;
}

template <typename DerivativeTraits>
double AnisotropicLiquid<DerivativeTraits>::doubleLayer_impl(
    const Element & /* e */,
    double /* factor */) const {
  PCMSOLVER_ERROR("Not implemented yet for AnisotropicLiquid");
  // return 0.0;
}

template <typename DerivativeTraits>
std::ostream & AnisotropicLiquid<DerivativeTraits>::printObject(std::ostream & os) {
  os << "Green's function type: anisotropic liquid" << std::endl;
  os << this->profile_;
  return os;
}

template class AnisotropicLiquid<Stencil>;
template class AnisotropicLiquid<AD_directional>;
template class AnisotropicLiquid<AD_gradient>;
template class AnisotropicLiquid<AD_hessian>;

IGreensFunction * createAnisotropicLiquid(const GreenData & data) {
  detail::buildAnisotropicLiquid build;
  return for_id<derivative_types, IGreensFunction>(build, data, data.howDerivative);
}
} // namespace green
} // namespace pcm
