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

#include "UniformDielectric.hpp"

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

#include "DerivativeTypes.hpp"
#include "DerivativeUtils.hpp"
#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "cavity/Element.hpp"
#include "dielectric_profile/Uniform.hpp"
#include "utils/ForId.hpp"

namespace pcm {
using cavity::Element;
using dielectric_profile::Uniform;
namespace green {
template <typename DerivativeTraits>
UniformDielectric<DerivativeTraits>::UniformDielectric(double eps)
    : GreensFunction<DerivativeTraits, Uniform>() {
  this->profile_ = Uniform(eps);
}

template <typename DerivativeTraits>
DerivativeTraits UniformDielectric<DerivativeTraits>::operator()(
    DerivativeTraits * sp,
    DerivativeTraits * pp) const {
  return 1 / (this->profile_.epsilon * distance(sp, pp));
}

template <typename DerivativeTraits>
double UniformDielectric<DerivativeTraits>::kernelD_impl(
    const Eigen::Vector3d & direction,
    const Eigen::Vector3d & p1,
    const Eigen::Vector3d & p2) const {
  return this->profile_.epsilon * (this->derivativeProbe(direction, p1, p2));
}

template <typename DerivativeTraits>
KernelS UniformDielectric<DerivativeTraits>::exportKernelS_impl() const {
  return pcm::bind(
      &UniformDielectric<DerivativeTraits>::kernelS, *this, pcm::_1, pcm::_2);
}

template <typename DerivativeTraits>
KernelD UniformDielectric<DerivativeTraits>::exportKernelD_impl() const {
  return pcm::bind(&UniformDielectric<DerivativeTraits>::kernelD,
                   *this,
                   pcm::_1,
                   pcm::_2,
                   pcm::_3);
}

template <typename DerivativeTraits>
DerivativeProbe UniformDielectric<DerivativeTraits>::exportDerivativeProbe_impl()
    const {
  return pcm::bind(&UniformDielectric<DerivativeTraits>::derivativeProbe,
                   *this,
                   pcm::_1,
                   pcm::_2,
                   pcm::_3);
}

template <typename DerivativeTraits>
double UniformDielectric<DerivativeTraits>::singleLayer_impl(const Element & e,
                                                             double factor) const {
  return (detail::diagonalSi(e.area(), factor) / this->profile_.epsilon);
}

template <typename DerivativeTraits>
double UniformDielectric<DerivativeTraits>::doubleLayer_impl(const Element & e,
                                                             double factor) const {
  return detail::diagonalDi(e.area(), e.sphere().radius, factor);
}

template <typename DerivativeTraits>
std::ostream & UniformDielectric<DerivativeTraits>::printObject(std::ostream & os) {
  os << "Green's function type: uniform dielectric" << std::endl;
  os << this->profile_;
  return os;
}

template class UniformDielectric<Stencil>;
template class UniformDielectric<AD_directional>;
template class UniformDielectric<AD_gradient>;
template class UniformDielectric<AD_hessian>;

IGreensFunction * createUniformDielectric(const GreenData & data) {
  detail::buildUniformDielectric build;
  return for_id<derivative_types, IGreensFunction>(build, data, data.howDerivative);
}
} // namespace green
} // namespace pcm
