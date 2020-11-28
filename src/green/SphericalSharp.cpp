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

#include "SphericalSharp.hpp"

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "cavity/Element.hpp"
#include "dielectric_profile/Sharp.hpp"
#include "utils/Legendre.hpp"
#include "utils/MathUtils.hpp"

namespace pcm {
using cavity::Element;
using dielectric_profile::Sharp;
namespace green {
template <typename DerivativeTraits>
SphericalSharp<DerivativeTraits>::SphericalSharp(double e,
                                                 double esolv,
                                                 double r,
                                                 const Eigen::Vector3d & o,
                                                 int l)
    : GreensFunction<DerivativeTraits, Sharp>(Sharp(e, esolv, r)),
      origin_(o),
      maxLGreen_(l) {}

template <typename DerivativeTraits>
double SphericalSharp<DerivativeTraits>::imagePotential(
    const Eigen::Vector3d & source,
    const Eigen::Vector3d & probe) const {
  DerivativeTraits sp[3], pp[3];
  sp[0] = source(0);
  sp[1] = source(1);
  sp[2] = source(2);
  pp[0] = probe(0);
  pp[1] = probe(1);
  pp[2] = probe(2);
  return this->imagePotential_impl(sp, pp)[0];
}

template <>
double SphericalSharp<Stencil>::imagePotential(const Eigen::Vector3d & source,
                                               const Eigen::Vector3d & probe) const {
  Stencil sp[3], pp[3];
  sp[0] = source(0);
  sp[1] = source(1);
  sp[2] = source(2);
  pp[0] = probe(0);
  pp[1] = probe(1);
  pp[2] = probe(2);
  return this->imagePotential_impl(sp, pp);
}

template <typename DerivativeTraits>
double SphericalSharp<DerivativeTraits>::imagePotentialDerivative(
    const Eigen::Vector3d & direction,
    const Eigen::Vector3d & p1,
    const Eigen::Vector3d & p2) const {
  DerivativeTraits t1[3], t2[3];
  t1[0] = p1(0);
  t1[1] = p1(1);
  t1[2] = p1(2);
  t2[0] = p2(0);
  t2[1] = p2(1);
  t2[2] = p2(2);
  t2[0][1] = direction(0);
  t2[1][1] = direction(1);
  t2[2][1] = direction(2);
  return this->imagePotential_impl(t1, t2)[1];
}

template <>
double SphericalSharp<Stencil>::imagePotentialDerivative(
    const Eigen::Vector3d & direction,
    const Eigen::Vector3d & p1,
    const Eigen::Vector3d & p2) const {
  return threePointStencil(std::bind(&SphericalSharp<Stencil>::imagePotential,
                                     this,
                                     std::placeholders::_1,
                                     std::placeholders::_2),
                           p2,
                           p1,
                           direction,
                           this->delta_);
}

template <typename DerivativeTraits>
DerivativeTraits SphericalSharp<DerivativeTraits>::operator()(
    DerivativeTraits * sp,
    DerivativeTraits * pp) const {
  DerivativeTraits distance =
      sqrt((sp[0] - pp[0]) * (sp[0] - pp[0]) + (sp[1] - pp[1]) * (sp[1] - pp[1]) +
           (sp[2] - pp[2]) * (sp[2] - pp[2]));
  return (1 / (this->profile_.epsilonSolvent * distance) -
          this->imagePotential_impl(sp, pp));
}

template <>
Stencil SphericalSharp<Stencil>::operator()(Stencil * sp, Stencil * pp) const {
  Stencil distance =
      sqrt((sp[0] - pp[0]) * (sp[0] - pp[0]) + (sp[1] - pp[1]) * (sp[1] - pp[1]) +
           (sp[2] - pp[2]) * (sp[2] - pp[2]));
  return (1 / (this->profile_.epsilonSolvent * distance) -
          this->imagePotential_impl(sp, pp));
}

template <typename DerivativeTraits>
double SphericalSharp<DerivativeTraits>::kernelD_impl(
    const Eigen::Vector3d & direction,
    const Eigen::Vector3d & p1,
    const Eigen::Vector3d & p2) const {
  return (this->profile_.epsilonSolvent * this->derivativeProbe(direction, p1, p2));
}

template <typename DerivativeTraits>
KernelS SphericalSharp<DerivativeTraits>::exportKernelS_impl() const {
  return std::bind(&SphericalSharp<DerivativeTraits>::kernelS,
                   *this,
                   std::placeholders::_1,
                   std::placeholders::_2);
}

template <typename DerivativeTraits>
KernelD SphericalSharp<DerivativeTraits>::exportKernelD_impl() const {
  return std::bind(&SphericalSharp<DerivativeTraits>::kernelD,
                   *this,
                   std::placeholders::_1,
                   std::placeholders::_2,
                   std::placeholders::_3);
}

template <typename DerivativeTraits>
DerivativeProbe SphericalSharp<DerivativeTraits>::exportDerivativeProbe_impl()
    const {
  return std::bind(&SphericalSharp<DerivativeTraits>::derivativeProbe,
                   *this,
                   std::placeholders::_1,
                   std::placeholders::_2,
                   std::placeholders::_3);
}

template <typename DerivativeTraits>
double SphericalSharp<DerivativeTraits>::singleLayer_impl(const Element & e,
                                                          double factor) const {
  // Diagonal of S inside the cavity
  double Sii_I = detail::diagonalSi(e.area(), factor);
  double image = this->imagePotential(e.center(), e.center());
  return (Sii_I / this->profile_.epsilonSolvent + image);
}

template <typename DerivativeTraits>
double SphericalSharp<DerivativeTraits>::doubleLayer_impl(const Element & e,
                                                          double factor) const {
  // Diagonal of D inside the cavity
  double Dii_I = detail::diagonalDi(e.area(), e.sphere().radius, factor);
  // "Diagonal" of the directional derivative of the image Green's function
  double image_grad =
      this->imagePotentialDerivative(e.normal(), e.center(), e.center());
  return (Dii_I + this->profile_.epsilonSolvent * image_grad);
}

template <typename DerivativeTraits>
DerivativeTraits SphericalSharp<DerivativeTraits>::imagePotential_impl(
    DerivativeTraits * sp,
    DerivativeTraits * pp) const {
  // Data from permittivity profile
  double radius = this->profile_.radius;
  double epsSolv = this->profile_.epsilonSolvent;
  double eps = this->profile_.epsilon;
  DerivativeTraits sp_origin[3], pp_origin[3];
  for (int i = 0; i < 3; ++i) {
    sp_origin[i] = sp[i] - origin_(i);
    pp_origin[i] = pp[i] - origin_(i);
  }
  DerivativeTraits sp_origin_norm = norm(sp_origin);
  DerivativeTraits pp_origin_norm = norm(pp_origin);
  // Angle between source and probe point
  DerivativeTraits cos_gamma =
      dot_product(sp_origin, pp_origin) / (sp_origin_norm * pp_origin_norm);

  DerivativeTraits r_img[3];
  for (int i = 0; i < 3; ++i) {
    r_img[i] = origin_(i) + pow(radius / pp_origin_norm, 2) * pp_origin[i];
  }
  // Distance between sp and r_img
  DerivativeTraits sp_image = distance(sp, r_img);
  DerivativeTraits q_img = radius / pp_origin_norm;
  double factor = (eps - epsSolv) / (eps + epsSolv);
  // DerivativeTraits G_img = DerivativeTraits(0.0);
  DerivativeTraits G_img = factor * (q_img / sp_image - q_img / sp_origin_norm);
  DerivativeTraits f_0 = radius / (sp_origin_norm * pp_origin_norm);
  DerivativeTraits f_L = f_0;
  for (int L = 1; L <= maxLGreen_; ++L) {
    f_L = f_L * radius * f_0;
    double C_0_L = (eps - epsSolv) * L / ((eps + epsSolv) * L + epsSolv);
    DerivativeTraits pl_x = Legendre::Pn<DerivativeTraits>(L, cos_gamma);
    G_img += f_L * (C_0_L - factor) * pl_x;
  }

  return G_img / epsSolv;
}

template <typename DerivativeTraits>
std::ostream & SphericalSharp<DerivativeTraits>::printObject(std::ostream & os) {
  Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, 0, ", ", "\n", "(", ")");
  os << "Green's function type: spherical sharp" << std::endl;
  os << this->profile_ << std::endl;
  os << "Sphere center        = " << this->origin_.transpose().format(CleanFmt)
     << std::endl;
  os << "Angular momentum (Green's function)    = " << this->maxLGreen_;
  return os;
}

template class SphericalSharp<Stencil>;
template class SphericalSharp<AD_directional>;
template class SphericalSharp<AD_gradient>;
template class SphericalSharp<AD_hessian>;
} // namespace green
} // namespace pcm
