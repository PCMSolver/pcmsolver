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

#include "SphericalDiffuse.hpp"

#include <iosfwd>

#include "Config.hpp"

#include <Eigen/Core>

// Has to be included here
#include "InterfacesImpl.hpp"

#include "GreenData.hpp"
#include "GreensFunction.hpp"
#include "cavity/Element.hpp"
#include "utils/Legendre.hpp"
#include "utils/MathUtils.hpp"

namespace pcm {
namespace green {
template <typename ProfilePolicy>
double SphericalDiffuse<ProfilePolicy>::coefficientCoulomb(
    const Eigen::Vector3d & source,
    const Eigen::Vector3d & probe) const {
  // Obtain coefficient for the separation of the Coulomb singularity
  return this->coefficient_impl(source, probe);
}

template <typename ProfilePolicy>
double SphericalDiffuse<ProfilePolicy>::Coulomb(
    const Eigen::Vector3d & source,
    const Eigen::Vector3d & probe) const {
  double r12 = (source - probe).norm();

  // Obtain coefficient for the separation of the Coulomb singularity
  return (1.0 / (this->coefficient_impl(source, probe) * r12));
}

template <typename ProfilePolicy>
double SphericalDiffuse<ProfilePolicy>::imagePotential(
    const Eigen::Vector3d & source,
    const Eigen::Vector3d & probe) const {
  // Obtain coefficient for the separation of the Coulomb singularity
  double Cr12 = this->coefficient_impl(source, probe);

  double gr12 = 0.0;
  for (int L = 1; L <= maxLGreen_; ++L) {
    gr12 += this->imagePotentialComponent_impl(L, source, probe, Cr12);
  }

  return gr12;
}

template <typename ProfilePolicy>
double SphericalDiffuse<ProfilePolicy>::coefficientCoulombDerivative(
    const Eigen::Vector3d & direction,
    const Eigen::Vector3d & p1,
    const Eigen::Vector3d & p2) const {
  return threePointStencil(
      pcm::bind(&SphericalDiffuse<ProfilePolicy>::coefficientCoulomb,
                this,
                pcm::_1,
                pcm::_2),
      p2,
      p1,
      direction,
      this->delta_);
}

template <typename ProfilePolicy>
double SphericalDiffuse<ProfilePolicy>::CoulombDerivative(
    const Eigen::Vector3d & direction,
    const Eigen::Vector3d & p1,
    const Eigen::Vector3d & p2) const {
  return threePointStencil(
      pcm::bind(&SphericalDiffuse<ProfilePolicy>::Coulomb, this, pcm::_1, pcm::_2),
      p2,
      p1,
      direction,
      this->delta_);
}

template <typename ProfilePolicy>
double SphericalDiffuse<ProfilePolicy>::imagePotentialDerivative(
    const Eigen::Vector3d & direction,
    const Eigen::Vector3d & p1,
    const Eigen::Vector3d & p2) const {
  return threePointStencil(
      pcm::bind(
          &SphericalDiffuse<ProfilePolicy>::imagePotential, this, pcm::_1, pcm::_2),
      p2,
      p1,
      direction,
      this->delta_);
}

template <typename ProfilePolicy>
pcm::tuple<double, double> SphericalDiffuse<ProfilePolicy>::epsilon(
    const Eigen::Vector3d & point) const {
  return this->profile_((point + this->origin_).norm());
}

template <typename ProfilePolicy>
void SphericalDiffuse<ProfilePolicy>::toFile(const std::string & prefix) {
  std::string tmp;
  prefix.empty() ? tmp = prefix : tmp = prefix + "-";
  writeToFile(zetaC_, tmp + "zetaC.dat");
  writeToFile(omegaC_, tmp + "omegaC.dat");
  for (int L = 1; L <= maxLGreen_; ++L) {
    writeToFile(zeta_[L], tmp + "zeta_" + pcm::to_string(L) + ".dat");
    writeToFile(omega_[L], tmp + "omega_" + pcm::to_string(L) + ".dat");
  }
}

template <typename ProfilePolicy>
Stencil SphericalDiffuse<ProfilePolicy>::operator()(Stencil * sp,
                                                    Stencil * pp) const {
  // Transfer raw arrays to Eigen vectors using the Map type
  Eigen::Map<Eigen::Matrix<double, 3, 1> > source(sp), probe(pp);

  // Obtain coefficient for the separation of the Coulomb singularity
  double Cr12 = this->coefficient_impl(source, probe);

  double gr12 = 0.0;
  for (int L = 0; L <= maxLGreen_; ++L) {
    gr12 += this->imagePotentialComponent_impl(L, source, probe, Cr12);
  }
  double r12 = (source - probe).norm();

  double gf = (1.0 / (Cr12 * r12) + gr12);
  return gf;
}

template <typename ProfilePolicy>
double SphericalDiffuse<ProfilePolicy>::kernelD_impl(
    const Eigen::Vector3d & direction,
    const Eigen::Vector3d & p1,
    const Eigen::Vector3d & p2) const {
  double eps_r2 = 0.0;
  // Shift p2 by origin_
  pcm::tie(eps_r2, pcm::ignore) = this->epsilon(p2);

  return (eps_r2 * this->derivativeProbe(direction, p1, p2));
}

template <typename ProfilePolicy>
KernelS SphericalDiffuse<ProfilePolicy>::exportKernelS_impl() const {
  return pcm::bind(
      &SphericalDiffuse<ProfilePolicy>::kernelS, *this, pcm::_1, pcm::_2);
}

template <typename ProfilePolicy>
KernelD SphericalDiffuse<ProfilePolicy>::exportKernelD_impl() const {
  return pcm::bind(
      &SphericalDiffuse<ProfilePolicy>::kernelD, *this, pcm::_1, pcm::_2, pcm::_3);
}

template <typename DerivativeTraits>
DerivativeProbe SphericalDiffuse<DerivativeTraits>::exportDerivativeProbe_impl()
    const {
  return pcm::bind(&SphericalDiffuse<DerivativeTraits>::derivativeProbe,
                   *this,
                   pcm::_1,
                   pcm::_2,
                   pcm::_3);
}

template <typename ProfilePolicy>
double SphericalDiffuse<ProfilePolicy>::singleLayer_impl(const Element & e,
                                                         double factor) const {
  // Diagonal of S inside the cavity
  double Sii_I = detail::diagonalSi(e.area(), factor);
  // "Diagonal" of Coulomb singularity separation coefficient
  double coulomb_coeff = this->coefficientCoulomb(e.center(), e.center());
  // "Diagonal" of the image Green's function
  double image = this->imagePotential(e.center(), e.center());
  return (Sii_I / coulomb_coeff + image);
}

template <typename ProfilePolicy>
double SphericalDiffuse<ProfilePolicy>::doubleLayer_impl(const Element & e,
                                                         double factor) const {
  // Diagonal of S inside the cavity
  double Sii_I = detail::diagonalSi(e.area(), factor);
  // Diagonal of D inside the cavity
  double Dii_I = detail::diagonalDi(e.area(), e.sphere().radius, factor);
  // "Diagonal" of Coulomb singularity separation coefficient
  double coulomb_coeff = this->coefficientCoulomb(e.center(), e.center());
  // "Diagonal" of the directional derivative of the Coulomb singularity
  // separation coefficient
  double coeff_grad =
      this->coefficientCoulombDerivative(e.normal(), e.center(), e.center()) /
      std::pow(coulomb_coeff, 2);
  // "Diagonal" of the directional derivative of the image Green's function
  double image_grad =
      this->imagePotentialDerivative(e.normal(), e.center(), e.center());

  double eps_r2 = 0.0;
  pcm::tie(eps_r2, pcm::ignore) = this->epsilon(e.center());

  return (eps_r2 * (Dii_I / coulomb_coeff - Sii_I * coeff_grad + image_grad));
}

template <typename ProfilePolicy>
std::ostream & SphericalDiffuse<ProfilePolicy>::printObject(std::ostream & os) {
  Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, 0, ", ", "\n", "(", ")");
  os << "Green's function type: spherical diffuse" << std::endl;
  os << this->profile_ << std::endl;
  os << "Sphere center        = " << this->origin_.transpose().format(CleanFmt)
     << std::endl;
  os << "Angular momentum (Green's function)    = " << this->maxLGreen_ << std::endl;
  os << "Angular momentum (Coulomb coefficient) = " << this->maxLC_;
  return os;
}

template <typename ProfilePolicy>
void SphericalDiffuse<ProfilePolicy>::initSphericalDiffuse() {
  using namespace detail;

  // Parameters for the numerical solution of the radial differential equation
  double r_0_ = 0.1; /*! Lower bound of the integration interval */
  double r_infinity_ = this->profile_.upperLimit() +
                       30.0; /*! Upper bound of the integration interval */

  double y_0_ = std::log(r_0_);
  double y_infinity_ = std::log(r_infinity_);
  double relative_width = this->profile_.relativeWidth();
  double observer_step_ =
      1.0e-2 * relative_width; /*! Time step between observer calls */

  IntegratorParameters params_(y_0_, y_infinity_, observer_step_);
  ProfileEvaluator eval_ =
      pcm::bind(&ProfilePolicy::operator(), this->profile_, pcm::_1);

  zetaC_ = RadialFunction<StateType, LnTransformedRadial, Zeta>(
      maxLC_, y_0_, y_infinity_, eval_, params_);
  omegaC_ = RadialFunction<StateType, LnTransformedRadial, Omega>(
      maxLC_, y_0_, y_infinity_, eval_, params_);
  zeta_.reserve(maxLGreen_ + 1);
  omega_.reserve(maxLGreen_ + 1);
  for (int L = 0; L <= maxLGreen_; ++L) {
    RadialFunction<StateType, LnTransformedRadial, Zeta> tmp_zeta_(
        L, y_0_, y_infinity_, eval_, params_);
    zeta_.push_back(tmp_zeta_);
  }
  for (int L = 0; L <= maxLGreen_; ++L) {
    RadialFunction<StateType, LnTransformedRadial, Omega> tmp_omega_(
        L, y_0_, y_infinity_, eval_, params_);
    omega_.push_back(tmp_omega_);
  }
}

template <typename ProfilePolicy>
double SphericalDiffuse<ProfilePolicy>::imagePotentialComponent_impl(
    int L,
    const Eigen::Vector3d & sp,
    const Eigen::Vector3d & pp,
    double Cr12) const {
  Eigen::Vector3d sp_shift = sp + this->origin_;
  Eigen::Vector3d pp_shift = pp + this->origin_;
  double r1 = sp_shift.norm();
  double r2 = pp_shift.norm();
  double cos_gamma = sp_shift.dot(pp_shift) / (r1 * r2);
  // Evaluate Legendre polynomial of order L
  // First of all clean-up cos_gamma, Legendre polynomials
  // are only defined for -1 <= x <= 1
  if (utils::numericalZero(cos_gamma - 1))
    cos_gamma = 1.0;
  if (utils::numericalZero(cos_gamma + 1))
    cos_gamma = -1.0;
  double pl_x = Legendre::Pn<double>(L, cos_gamma);

  /* Zeta and Omega are now on a logarithmic scale We need to pass the
         arguments in the correct scale and then use the chain rule to get
         the proper derivative */
  double y1 = std::log(r1);
  double y2 = std::log(r2);

  /* Sample zeta_[L] */
  double zeta1 = 0.0, zeta2 = 0.0, d_zeta2 = 0.0;
  /* Value of zeta_[L] at point with index 1 */
  pcm::tie(zeta1, pcm::ignore) = zeta_[L](y1);
  /* Value of zeta_[L] and its first derivative at point with index 2 */
  pcm::tie(zeta2, d_zeta2) = zeta_[L](y2);

  /* Sample omega_[L] */
  double omega1 = 0.0, omega2 = 0.0, d_omega2 = 0.0;
  /* Value of omega_[L] at point with index 1 */
  pcm::tie(omega1, pcm::ignore) = omega_[L](y1);
  /* Value of omega_[L] and its first derivative at point with index 2 */
  pcm::tie(omega2, d_omega2) = omega_[L](y2);

  double eps_r2 = 0.0;
  pcm::tie(eps_r2, pcm::ignore) = this->profile_(pp_shift.norm());

  /* Evaluation of the Wronskian and the denominator */
  double denominator = (d_zeta2 - d_omega2) * r2 * eps_r2;

  double gr12 = 0.0;
  if (r1 < r2) {
    gr12 = std::exp(zeta1 - zeta2) * (2 * L + 1) / denominator;
    double f_L = r1 / r2;
    for (int i = 1; i < L; ++i) {
      f_L *= r1 / r2;
    }
    gr12 = (gr12 - f_L / (r2 * Cr12)) * pl_x;
  } else {
    gr12 = std::exp(omega1 - omega2) * (2 * L + 1) / denominator;
    double f_L = r2 / r1;
    for (int i = 1; i < L; ++i) {
      f_L *= r2 / r1;
    }
    gr12 = (gr12 - f_L / (r1 * Cr12)) * pl_x;
  }
  return gr12;
}

template <typename ProfilePolicy>
double SphericalDiffuse<ProfilePolicy>::coefficient_impl(
    const Eigen::Vector3d & sp,
    const Eigen::Vector3d & pp) const {
  double r1 = (sp + this->origin_).norm();
  double r2 = (pp + this->origin_).norm();

  double y1 = std::log(r1);
  double y2 = std::log(r2);

  /* Sample zetaC_ */
  double zeta1 = 0.0, zeta2 = 0.0, d_zeta2 = 0.0;
  /* Value of zetaC_ at point with index 1 */
  pcm::tie(zeta1, pcm::ignore) = zetaC_(y1);
  /* Value of zetaC_ and its first derivative at point with index 2 */
  pcm::tie(zeta2, d_zeta2) = zetaC_(y2);

  /* Sample omegaC_ */
  double omega1 = 0.0, omega2 = 0.0, d_omega2 = 0.0;
  /* Value of omegaC_ at point with index 1 */
  pcm::tie(omega1, pcm::ignore) = omegaC_(y1);
  /* Value of omegaC_ and its first derivative at point with index 2 */
  pcm::tie(omega2, d_omega2) = omegaC_(y2);

  double tmp = 0.0, coeff = 0.0;
  double eps_r2 = 0.0;
  pcm::tie(eps_r2, pcm::ignore) = this->profile_(r2);

  /* Evaluation of the Wronskian and the denominator */
  double denominator = (d_zeta2 - d_omega2) * r2 * eps_r2;

  if (r1 < r2) {
    double f_L = r1 / r2;
    for (int i = 1; i < maxLC_; ++i) {
      f_L *= r1 / r2;
    }
    tmp = std::exp(zeta1 - zeta2) * (2 * maxLC_ + 1) / denominator;
    coeff = f_L / (tmp * r2);
  } else {
    double f_L = r2 / r1;
    for (int i = 1; i < maxLC_; ++i) {
      f_L *= r2 / r1;
    }
    tmp = std::exp(omega1 - omega2) * (2 * maxLC_ + 1) / denominator;
    coeff = f_L / (tmp * r1);
  }

  return coeff;
}

using dielectric_profile::OneLayerTanh;
template class SphericalDiffuse<OneLayerTanh>;

using dielectric_profile::OneLayerErf;
template class SphericalDiffuse<OneLayerErf>;

using dielectric_profile::OneLayerLog;
template class SphericalDiffuse<OneLayerLog>;

} // namespace green
} // namespace pcm
