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

#ifndef ANALYTICEVALUATE_HPP
#define ANALYTICEVALUATE_HPP

#include <cmath>

#include <Eigen/Core>

#include "utils/Legendre.hpp"
#include "utils/MathUtils.hpp"
#include "utils/Stencils.hpp"

/*! \brief Analytic evaluation of vacuum Green's function and its derivatives
 *
 *  \f[
 *  \begin{align}
 *     G(\vect{r},\vect{r}^\prime) &= \frac{1}{|\vect{r}-\vect{r}^\prime|} \\
 *     \pderiv{}{{\vect{n}_{\vect{r}^\prime}}}G(\vect{r},\vect{r}^\prime) &=
 *\frac{(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}^\prime}}{|\vect{r}-\vect{r}^\prime|^3} \\
 *     \pderiv{}{{\vect{n}_{\vect{r}}}}G(\vect{r},\vect{r}^\prime) &=
 *-\frac{(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}}}{|\vect{r}-\vect{r}^\prime|^3} \\
 *     \frac{\partial^2}{\partial{\vect{n}_{\vect{r}}}\partial{\vect{n}_{\vect{r}^\prime}}}G(\vect{r},\vect{r}^\prime)
 *&=
 *     \frac{\vect{n}_{\vect{r}}\cdot
 *\vect{n}_{\vect{r}^\prime}}{|\vect{r}-\vect{r}^\prime|^3}
 *     -3\frac{[(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}}][(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}^\prime}]}{|\vect{r}-\vect{r}^\prime|^5}
 *  \end{align}
 *  \f]
 */
inline Eigen::Array4d analyticVacuum(const Eigen::Vector3d & spNormal,
                                     const Eigen::Vector3d & sp,
                                     const Eigen::Vector3d & ppNormal,
                                     const Eigen::Vector3d & pp) {
  Eigen::Array4d result = Eigen::Array4d::Zero();
  double distance = (sp - pp).norm();
  double distance_3 = std::pow(distance, 3);
  double distance_5 = std::pow(distance, 5);

  // Value of the function
  result(0) = 1.0 / distance;
  // Value of the directional derivative wrt probe
  result(1) = (sp - pp).dot(ppNormal) / distance_3;
  // Directional derivative wrt source
  result(2) = -(sp - pp).dot(spNormal) / distance_3;
  // Value of the Hessian
  result(3) = spNormal.dot(ppNormal) / distance_3 -
              3 * ((sp - pp).dot(spNormal)) * ((sp - pp).dot(ppNormal)) / distance_5;

  return result;
}

/*! \brief Analytic evaluation of uniform dielectric Green's function and its
 *derivatives
 *
 *  \f[
 *  \begin{align}
 *     G(\vect{r},\vect{r}^\prime) &= \frac{1}{|\diel\vect{r}-\vect{r}^\prime|} \\
 *     \pderiv{}{{\vect{n}_{\vect{r}^\prime}}}G(\vect{r},\vect{r}^\prime) &=
 *\frac{(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}^\prime}}{\diel|\vect{r}-\vect{r}^\prime|^3} \\
 *     \pderiv{}{{\vect{n}_{\vect{r}}}}G(\vect{r},\vect{r}^\prime) &=
 *-\frac{(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}}}{\diel|\vect{r}-\vect{r}^\prime|^3} \\
 *     \frac{\partial^2}{\partial{\vect{n}_{\vect{r}}}\partial{\vect{n}_{\vect{r}^\prime}}}G(\vect{r},\vect{r}^\prime)
 *&=
 *     \frac{\vect{n}_{\vect{r}}\cdot
 *\vect{n}_{\vect{r}^\prime}}{\diel|\vect{r}-\vect{r}^\prime|^3}
 *     -3\frac{[(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}}][(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}^\prime}]}{\diel|\vect{r}-\vect{r}^\prime|^5}
 *  \end{align}
 *  \f]
 */
inline Eigen::Array4d analyticUniformDielectric(double eps,
                                                const Eigen::Vector3d & spNormal,
                                                const Eigen::Vector3d & sp,
                                                const Eigen::Vector3d & ppNormal,
                                                const Eigen::Vector3d & pp) {
  Eigen::Array4d result = Eigen::Array4d::Zero();
  double distance = (sp - pp).norm();
  double distance_3 = std::pow(distance, 3);
  double distance_5 = std::pow(distance, 5);

  // Value of the function
  result(0) = 1.0 / (eps * distance);
  // Value of the directional derivative wrt probe
  result(1) = (sp - pp).dot(ppNormal) / (eps * distance_3);
  // Directional derivative wrt source
  result(2) = -(sp - pp).dot(spNormal) / (eps * distance_3);
  // Value of the Hessian
  result(3) =
      spNormal.dot(ppNormal) / (eps * distance_3) -
      3 * ((sp - pp).dot(spNormal)) * ((sp - pp).dot(ppNormal)) / (eps * distance_5);

  return result;
}

/*! \brief Analytic evaluation of ionic liquid Green's function and its derivatives
 *
 *  \f[
 *  \begin{align}
 *  \phantom{G(\vect{r},\vect{r}^\prime)}
 *  &\begin{aligned}
 *    G(\vect{r},\vect{r}^\prime) =
 *\frac{\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|}
 *  \end{aligned}\\
 *  &\begin{aligned}
 *    \pderiv{}{{\vect{n}_{\vect{r}^\prime}}}G(\vect{r},\vect{r}^\prime) =
 *    \frac{(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}^\prime}\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^3}
 *    +\kappa\frac{(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}^\prime}\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^2}
 *   \end{aligned}\\
 *   &\begin{aligned}
 *     \pderiv{}{{\vect{n}_{\vect{r}}}}G(\vect{r},\vect{r}^\prime) =
 *     -\frac{(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}}\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^3}
 *     -\kappa\frac{(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}}\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^2}
 *    \end{aligned}\\
 *   &\begin{aligned}
 *     \frac{\partial^2}{\partial{\vect{n}_{\vect{r}}}\partial{\vect{n}_{\vect{r}^\prime}}}G(\vect{r},\vect{r}^\prime)
 *&=
 *     \frac{\vect{n}_{\vect{r}}\cdot
 *\vect{n}_{\vect{r}^\prime}\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^3}
 *     -\kappa\frac{[(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}}][(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}^\prime}]\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^4}\\
 *     &-3\frac{[(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}}][(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}^\prime}]\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^5}
 *     + \kappa\frac{\vect{n}_{\vect{r}}\cdot
 *\vect{n}_{\vect{r}^\prime}\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^2}
 *\\
 *     &-\kappa^2\frac{[(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}}][(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}^\prime}]\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^3}\\
 *     &-2\kappa\frac{[(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}}][(\vect{r}-\vect{r}^\prime)\cdot
 *\vect{n}_{\vect{r}^\prime}]\mathrm{e}^{-\kappa|\vect{r}-\vect{r}^\prime|}}{4\pi\diel|\vect{r}-\vect{r}^\prime|^4}
 *    \end{aligned}
 *  \end{align}
 *  \f]
 */
inline Eigen::Array4d analyticIonicLiquid(double eps,
                                          double k,
                                          const Eigen::Vector3d & spNormal,
                                          const Eigen::Vector3d & sp,
                                          const Eigen::Vector3d & ppNormal,
                                          const Eigen::Vector3d & pp) {
  Eigen::Array4d result = Eigen::Array4d::Zero();
  double distance = (sp - pp).norm();
  double distance_3 = std::pow(distance, 3);
  double distance_5 = std::pow(distance, 5);

  // Value of the function
  result(0) = std::exp(-k * distance) / (eps * distance);
  // Value of the directional derivative wrt probe
  result(1) = (sp - pp).dot(ppNormal) * (1 + k * distance) *
              std::exp(-k * distance) / (eps * distance_3);
  // Directional derivative wrt source
  result(2) = -(sp - pp).dot(spNormal) * (1 + k * distance) *
              std::exp(-k * distance) / (eps * distance_3);
  // Value of the Hessian
  result(3) = spNormal.dot(ppNormal) * (1 + k * distance) * std::exp(-k * distance) /
                  (eps * distance_3) -
              std::pow(k, 2) * (sp - pp).dot(spNormal) * (sp - pp).dot(ppNormal) *
                  std::exp(-k * distance) / (eps * distance_3) -
              3 * (sp - pp).dot(spNormal) * (sp - pp).dot(ppNormal) *
                  (1 + k * distance) * std::exp(-k * distance) / (eps * distance_5);

  return result;
}

/*! \brief Analytic evaluation of anisotropic liquid Green's function and its
 *derivatives
 *
 *  \f[
 *  \begin{align}
 *   \phantom{G(\vect{r},\vect{r}^\prime)}
 *   &\begin{aligned}
 *     G(\vect{r},\vect{r}^\prime) =
 *\frac{1}{4\pi\sqrt{\det{\bm{\diel}}}\sqrt{(\vect{r}-\vect{r}^\prime)^t\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)}}
 *     \end{aligned}\\
 *   &\begin{aligned}
 *     \pderiv{}{{\vect{n}_{\vect{r}^\prime}}}G(\vect{r},\vect{r}^\prime) =
 *\frac{[\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)]\cdot
 *\vect{n}_{\vect{r}^\prime}}{4\pi\sqrt{\det{\bm{\diel}}}
 *     [(\vect{r}-\vect{r}^\prime)^t\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)]^{\frac{3}{2}}}
 *     \end{aligned}\\
 *   &\begin{aligned}
 *     \pderiv{}{{\vect{n}_{\vect{r}}}}G(\vect{r},\vect{r}^\prime) =
 *-\frac{[\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)]\cdot
 *\vect{n}_{\vect{r}}}{4\pi\sqrt{\det{\bm{\diel}}}
 *     [(\vect{r}-\vect{r}^\prime)^t\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)]^{\frac{3}{2}}}
 *     \end{aligned}\\
 *   &\begin{aligned}
 *     \frac{\partial^2}{\partial{\vect{n}_{\vect{r}}}\partial{\vect{n}_{\vect{r}^\prime}}}G(\vect{r},\vect{r}^\prime)
 *&=
 *     \frac{\vect{n}_{\vect{r}}\cdot[\bm{\diel}^{-1}
 *\vect{n}_{\vect{r}^\prime}]}{4\pi\sqrt{\det{\bm{\diel}}}
 *     [(\vect{r}-\vect{r}^\prime)^t\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)]^{\frac{3}{2}}}
 *\\
 *     &-3\frac{\lbrace[\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)]\cdot
 *\vect{n}_{\vect{r}}\rbrace\lbrace[\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)]\cdot
 *\vect{n}_{\vect{r}^\prime}\rbrace}{4\pi\sqrt{\det{\bm{\diel}}}
 *     [(\vect{r}-\vect{r}^\prime)^t\bm{\diel}^{-1}(\vect{r}-\vect{r}^\prime)]^{\frac{5}{2}}}
 *    \end{aligned}
 *   \end{align}
 *   \f]
 */
inline Eigen::Array4d analyticAnisotropicLiquid(const Eigen::Vector3d & epsilon,
                                                const Eigen::Vector3d & euler,
                                                const Eigen::Vector3d & spNormal,
                                                const Eigen::Vector3d & sp,
                                                const Eigen::Vector3d & ppNormal,
                                                const Eigen::Vector3d & pp) {
  Eigen::Array4d result = Eigen::Array4d::Zero();
  Eigen::Matrix3d epsilonInv, R;
  pcm::utils::eulerRotation(R, euler);
  Eigen::Vector3d scratch;
  scratch << (1.0 / epsilon(0)), (1.0 / epsilon(1)), (1.0 / epsilon(2));
  epsilonInv = R * scratch.asDiagonal() * R.transpose();
  double detEps = epsilon(0) * epsilon(1) * epsilon(2);
  Eigen::Vector3d right = epsilonInv * (sp - pp);
  Eigen::Vector3d left = (sp - pp).transpose();
  double distance = std::sqrt(left.dot(right));
  double distance_3 = std::pow(distance, 3);
  double distance_5 = std::pow(distance, 5);

  // Value of the function
  result(0) = 1.0 / (std::sqrt(detEps) * distance);
  // Value of the directional derivative wrt probe
  result(1) = right.dot(ppNormal) / (std::sqrt(detEps) * distance_3);
  // Directional derivative wrt source
  result(2) = -right.dot(spNormal) / (std::sqrt(detEps) * distance_3);
  // Value of the Hessian
  Eigen::Vector3d eps_ppNormal = epsilonInv * ppNormal;
  result(3) = spNormal.dot(eps_ppNormal) / (std::sqrt(detEps) * distance_3) -
              3 * (right.dot(spNormal)) * (right.dot(ppNormal)) /
                  (std::sqrt(detEps) * distance_5);

  return result;
}

inline double imagePotential(double eps,
                             double epsSolv,
                             double radius,
                             const Eigen::Vector3d & origin,
                             int maxL,
                             const Eigen::Vector3d & sp,
                             const Eigen::Vector3d & pp) {
  Eigen::Vector3d sp_origin = sp - origin;
  double sp_origin_norm = sp_origin.norm();
  Eigen::Vector3d pp_origin = pp - origin;
  double pp_origin_norm = pp_origin.norm();
  double cos_gamma =
      sp_origin.dot(pp_origin) / (sp_origin.norm() * pp_origin.norm());
  // Clean-up cos_gamma, Legendre polynomials are only defined for -1 <= x <= 1
  if (pcm::utils::numericalZero(cos_gamma - 1))
    cos_gamma = 1.0;
  if (pcm::utils::numericalZero(cos_gamma + 1))
    cos_gamma = -1.0;
  // Image charge position
  Eigen::Vector3d r_img = origin + std::pow(radius / pp_origin_norm, 2) * pp_origin;
  double sp_image_norm = (sp - r_img).norm();
  // Image charge
  double q_img = radius / pp_origin_norm;
  // Permittivity factor
  double factor = (eps - epsSolv) / (eps + epsSolv);

  // Image Green's function
  double G_img = factor * (q_img / sp_image_norm - q_img / sp_origin_norm);
  // Image Green's function
  // Accumulate Legendre polynomial expansion of image potential
  double f_0 = radius / (sp_origin_norm * pp_origin_norm);
  double f_l = f_0;
  for (int l = 1; l <= maxL; ++l) {
    f_l = f_l * radius * f_0;
    double C_0_l = (eps - epsSolv) * l / ((eps + epsSolv) * l + epsSolv);
    double pl_x = Legendre::Pn<double>(l, cos_gamma);
    G_img += f_l * (C_0_l - factor) * pl_x;
  }

  return G_img / epsSolv;
}

inline double derivativeImagePotential(double eps,
                                       double epsSolv,
                                       double radius,
                                       const Eigen::Vector3d & origin,
                                       const Eigen::Vector3d & sp,
                                       const Eigen::Vector3d & ppNormal,
                                       const Eigen::Vector3d & pp) {
  Eigen::Vector3d sp_origin = sp - origin;
  double sp_origin_norm = sp_origin.norm();
  Eigen::Vector3d pp_origin = pp - origin;
  double pp_origin_norm = pp_origin.norm();
  double cos_gamma = sp_origin.dot(pp_origin) / (sp_origin_norm * pp_origin_norm);
  // Clean-up cos_gamma, Legendre polynomials are only defined for -1 <= x <= 1
  if (pcm::utils::numericalZero(cos_gamma - 1))
    cos_gamma = 1.0;
  if (pcm::utils::numericalZero(cos_gamma + 1))
    cos_gamma = -1.0;
  double pp_origin_norm_3 = std::pow(pp_origin_norm, 3);

  Eigen::Vector3d r_img = origin + std::pow(radius / pp_origin_norm, 2) * pp_origin;
  Eigen::Vector3d pp_image = pp - r_img;
  double pp_image_norm_3 = std::pow(pp_image.norm(), 3);
  double factor = (eps - epsSolv) / (eps + epsSolv);

  double der_G_img = factor * (radius / sp_origin_norm) *
                     (pp_origin.dot(ppNormal) / pp_origin_norm_3 -
                      pp_image.dot(ppNormal) / pp_image_norm_3);
  /*
  // Accumulate Legendre polynomial expansion of image potential
  double f_0 = radius / (sp_origin_norm * pp_origin_norm);
  double f_l = f_0;
  double pp_origin_norm_l_3 = pp_origin_norm_3; // To accumulate
  (pp-origin).norm()^(l+3)
  double pp_origin_norm_l_1 = pp_origin_norm;   // To accumulate
  (pp-origin).norm()^(l+1)
  for (int l = 1; l <= 200; ++l) {
      f_l = f_l * radius * f_0;
      pp_origin_norm_l_3 *= pp_origin_norm;
      pp_origin_norm_l_1 *= pp_origin_norm;

      double pl_x = Legendre::Pn<double>(l, cos_gamma); // P_l(cos_gamma)
      double pl_1_x = Legendre::Pn<double>(l+1, cos_gamma); // P_(l+1)(cos_gamma)
      double cos_denom = std::pow(cos_gamma, 2) - 1;

      double tmp_a = ((l+1) * pl_x * pp_origin.dot(ppNormal)) / pp_origin_norm_l_3;
      double tmp_b = ((l+1) * (pl_1_x - cos_gamma * pl_x) ) / (pp_origin_norm_l_1 *
  cos_denom);
      double tmp_c = sp_origin.dot(ppNormal) / (sp_origin_norm * pp_origin_norm);
      double tmp_d = (pp_origin.dot(sp_origin) * pp_origin.dot(ppNormal)) /
  (sp_origin_norm * pp_origin_norm_3);
      double tmp_e = tmp_b * (tmp_c - tmp_d);

      double C_0_l = (eps - epsSolv) * l / ((eps + epsSolv) * l + epsSolv);
      double C_l = C_0_l - factor;

      der_G_img += f_l * C_l * (-tmp_a + tmp_e);
  }
  */

  return der_G_img / epsSolv;
}

/*! \brief Analytic evaluation of spherical sharp Green's function and its
 * derivatives
 *  Derivation details in J. Chem. Phys. 139, 0224105 (2013)
 *  TODO Now using a 7-point stencil for the derivatives. Should double-check formula
 * given in paper and code that up.
 */
inline Eigen::Array4d analyticSphericalSharp(double eps,
                                             double epsSolv,
                                             double radius,
                                             const Eigen::Vector3d & origin,
                                             int maxL,
                                             const Eigen::Vector3d & spNormal,
                                             const Eigen::Vector3d & sp,
                                             const Eigen::Vector3d & ppNormal,
                                             const Eigen::Vector3d & pp) {
  Eigen::Array4d result = Eigen::Array4d::Zero();
  double distance = (sp - pp).norm();
  double distance_3 = std::pow(distance, 3);

  double G_img = imagePotential(eps, epsSolv, radius, origin, maxL, sp, pp);
  // Value of the function
  result(0) = 1.0 / (epsSolv * distance) - G_img;

  double d_probe_G_img = sevenPointStencil(
      pcm::bind(
          imagePotential, eps, epsSolv, radius, origin, maxL, pcm::_1, pcm::_2),
      pp,
      sp,
      ppNormal,
      1.0e-04);
  // Value of the directional derivative wrt probe
  result(1) = (sp - pp).dot(ppNormal) / (epsSolv * distance_3) - d_probe_G_img;

  double d_source_G_img = sevenPointStencil(
      pcm::bind(
          imagePotential, eps, epsSolv, radius, origin, maxL, pcm::_1, pcm::_2),
      sp,
      pp,
      spNormal,
      1.0e-04);
  // Directional derivative wrt source
  result(2) = -(sp - pp).dot(spNormal) / (epsSolv * distance_3) - d_source_G_img;
  // Value of the Hessian
  result(3) = 0.0;

  return result;
}

#endif // ANALYTICEVALUATE_HPP
