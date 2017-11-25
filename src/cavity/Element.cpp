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

#include "Element.hpp"

#include <cmath>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>

#include <boost/math/special_functions/sign.hpp>

#include "utils/Sphere.hpp"

namespace pcm {
namespace cavity {
void Element::spherical_polygon(Eigen::Vector3d & t_,
                                Eigen::Vector3d & b_,
                                std::vector<double> & theta,
                                std::vector<double> & phi,
                                std::vector<double> & phinumb,
                                std::vector<int> & numb) const {
  Eigen::Vector3d sph_center = sphere_.center;
  double radius = sphere_.radius;
  // Calculate the azimuthal and polar angles for the tessera vertices:
  // we use the normal, tangent and bitangent as a local reference frame
  for (int i = 0; i < nVertices_; ++i) {
    Eigen::Vector3d vertex_normal = vertices_.col(i) - sph_center;
    // The cosine of the polar angle is given as the dot product of the normal at the
    // vertex and the
    // normal at the tessera center: R\cos\theta
    double cos_theta = vertex_normal.dot(normal_) / radius;
    if (cos_theta >= 1.0)
      cos_theta = 1.0;
    if (cos_theta <= -1.0)
      cos_theta = -1.0;
    theta[i] = std::acos(cos_theta);
    // The cosine of the azimuthal angle is given as the dot product of the normal at
    // the vertex and the
    // tangent at the tessera center divided by the sine of the polar angle:
    // R\sin\theta\cos\phi
    double cos_phi = vertex_normal.dot(t_) / (radius * std::sin(theta[i]));
    if (cos_phi >= 1.0)
      cos_phi = 1.0;
    if (cos_phi <= -1.0)
      cos_phi = -1.0;
    phi[i] = std::acos(cos_phi);
    // The sine of the azimuthal angle is given as the dot product of the normal at
    // the vertex and the
    // bitangent at the tessera center divided by the sine of the polar angle:
    // R\sin\theta\sin\phi
    double sin_phi = vertex_normal.dot(b_) / (radius * std::sin(theta[i]));
    if (sin_phi <= 0.0)
      phi[i] = 2 * M_PI - phi[i];
  }
  for (int i = 1; i < nVertices_; ++i) {
    phi[i] = phi[i] - phi[0];
    if (phi[i] < 0.0)
      phi[i] = 2 * M_PI + phi[i];
  }
  // Rewrite tangent as linear combination of original tangent and bitangent
  // then recalculate bitangent so that it's orthogonal to the tangent
  t_ = t_ * std::cos(phi[0]) + b_ * std::sin(phi[0]);
  b_ = normal_.cross(t_);
  // Populate numb and phinumb arrays
  phi[0] = 0.0;
  numb[0] = 0;
  numb[1] = 1;
  phinumb[0] = phi[0];
  phinumb[1] = phi[1];
  for (int i = 2; i < nVertices_; ++i) { // This loop is 2-based
    for (int j = 1; j < i; ++j) {        // This loop is 1-based
      if (phi[i] < phinumb[j]) {
        for (int k = 0; k < (i - j); ++k) {
          numb[i - k] = numb[i - k - 1];
          phinumb[i - k] = phinumb[i - k - 1];
        }
        numb[j] = i;
        phinumb[j] = phi[i];
        goto jump;
      }
    }
    numb[i] = i;
    phinumb[i] = phi[i];
  jump:; // Do nothing...
  }
  numb[nVertices_] = numb[0];
  phinumb[nVertices_] = 2 * M_PI;
}

namespace detail {
void tangent_and_bitangent(const Eigen::Vector3d & n_,
                           Eigen::Vector3d & t_,
                           Eigen::Vector3d & b_) {
  double rmin = 0.99;
  double n0 = n_(0), n1 = n_(1), n2 = n_(2);
  if (std::abs(n0) <= rmin) {
    rmin = std::abs(n0);
    t_(0) = 0.0;
    t_(1) = -n2 / std::sqrt(1.0 - std::pow(n0, 2));
    t_(2) = n1 / std::sqrt(1.0 - std::pow(n0, 2));
  }
  if (std::abs(n1) <= rmin) {
    rmin = std::abs(n1);
    t_(0) = n2 / std::sqrt(1.0 - std::pow(n1, 2));
    t_(1) = 0.0;
    t_(2) = -n0 / std::sqrt(1.0 - std::pow(n1, 2));
  }
  if (std::abs(n2) <= rmin) {
    rmin = std::abs(n2);
    t_(0) = n1 / std::sqrt(1.0 - std::pow(n2, 2));
    t_(1) = -n0 / std::sqrt(1.0 - std::pow(n2, 2));
    t_(2) = 0.0;
  }
  b_ = n_.cross(t_);
  // Check that the calculated Frenet-Serret frame is left-handed (levogiro)
  // by checking that the determinant of the matrix whose columns are the normal,
  // tangent and bitangent vectors has determinant 1 (the system is orthonormal!)
  Eigen::Matrix3d M;
  M.col(0) = n_;
  M.col(1) = t_;
  M.col(2) = b_;
  if (boost::math::sign(M.determinant()) != 1) {
    PCMSOLVER_ERROR("Frenet-Serret local frame is not left-handed!");
  }
}
} // namespace detail
} // namespace cavity
} // namespace pcm
