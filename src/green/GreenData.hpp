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

#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

namespace pcm {
/*! @struct GreenData
 *  @brief Contains all data defined from user input in the green section.
 */
struct GreenData {
  /*! The Green's function type inside the cavity.
   * It encodes the Green's function type, derivative calculation strategy and
   * dielectric profile: TYPE_DERIVATIVE_PROFILE
   */
  std::string greensFunctionType;
  /*! The permittivity */
  double epsilon;
  /*! Inverse of the Debye length */
  double kappa;
  /*! Diagonal values of the permittivity tensor with respect to the lab frame */
  Eigen::Vector3d epsilonTensor;
  /*! Euler angles giving the rotation of the solvent orientation with respect to the
   * lab frame */
  Eigen::Vector3d eulerAngles;
  /*! Permittivity inside the interface */
  double epsilon1;
  /*! Permittivity outside the interface */
  double epsilon2;
  /*! Center of the diffuse/sharp layer aka the radius of diffuse/sharp sphere */
  double center;
  /*! Width of the diffuse layer */
  double width;
  /*! Origin of the dielectric diffuse/sharp layer aka the center of the
   * diffuse/sharp sphere.
   */
  Eigen::Vector3d origin;
  /*! Maximum angular momentum in the spherical diffuse/sharp Green's function
   * summation.
   */
  int maxL;

  GreenData(const std::string & type,
            double eps = 1.0,
            double k = 0.0,
            const Eigen::Vector3d & epstens = Eigen::Vector3d::Zero(),
            const Eigen::Vector3d & euler = Eigen::Vector3d::Zero(),
            double e1 = 1.0,
            double e2 = 1.0,
            double c = 100.0,
            double w = 5.0,
            const Eigen::Vector3d & o = Eigen::Vector3d::Zero(),
            int l = 50)
      : greensFunctionType(type),
        epsilon(eps),
        kappa(k),
        epsilonTensor(epstens),
        eulerAngles(euler),
        epsilon1(e1),
        epsilon2(e2),
        center(c),
        width(w),
        origin(o),
        maxL(l) {}
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
} // namespace pcm
