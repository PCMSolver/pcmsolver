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

#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

/*! \file Sphere.hpp */

namespace pcm {
namespace utils {
/*! \struct Sphere
 *  \brief POD describing a sphere.
 *  \author Roberto Di Remigio
 *  \date 2011, 2016
 */
struct Sphere {
  Sphere() {}
  Sphere(const Eigen::Vector3d & c, double r) : center(c), radius(r) {}
  ~Sphere() {}
  /// Scale sphere to other units
  void scale(double scaling) {
    center *= scaling;
    radius *= scaling;
  }
  Eigen::Vector3d center;
  double radius;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See
                                     http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
                                     */
      friend std::ostream &
      operator<<(std::ostream & os, Sphere & sph) {
    os << "Sphere radius " << sph.radius << std::endl;
    os << "Sphere center " << sph.center.transpose();

    return os;
  }
};
} // namespace utils

/*! \fn inline void transfer_spheres(const std::vector<Sphere> & spheres,
 *Eigen::Matrix3Xd & sphereCenter, Eigen::VectorXd & sphereRadius)
 *  \brief Transfer info from std::vector<Sphere> to Eigen objects.
 *  \param[in] spheres list of spheres as std::vector<Sphere>
 *  \param[out] sphereCenter sphere centers as Eigen::Matrix3Xd (xyz * nSpheres)
 *  \param[out] sphereRadius sphere radii as Eigen::VectorXd
 *
 *  This is used in the ICavity.hpp constructor
 */
void transfer_spheres(const std::vector<utils::Sphere> & spheres,
                      Eigen::Matrix3Xd & sphereCenter,
                      Eigen::VectorXd & sphereRadius);
} // namespace pcm
