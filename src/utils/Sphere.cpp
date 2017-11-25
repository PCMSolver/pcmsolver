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

#include "Sphere.hpp"

#include <ostream>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

namespace pcm {
using utils::Sphere;

std::ostream & operator<<(std::ostream & os, Sphere & sph) {
  os << "Sphere radius " << sph.radius << std::endl;
  os << "Sphere center\n" << sph.center.transpose();

  return os;
}

void transfer_spheres(const std::vector<Sphere> & spheres,
                      Eigen::Matrix3Xd & sphereCenter,
                      Eigen::VectorXd & sphereRadius) {
  size_t nSpheres = spheres.size();
  sphereCenter.resize(Eigen::NoChange, nSpheres);
  sphereRadius.resize(nSpheres);
  for (size_t i = 0; i < nSpheres; ++i) {
    sphereCenter.col(i) = spheres[i].center;
    sphereRadius(i) = spheres[i].radius;
  }
}
} // namespace pcm
