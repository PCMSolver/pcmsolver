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

#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

/*! @struct GreenData
 *  @brief Contains all data defined from user input in the green section.
 */

namespace pcm {
struct GreenData {
  /*! The way the derivatives of the Green's function are evaluated */
  int howDerivative;
  /*! Dielectric profile type */
  int howProfile;
  /*! The permittivity */
  double epsilon;
  /*! Inverse of the Debye length */
  double kappa;
  /*! Diagonal values of the permittivity tensor with respect to the lab frame */
  Eigen::Vector3d epsilonTensor;
  /*! Euler angles giving the rotation of the solvent orientation with respect to the
   * lab frame */
  Eigen::Vector3d eulerAngles;
  /*! Real part of the permittivity of a metal sphere */
  double epsilonReal;
  /*! Imaginary part of the permittivity of a metal sphere */
  double epsilonImaginary;
  /*! Coordinates of the metal sphere center */
  std::vector<double> NPspheres;
  /*! Radius of the the metal sphere */
  double NPradii;
  /*! Permittivity inside the interface */
  double epsilon1;
  /*! Permittivity outside the interface */
  double epsilon2;
  /*! Center of the diffuse layer */
  double center;
  /*! Width of the diffuse layer */
  double width;
  /*! Origin of the dielectric layer */
  Eigen::Vector3d origin;
  /*! Maximum angular momentum */
  int maxL;
  /*! Whether the structure was initialized with user input or not */
  bool empty;

  GreenData() { empty = true; }
  GreenData(int how_d,
            int how_p,
            double _epsilon = 1.0,
            double _kappa = 0.0,
            const Eigen::Vector3d & epstens = Eigen::Vector3d::Zero(),
            const Eigen::Vector3d & euler = Eigen::Vector3d::Zero(),
            double _epsReal = 0.0,
            double _epsImaginary = 0.0,
            const std::vector<double> & _sphere = std::vector<double>(),
            double _sphRadius = 0.0,
            double _e1 = 1.0,
            double _e2 = 1.0,
            double _c = 100.0,
            double _w = 5.0,
            const Eigen::Vector3d & _o = Eigen::Vector3d::Zero(),
            int l = 30)
      : howDerivative(how_d),
        howProfile(how_p),
        epsilon(_epsilon),
        kappa(_kappa),
        epsilonTensor(epstens),
        eulerAngles(euler),
        epsilonReal(_epsReal),
        epsilonImaginary(_epsImaginary),
        NPspheres(_sphere),
        NPradii(_sphRadius),
        epsilon1(_e1),
        epsilon2(_e2),
        center(_c),
        width(_w),
        origin(_o),
        maxL(l) {
    empty = false;
  }
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See
                                     http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
                                     */
};
} // namespace pcm
