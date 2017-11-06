/**
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

#include "Config.hpp"

#include <Eigen/Core>

#include "utils/MathUtils.hpp"

/*! \file Anisotropic.hpp
 *  \class Anisotropic
 *  \brief describes a medium with anisotropy, i.e. liquid crystal
 *  \author Roberto Di Remigio
 *  \date 2014
 */

namespace pcm {
namespace dielectric_profile {
class Anisotropic __final {
private:
  /// Diagonal of the permittivity tensor in the lab-fixed frame
  Eigen::Vector3d epsilonLab_;
  /// Euler angles (in degrees) relating molecule-fixed and lab-fixed frames
  Eigen::Vector3d eulerAngles_;
  /// Permittivity tensor in molecule-fixed frame
  Eigen::Matrix3d epsilon_;
  /// Inverse of the permittivity tensor in molecule-fixed frame
  Eigen::Matrix3d epsilonInv_;
  /// molecule-fixed to lab-fixed frames rotation matrix
  Eigen::Matrix3d R_;
  /// Determinant of the permittivity tensor
  double detEps_;
  /*! Initializes some internals: molecule-fixed to lab-fixed frame rotation matrix,
   * permittivity tensor in molecule-fixed frame and its inverse
   */
  void build() {
    // 1. construct rotation matrix from Euler angles
    utils::eulerRotation(R_, eulerAngles_);
    // 2. Apply the rotation matrix: epsilon_ = R_^t * epsilonLab_ * R_
    epsilon_ = R_.transpose() * epsilonLab_.asDiagonal() * R_;
    // 3. Obtain epsilonInv_ = R_ * epsilonLab_^-1 * R_^t
    Eigen::Vector3d scratch;
    scratch << (1.0 / epsilonLab_(0)), (1.0 / epsilonLab_(1)),
        (1.0 / epsilonLab_(2));
    epsilonInv_ = R_ * scratch.asDiagonal() * R_.transpose();
    // 4. As a __final step, calculate the determinant
    detEps_ = epsilonLab_(0) * epsilonLab_(1) * epsilonLab_(2);
  }

public:
  Anisotropic()
      : epsilonLab_(Eigen::Vector3d::Ones()), eulerAngles_(Eigen::Vector3d::Zero()) {
    this->build();
  }
  /*!
   * \param[in] eigen_eps eigenvalues of the permittivity tensors
   * \param[in] euler_ang Euler angles in degrees
   */
  Anisotropic(const Eigen::Vector3d & eigen_eps, const Eigen::Vector3d & euler_ang)
      : epsilonLab_(eigen_eps), eulerAngles_(euler_ang) {
    this->build();
  }
  const Eigen::Matrix3d & epsilon() const { return epsilon_; }
  const Eigen::Matrix3d & epsilonInv() const { return epsilonInv_; }
  double detEps() const { return detEps_; }
  friend std::ostream & operator<<(std::ostream & os, Anisotropic & arg) {
    os << "Permittivity tensor diagonal (lab frame)   = "
       << arg.epsilonLab_.transpose() << std::endl;
    os << "Euler angles (molecule-to-lab frame)       = "
       << arg.eulerAngles_.transpose() << std::endl;
    os << "Permittivity tensor (molecule-fixed frame) =\n" << arg.epsilon_;
    return os;
  }
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW /* See
                                     http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
                                     */
};
} // namespace dielectric_profile
} // namespace pcm
