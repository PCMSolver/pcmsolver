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

#include <iosfwd>
#include <utility>

#include "Config.hpp"

/*! \file MembraneTanh.hpp */

namespace pcm {
namespace dielectric_profile {
/*! \class MembraneTanh
 *  \brief A dielectric profile mimicking a memmbrane as in Eq. 30 in
 *  \cite Frediani2004a
 *  \author Roberto Di Remigio
 *  \date 2015
 */
class MembraneTanh {
private:
  double epsilon1_;
  double epsilon2_;
  double epsilon3_;
  double width12_;
  double width23_;
  double center12_;
  double center23_;
  /*! Domain of the permittivity function
   * This is formally \f$ [0, +\infty) \f$, for all practical purposes
   * the permittivity function is equal to the epsilon3_ already at 6.0 * width23_
   * Thus the upper limit in the domain_ is initialized as center23_ + 12.0 *
   * width23_
   */
  std::pair<double, double> domain_;
  /*! Returns value of dielectric profile at given point
   *  \param[in] point where to evaluate the profile
   *  \note We return epsilon3_ when the sampling point is outside the upper limit.
   */
  double value(double point) const {
    double retval = 0.0;
    if (point < domain_.first) {
      retval = epsilon1_;
    } else if (point > domain_.second) {
      retval = epsilon3_;
    } else {
      double eps_13 = (epsilon1_ + epsilon3_) / 2.0;
      double eps_21 = (epsilon2_ - epsilon1_) / 2.0;
      double eps_23 = (epsilon2_ - epsilon3_) / 2.0;
      double tanh_r_12 = std::tanh((point - center12_) / width12_);
      double tanh_r_23 = std::tanh((point - center23_) / width23_);
      retval = eps_13 + eps_21 * tanh_r_12 - eps_23 * tanh_r_23;
    }
    return retval;
  }
  /*! Returns value of derivative of dielectric profile at given point
   *  \param[in] point where to evaluate the derivative
   *  \note We return 0.0 (derivative of the constant value epsilon3_) when the
   * sampling point is outside the upper limit.
   */
  double derivative(double point) const {
    double retval = 0.0;
    if (point < domain_.first || point > domain_.second) {
      retval = 0.0;
    } else {
      double factor_21 = (epsilon2_ - epsilon1_) / (2.0 * width12_);
      double factor_23 = (epsilon2_ - epsilon3_) / (2.0 * width23_);
      double tanh_r_12 = std::tanh((point - center12_) / width12_);
      double tanh_r_23 = std::tanh((point - center23_) / width23_);
      retval = factor_21 * (1 - std::pow(tanh_r_12, 2)) -
               factor_23 * (1 - std::pow(tanh_r_23, 2));
    }
    return retval;
  }
  std::ostream & printObject(std::ostream & os) {
    os << "Profile functional form: tanh" << std::endl;
    os << "Permittivity left-side  = " << epsilon1_ << std::endl;
    os << "Permittivity middle     = " << epsilon2_ << std::endl;
    os << "Permittivity right-side = " << epsilon3_ << std::endl;
    os << "Profile width, 1st layer = " << width12_ << " AU" << std::endl;
    os << "Profile width, 2nd layer = " << width23_ << " AU" << std::endl;
    os << "Profile center, 1st layer = " << center12_ << " AU" << std::endl;
    os << "Profile center, 2nd layer = " << center23_ << " AU";
    return os;
  }

public:
  MembraneTanh() {}
  /*! Constructor for a two-layer interface (membrane)
   * \param[in] e1 left-side dielectric constant
   * \param[in] e2 middle portion dielectric constant
   * \param[in] e3 right-side dielectric constant
   * \param[in] w12 width of the first interface layer
   * \param[in] w23 width of the second interface layer
   * \param[in] c12 center of the first diffuse layer
   * \param[in] c23 center of the second diffuse layer
   */
  MembraneTanh(double e1,
               double e2,
               double e3,
               double w12,
               double w23,
               double c12,
               double c23)
      : epsilon1_(e1),
        epsilon2_(e2),
        epsilon3_(e3),
        width12_(w12 / 6.0),
        width23_(w23 / 6.0),
        center12_(c12),
        center23_(c23),
        domain_(std::make_pair(0.0, center23_ + 12.0 * width23_)) {}
  /*! Returns a tuple holding the permittivity and its derivative
   *  \param[in]   r evaluation point
   */
  std::tuple<double, double> operator()(const double r) const {
    return std::make_tuple(value(r), derivative(r));
  }
  double upperLimit() const { return domain_.second; }
  friend std::ostream & operator<<(std::ostream & os, MembraneTanh & th) {
    return th.printObject(os);
  }
};
} // namespace dielectric_profile
} // namespace pcm
