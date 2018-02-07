/*
 * PCMSolver, an API for the Polarizable Continuum Model
 * Copyright (C) 2018 Roberto Di Remigio, Luca Frediani and contributors.
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

/*! \file OneLayerTanh.hpp */

namespace pcm {
namespace dielectric_profile {
/*! \class OneLayerTanh
 *  \brief A tanh dielectric profile as in \cite Frediani2004a
 *  \author Roberto Di Remigio
 *  \date 2014
 *  \note The parameter given from user input for width_ is divided by 6.0 in
 *  the constructor to keep consistency with \cite Frediani2004a
 */
class OneLayerTanh {
private:
  /// Dielectric constant on the left of the interface
  double epsilon1_;
  /// Dielectric constant one the right of the interface
  double epsilon2_;
  /// Width of the transition layer
  double width_;
  /// Center of the transition layer
  double center_;
  /*! Domain of the permittivity function
   * This is formally \f$ [0, +\infty) \f$, for all practical purposes
   * the permittivity function is equal to the epsilon2_ already at 6.0 * width_
   * Thus the upper limit in the domain_ is initialized as center_ + 12.0 * width_
   */
  std::pair<double, double> domain_;
  /*! Returns value of dielectric profile at given point
   *  \param[in] point where to evaluate the profile
   *  \note We return epsilon2_ when the sampling point is outside the upper limit.
   */
  double value(double point) const {
    double retval = 0.0;
    if (point < domain_.first) {
      retval = epsilon1_;
    } else if (point > domain_.second) {
      retval = epsilon2_;
    } else {
      double epsPlus = (epsilon1_ + epsilon2_) / 2.0;
      double epsMinus = (epsilon2_ - epsilon1_) / 2.0;
      double tanh_r = std::tanh((point - center_) / width_);
      retval = epsPlus + epsMinus * tanh_r;
    }
    return retval;
  }
  /*! Returns value of derivative of dielectric profile at given point
   *  \param[in] point where to evaluate the derivative
   *  \note We return 0.0 (derivative of the constant value epsilon2_) when the
   * sampling point is outside the upper limit.
   */
  double derivative(double point) const {
    double retval = 0.0;
    if (point < domain_.first || point > domain_.second) {
      retval = 0.0;
    } else {
      double factor = (epsilon2_ - epsilon1_) / (2.0 * width_);
      double tanh_r = std::tanh((point - center_) / width_);
      retval = factor * (1 - std::pow(tanh_r, 2));
    }
    return retval;
  }
  std::ostream & printObject(std::ostream & os) {
    os << "Profile functional form: tanh" << std::endl;
    os << "Permittivity inside  = " << epsilon1_ << std::endl;
    os << "Permittivity outside = " << epsilon2_ << std::endl;
    os << "Profile width        = " << width_ << " AU" << std::endl;
    os << "Profile center       = " << center_ << " AU";
    return os;
  }

public:
  OneLayerTanh() {}
  OneLayerTanh(double e1, double e2, double w, double c)
      : epsilon1_(e1),
        epsilon2_(e2),
        width_(w / 6.0),
        center_(c),
        domain_(std::make_pair(0.0, center_ + 12.0 * width_)) {}
  /*! Returns a tuple holding the permittivity and its derivative
   *  \param[in]   r evaluation point
   */
  pcm::tuple<double, double> operator()(const double r) const {
    return pcm::make_tuple(value(r), derivative(r));
  }
  double upperLimit() const { return domain_.second; }
  friend std::ostream & operator<<(std::ostream & os, OneLayerTanh & th) {
    return th.printObject(os);
  }
};
} // namespace dielectric_profile
} // namespace pcm
