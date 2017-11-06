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

/*! \file OneLayerTanh.hpp
 *  \class OneLayerTanh
 *  \brief A tanh dielectric profile as in \cite Frediani2004a
 *  \author Roberto Di Remigio
 *  \date 2014
 *  \note The parameter given from user input for width_ is divided by 6.0 in
 *  the constructor to keep consistency with \cite Frediani2004a
 */

namespace pcm {
namespace dielectric_profile {
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
  /*! Returns value of dielectric profile at given point
   *  \param[in] point where to evaluate the profile
   */
  double value(double point) const {
    double epsPlus = (epsilon1_ + epsilon2_) / 2.0;
    double epsMinus = (epsilon2_ - epsilon1_) / 2.0;
    double tanh_r = std::tanh((point - center_) / width_);
    return (epsPlus + epsMinus * tanh_r); // epsilon(r)
  }
  /*! Returns value of derivative of dielectric profile at given point
   *  \param[in] point where to evaluate the derivative
   */
  double derivative(double point) const {
    double factor = (epsilon2_ - epsilon1_) / (2.0 * width_);
    double tanh_r = std::tanh((point - center_) / width_);
    return (factor * (1 - std::pow(tanh_r, 2))); // first derivative of epsilon(r)
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
      : epsilon1_(e1), epsilon2_(e2), width_(w / 6.0), center_(c) {}
  /*! Returns a tuple holding the permittivity and its derivative
   *  \param[in]   r evaluation point
   */
  pcm::tuple<double, double> operator()(const double r) const {
    return pcm::make_tuple(value(r), derivative(r));
  }
  double epsilon1() const { return epsilon1_; }
  double epsilon2() const { return epsilon2_; }
  double width() const { return width_; }
  double center() const { return center_; }
  friend std::ostream & operator<<(std::ostream & os, OneLayerTanh & th) {
    return th.printObject(os);
  }
};
} // namespace dielectric_profile
} // namespace pcm
