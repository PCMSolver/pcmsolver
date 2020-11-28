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

#include <cmath>
#include <iosfwd>

#include "Config.hpp"

/*! \file OneLayerLog.hpp */

namespace pcm {
namespace dielectric_profile {
/*!  \class OneLayerLog
 *  \brief A dielectric profile based on the Harrison and Fosso-Tande work
 *  \cite Fosso-Tande2013
 *  \author Luca Frediani
 *  \date 2017
 */
class OneLayerLog {
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
   */
  double value(double point) const {
    double epsLog = std::log(epsilon2_ / epsilon1_);
    double val = (1.0 + std::erf((point - center_) / width_)) / 2.0;
    double retval = epsilon1_ * std::exp(epsLog * val); // epsilon(r)
    return retval;
  }
  /*! Returns value of derivative of dielectric profile at given point
   *  \param[in] point where to evaluate the derivative
   */
  double derivative(double point) const {
    double functionValue = value(point);
    double epsLog = std::log(epsilon2_ / epsilon1_);
    double factor = epsLog / (width_ * std::sqrt(M_PI));
    double t = (point - center_) / width_;
    double val = std::exp(-std::pow(t, 2));
    return functionValue * factor * val; // first derivative of epsilon(r)
  }
  std::ostream & printObject(std::ostream & os) {
    os << "Profile functional form: log" << std::endl;
    os << "Permittivity left/inside   = " << epsilon1_ << std::endl;
    os << "Permittivity right/outside = " << epsilon2_ << std::endl;
    os << "Profile width        = " << width_ << " AU" << std::endl;
    os << "Profile center       = " << center_ << " AU";
    return os;
  }

public:
  OneLayerLog() {}
  OneLayerLog(double e1, double e2, double w, double c)
      : epsilon1_(e1),
        epsilon2_(e2),
        width_(w / 6.0),
        center_(c),
        domain_(std::make_pair(0.0, center_ + 12.0 * width_)) {}
  /*! Returns a tuple holding the permittivity and its derivative
   *  \param[in]   r evaluation point
   */
  std::tuple<double, double> operator()(const double r) const {
    return std::make_tuple(value(r), derivative(r));
  }
  double upperLimit() const { return domain_.second; }
  double relativeWidth() const {
    return width_ / std::abs(domain_.second - domain_.first);
  }
  friend std::ostream & operator<<(std::ostream & os, OneLayerLog & th) {
    return th.printObject(os);
  }
};
} // namespace dielectric_profile
} // namespace pcm
