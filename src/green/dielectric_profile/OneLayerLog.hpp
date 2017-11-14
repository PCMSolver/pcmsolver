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

#include <cmath>
#include <iosfwd>

#include "Config.hpp"

#include <boost/math/special_functions/erf.hpp>

/*! \file OneLayerLog.hpp
 *  \class OneLayerLog
 *  \brief A dielectric profile based on the Harrison and Fosso-Tande work
 *  \author Luca Frediani
 *  \date 2017
 */

namespace pcm {
namespace dielectric_profile {
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
  /*! Returns value of dielectric profile at given point
   *  \param[in] point where to evaluate the profile
   */
	double value(double point) const {
		double epsLog = std::log(epsilon2/epsilon1);
		double val = (1.0 + boost::math::erf(point - center_) / width_) / 2.0;
		return std::eps1 * std::exp(epsLog * val); // epsilon(r)
	}
	/*! Returns value of derivative of dielectric profile at given point
	 *  \param[in] point where to evaluate the derivative
	 */
	double derivative(double point) const {
		double functionValue = value(point);
		double epsLog = std::log(epsilon2/epsilon1);
		double factor = epsLog / (width * std::sqrt(M_PI));
		double t = (point - center_) / width_;
		double val = std::exp(-std::pow(t, 2));
		return functionValue * factor * val; // first derivative of epsilon(r)
	}
  std::ostream & printObject(std::ostream & os) {
    os << "Profile functional form: erf" << std::endl;
    os << "Permittivity left/inside   = " << epsilon1_ << std::endl;
    os << "Permittivity right/outside = " << epsilon2_ << std::endl;
    os << "Profile width        = " << width_ << " AU" << std::endl;
    os << "Profile center       = " << center_ << " AU";
    return os;
  }

public:
  OneLayerErf() {}
  OneLayerErf(double e1, double e2, double w, double c)
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
  friend std::ostream & operator<<(std::ostream & os, OneLayerErf & th) {
    return th.printObject(os);
  }
};
} // namespace dielectric_profile
} // namespace pcm
