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

/*! \file MembraneTanh.hpp
 *  \class MembraneTanh
 *  \brief A dielectric profile mimicking a memmbrane as in Eq. 30 in \cite
 * Frediani2004a
 *  \author Roberto Di Remigio
 *  \date 2015
 */

namespace pcm {
namespace dielectric_profile {
class MembraneTanh {
private:
  double epsilon1_;
  double epsilon2_;
  double epsilon3_;
  double width12_;
  double width23_;
  double center12_;
  double center23_;
  double value(double point) const {
    double eps_13 = (epsilon1_ + epsilon3_) / 2.0;
    double eps_21 = (epsilon2_ - epsilon1_) / 2.0;
    double eps_23 = (epsilon2_ - epsilon3_) / 2.0;
    double tanh_r_12 = std::tanh((point - center12_) / width12_);
    double tanh_r_23 = std::tanh((point - center23_) / width23_);
    return (eps_13 + eps_21 * tanh_r_12 - eps_23 * tanh_r_23); // epsilon(r)
  }
  double derivative(double point) const {
    double factor_21 = (epsilon2_ - epsilon1_) / (2.0 * width12_);
    double factor_23 = (epsilon2_ - epsilon3_) / (2.0 * width23_);
    double tanh_r_12 = std::tanh((point - center12_) / width12_);
    double tanh_r_23 = std::tanh((point - center23_) / width23_);
    return (factor_21 * (1 - std::pow(tanh_r_12, 2)) -
            factor_23 *
                (1 - std::pow(tanh_r_23, 2))); // first derivative of epsilon(r)
  }

public:
  MembraneTanh() {}
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
        width12_(w12),
        width23_(w23),
        center12_(c12),
        center23_(c23) {}
  /*! The permittivity profile of the transition layer
   *  \param[out]  e the value of the dielectric constant at point r
   *  \param[out] de the value of the derivative of the dielectric constant
   *                 at point r
   *  \param[in]   r evaluation point
   */
  void operator()(double & e, double & de, const double r) const {
    e = value(r);
    de = derivative(r);
  }
  double epsilon1() const { return epsilon1_; }
  double epsilon2() const { return epsilon2_; }
  double epsilon3() const { return epsilon3_; }
  double width12() const { return width12_; }
  double width23() const { return width23_; }
  double center12() const { return center12_; }
  double center23() const { return center23_; }
};
} // namespace dielectric_profile
} // namespace pcm
