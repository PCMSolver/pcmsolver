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

#include "Config.hpp"

#include <boost/any.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/variant.hpp>

#include "Anisotropic.hpp"
#include "MembraneTanh.hpp"
#include "Metal.hpp"
#include "OneLayerErf.hpp"
#include "OneLayerTanh.hpp"
#include "Sharp.hpp"
#include "Uniform.hpp"
#include "Yukawa.hpp"

namespace pcm {
namespace dielectric_profile {
/*! All possible profile types */
typedef boost::mpl::vector<Uniform,
                           Yukawa,
                           Anisotropic,
                           OneLayerTanh,
                           OneLayerErf,
                           MembraneTanh,
                           Metal,
                           Sharp>
    profile_types;

/*! One-layer diffuse profile types */
typedef boost::mpl::vector<OneLayerTanh, OneLayerErf> onelayer_diffuse_profile_types;

/*! Two-layer (aka membrane-like) diffuse profile types */
typedef boost::mpl::vector<MembraneTanh> membrane_diffuse_profile_types;

/*! An object of Permittivity type can have one of the types in the type sequence
 * profile_types */
typedef boost::make_variant_over<profile_types>::type Permittivity;

class isUniform : public boost::static_visitor<bool> {
public:
  bool operator()(const Uniform & /* arg */) const { return true; }
  bool operator()(const Yukawa & /* arg */) const { return false; }
  bool operator()(const Anisotropic & /* arg */) const { return false; }
  bool operator()(const OneLayerTanh & /* arg */) const { return false; }
  bool operator()(const OneLayerErf & /* arg */) const { return false; }
  bool operator()(const MembraneTanh & /* arg */) const { return false; }
  bool operator()(const Metal & /* arg */) const { return false; }
  bool operator()(const Sharp & /* arg */) const { return false; }
};

inline bool uniform(const Permittivity & arg) {
  return boost::apply_visitor(isUniform(), arg);
}

class epsilonValue : public boost::static_visitor<boost::any> {
public:
  double operator()(const Uniform & arg) const { return arg.epsilon; }
  pcm::tuple<double, double> operator()(const Yukawa & arg) const {
    return pcm::make_tuple(arg.epsilon, arg.kappa);
  }
  bool operator()(const Anisotropic & /* arg */) const { return false; }
  bool operator()(const OneLayerTanh & /* arg */) const { return false; }
  bool operator()(const OneLayerErf & /* arg */) const { return false; }
  bool operator()(const MembraneTanh & /* arg */) const { return false; }
  std::complex<double> operator()(const Metal & arg) const { return arg.epsilon; }
  double operator()(const Sharp & arg) const { return arg.epsilon; }
};

inline double epsilon(const Permittivity & arg) {
  return boost::any_cast<double>(boost::apply_visitor(epsilonValue(), arg));
}
} // namespace dielectric_profile
} // namespace pcm
