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

#include "catch.hpp"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include <Eigen/Core>
#include <boost/math/special_functions/erf.hpp>

#include "green/dielectric_profile/OneLayerErf.hpp"
#include "green/dielectric_profile/OneLayerTanh.hpp"

double tanh_value(double point, double e1, double e2, double w, double c);
double tanh_deriv(double point, double e1, double e2, double w, double c);

double erf_value(double point, double e1, double e2, double w, double c);
double erf_deriv(double point, double e1, double e2, double w, double c);

double distribution(double fMin, double fMax);

using namespace pcm::dielectric_profile;

SCENARIO("Diffuse permittivity single layers", "[dielectric_profile][one_layer]") {
  GIVEN("The parameters for a one-layer diffuse permittivity profile") {
    std::srand(std::time(0));
    double eps1 = 78.39;
    double eps2 = 1.0;
    double center = 50.0;
    double width = 10.0;

    /*! \class OneLayerTanh
     *  \test \b TanhOneLayerTest tests the evaluation of the one layer hyperbolic
     * tangent profile against the analytic value
     */
    WHEN("the sigmoidal profile is modelled by the hyperbolic tangent function") {
      OneLayerTanh diffuse(eps1, eps2, width, center);
      double value = 0.0, deriv = 0.0;
      double point = distribution(0.0, 100.0);
      pcm::tie(value, deriv) = diffuse(point);
      THEN("the value of the profile at a random point is") {
        double analytic = tanh_value(point, eps1, eps2, width, center);
        INFO(" The evaluation point is: " << point);
        REQUIRE(value == Approx(analytic));
      }
      AND_THEN("the value of the first derivative at a random point is") {
        double analyticDeriv = tanh_deriv(point, eps1, eps2, width, center);
        INFO(" The evaluation point is: " << point);
        REQUIRE(deriv == Approx(analyticDeriv));
      }
    }

    /*! \class OneLayerErf
     *  \test \b TanhOneLayerTest tests the evaluation of the one layer hyperbolic
     * tangent profile against the analytic value
     */
    WHEN("the sigmoidal profile is modelled by the error function") {
      OneLayerErf diffuse(eps1, eps2, width, center);
      double value = 0.0, deriv = 0.0;
      double point = distribution(0.0, 100.0);
      pcm::tie(value, deriv) = diffuse(point);
      THEN("the value of the profile at a random point is") {
        double analytic = erf_value(point, eps1, eps2, width, center);
        INFO(" The evaluation point is: " << point);
        REQUIRE(value == Approx(analytic));
      }
      AND_THEN("the value of the first derivative at a random point is") {
        double analyticDeriv = erf_deriv(point, eps1, eps2, width, center);
        INFO(" The evaluation point is: " << point);
        REQUIRE(deriv == Approx(analyticDeriv));
      }
    }
  }
}

double tanh_value(double point, double e1, double e2, double w, double c) {
  w /= 6.0;
  double epsPlus = (e1 + e2) / 2.0;
  double epsMinus = (e2 - e1) / 2.0;
  double tanh_r = std::tanh((point - c) / w);
  return (epsPlus + epsMinus * tanh_r);
}
double tanh_deriv(double point, double e1, double e2, double w, double c) {
  w /= 6.0;
  double factor = (e2 - e1) / (2.0 * w);
  double tanh_r = std::tanh((point - c) / w);
  return (factor * (1 - std::pow(tanh_r, 2)));
}

double erf_value(double point, double e1, double e2, double w, double c) {
  w /= 6.0;
  double epsPlus = (e1 + e2) / 2.0;
  double epsMinus = (e2 - e1) / 2.0;
  double val = boost::math::erf((point - c) / w);
  return (epsPlus + epsMinus * val);
}
double erf_deriv(double point, double e1, double e2, double w, double c) {
  w /= 6.0;
  double factor = (e2 - e1) / (w * std::sqrt(M_PI));
  double t = (point - c) / w;
  double val = std::exp(-std::pow(t, 2));
  return (factor * val);
}

double distribution(double fMin, double fMax) {
  double f = (double)std::rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}
