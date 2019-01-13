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

#include "catch.hpp"

#include <cmath>
#include <iostream>

#include <Eigen/Core>

#include "AnalyticEvaluate.hpp"
#include "green/DerivativeTypes.hpp"
#include "green/UniformDielectric.hpp"

using namespace pcm;
using green::UniformDielectric;

TEST_CASE(
    "Evaluation of the uniform dielectric Green's function and its derivatives",
    "[green][green_uniform_dielectric]") {
  double epsilon = 60.0;
  Eigen::Vector3d source = Eigen::Vector3d::Random();
  Eigen::Vector3d sourceNormal = source + Eigen::Vector3d::Random();
  sourceNormal.normalize();
  Eigen::Vector3d probe = Eigen::Vector3d::Random();
  Eigen::Vector3d probeNormal = probe + Eigen::Vector3d::Random();
  probeNormal.normalize();
  Eigen::Array4d result =
      analyticUniformDielectric(epsilon, sourceNormal, source, probeNormal, probe);
  /*! \class UniformDielectric
   *  \test \b UniformDielectricTest_numerical tests the numerical evaluation of the
   * UniformDielectric Green's function against analytical result
   */
  SECTION("Numerical derivative") {
    UniformDielectric<Stencil> gf(epsilon);
    double value = result(0);
    double gf_value = gf.kernelS(source, probe);
    REQUIRE(value == Approx(gf_value));

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    REQUIRE(derProbe == Approx(gf_derProbe));

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    REQUIRE(derSource == Approx(gf_derSource));
  }

  /*! \class UniformDielectric
   *  \test \b UniformDielectricTest_directional_AD tests the automatic evaluation
   * (directional derivative only)
   *  of the UniformDielectric Green's function against analytical result
   */
  SECTION("Directional derivative via AD") {
    UniformDielectric<> gf(epsilon);
    double value = result(0);
    double gf_value = gf.kernelS(source, probe);
    REQUIRE(value == Approx(gf_value));

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    REQUIRE(derProbe == Approx(gf_derProbe));

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    REQUIRE(derSource == Approx(gf_derSource));
  }

  /*! \class UniformDielectric
   *  \test \b UniformDielectricTest_gradient_AD tests the automatic evaluation (full
   * gradient)
   *  of the UniformDielectric Green's function against analytical result
   */
  SECTION("Gradient via AD") {
    UniformDielectric<AD_gradient> gf(epsilon);
    double value = result(0);
    double gf_value = gf.kernelS(source, probe);
    REQUIRE(value == Approx(gf_value));

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    REQUIRE(derProbe == Approx(gf_derProbe));

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    REQUIRE(derSource == Approx(gf_derSource));
  }

  /*! \class UniformDielectric
   *  \test \b UniformDielectricTest_hessian_AD tests the automatic evaluation (full
   * hessian)
   *  of the UniformDielectric Green's function against analytical result
   */
  SECTION("Hessian via AD") {
    UniformDielectric<AD_hessian> gf(epsilon);
    double value = result(0);
    double gf_value = gf.kernelS(source, probe);
    REQUIRE(value == Approx(gf_value));

    double derProbe = result(1);
    double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
    REQUIRE(derProbe == Approx(gf_derProbe));

    double derSource = result(2);
    double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
    REQUIRE(derSource == Approx(gf_derSource));

    /*	double hessian = result(4);
        double gf_hessian = gf.hessian(sourceNormal, source, probeNormal, probe);
        REQUIRE(hessian == Approx(gf_hessian));
        */
  }
}
