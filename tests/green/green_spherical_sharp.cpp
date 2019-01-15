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
#include "green/SphericalSharp.hpp"

using namespace pcm;
using green::SphericalSharp;

TEST_CASE("Evaluation of the spherical sharp Green's function and its derivatives",
          "[green][green_spherical_sharp]") {
  double epsSolv = 80.0;
  double epsNP = 40.0;
  double NPradius = 50.0;
  int maxL = 200;
  Eigen::Vector3d NPcenter = Eigen::Vector3d::Zero();
  Eigen::Vector3d source = Eigen::Vector3d::Random() * 100.0;
  Eigen::Vector3d sourceNormal = source + Eigen::Vector3d::Random();
  sourceNormal.normalize();
  Eigen::Vector3d probe = Eigen::Vector3d::Random() * 100.0;
  Eigen::Vector3d probeNormal = probe + Eigen::Vector3d::Random();
  probeNormal.normalize();
  Eigen::Array4d result = analyticSphericalSharp(epsNP,
                                                 epsSolv,
                                                 NPradius,
                                                 NPcenter,
                                                 maxL,
                                                 sourceNormal,
                                                 source,
                                                 probeNormal,
                                                 probe);

  /*! \class SphericalSharp
   *  \test \b SphericalSharpTest_numerical tests the numerical evaluation of the
   * SphericalSharp Green's function against analytical result
   */
  SECTION("Numerical derivative") {
    SphericalSharp<Stencil> gf(epsNP, epsSolv, NPradius, NPcenter, maxL);
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

  /*! \class SphericalSharp
   *  \test \b SphericalSharpTest_directional_AD tests the automatic evaluation
   * (directional derivative only)
   *  of the SphericalSharp Green's function against analytical result
   */
  SECTION("Directional derivative via AD") {
    SphericalSharp<> gf(epsNP, epsSolv, NPradius, NPcenter, maxL);
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

  /*! \class SphericalSharp
   *  \test \b SphericalSharpTest_gradient_AD tests the automatic evaluation (full
   * gradient)
   *  of the SphericalSharp Green's function against analytical result
   */
  SECTION("Gradient via AD") {
    SphericalSharp<AD_gradient> gf(epsNP, epsSolv, NPradius, NPcenter, maxL);
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

  /*! \class SphericalSharp
   *  \test \b SphericalSharpTest_hessian_AD tests the automatic evaluation (full
   * hessian)
   *  of the SphericalSharp Green's function against analytical result
   */
  SECTION("Hessian via AD") {
    SphericalSharp<AD_hessian> gf(epsNP, epsSolv, NPradius, NPcenter, maxL);
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
