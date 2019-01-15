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
#include "green/AnisotropicLiquid.hpp"
#include "green/DerivativeTypes.hpp"

using namespace pcm;
using green::AnisotropicLiquid;

SCENARIO("Evaluation of the anisotropic liquid Green's function and its derivatives",
         "[green][green_anisotropic_liquid]") {
  GIVEN("A liquid with an anisotropic permittivity tensor") {
    Eigen::Vector3d epsilon = (Eigen::Vector3d() << 2.0, 80.0, 15.0).finished();
    Eigen::Vector3d euler = (Eigen::Vector3d() << 6.0, 40.0, 15.0).finished();
    Eigen::Vector3d source = Eigen::Vector3d::Random();
    Eigen::Vector3d sourceNormal = source + Eigen::Vector3d::Random();
    sourceNormal.normalize();
    Eigen::Vector3d probe = Eigen::Vector3d::Random();
    Eigen::Vector3d probeNormal = probe + Eigen::Vector3d::Random();
    probeNormal.normalize();
    Eigen::Array4d result = analyticAnisotropicLiquid(
        epsilon, euler, sourceNormal, source, probeNormal, probe);

    /*! \class AnisotropicLiquid
     *  \test \b AnisotropicLiquidTest_numerical tests the numerical evaluation of
     * the AnisotropicLiquid Green's function against analytical result
     */
    WHEN("the derivatives are evaluated numerically") {
      AnisotropicLiquid<Stencil> gf(epsilon, euler);
      THEN("the value of the Green's function is") {
        double value = result(0);
        double gf_value = gf.kernelS(source, probe);
        REQUIRE(value == Approx(gf_value));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "probe point is") {
        double derProbe = result(1);
        double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "source point is") {
        double derSource = result(2);
        double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
        REQUIRE(derSource == Approx(gf_derSource));
      }
    }

    /*! \class AnisotropicLiquid
     *  \test \b AnisotropicLiquidTest_directional_AD tests the automatic evaluation
     * (directional derivative only)
     *  of the AnisotropicLiquid Green's function against analytical result
     */
    WHEN("the derivatives are evaluated via AD") {
      AnisotropicLiquid<> gf(epsilon, euler);
      THEN("the value of the Green's function is") {
        double value = result(0);
        double gf_value = gf.kernelS(source, probe);
        REQUIRE(value == Approx(gf_value));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "probe point is") {
        double derProbe = result(1);
        double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "source point is") {
        double derSource = result(2);
        double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
        REQUIRE(derSource == Approx(gf_derSource));
      }
    }

    /*! \class AnisotropicLiquid
     *  \test \b AnisotropicLiquidTest_gradient_AD tests the automatic evaluation
     * (full gradient)
     *  of the AnisotropicLiquid Green's function against analytical result
     */
    WHEN("the derivatives are evaluated via AD using the full gradient") {
      AnisotropicLiquid<AD_gradient> gf(epsilon, euler);
      THEN("the value of the Green's function is") {
        double value = result(0);
        double gf_value = gf.kernelS(source, probe);
        REQUIRE(value == Approx(gf_value));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "probe point is") {
        double derProbe = result(1);
        double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "source point is") {
        double derSource = result(2);
        double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
        REQUIRE(derSource == Approx(gf_derSource));
      }
    }

    /*! \class AnisotropicLiquid
     *  \test \b AnisotropicLiquidTest_hessian_AD tests the automatic evaluation
     * (full hessian)
     *  of the AnisotropicLiquid Green's function against analytical result
     */
    WHEN("the derivatives are evaluated via AD using the full hessian") {
      AnisotropicLiquid<AD_hessian> gf(epsilon, euler);
      THEN("the value of the Green's function is") {
        double value = result(0);
        double gf_value = gf.kernelS(source, probe);
        REQUIRE(value == Approx(gf_value));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "probe point is") {
        double derProbe = result(1);
        double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "source point is") {
        double derSource = result(2);
        double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
        REQUIRE(derSource == Approx(gf_derSource));
      }
      /*
         AND_THEN("the value of the Green's function hessian is")
         {
         double hessian = result(4);
         double gf_hessian = gf.hessian(sourceNormal, source, probeNormal, probe);
         REQUIRE(hessian == Approx(gf_hessian));
         }
         */
    }
  }

  // Define a uniform dielectric as an anisotropic dielectric and test the evaluation
  // of the
  // Green's function and derivatives against the analytic result for the uniform
  // dielectric
  GIVEN("A liquid with an isotropic permittivity tensor") {
    Eigen::Vector3d epsilon = (Eigen::Vector3d() << 80.0, 80.0, 80.0).finished();
    Eigen::Vector3d euler = (Eigen::Vector3d() << 0.0, 0.0, 0.0).finished();
    Eigen::Vector3d source = Eigen::Vector3d::Random();
    Eigen::Vector3d sourceNormal = source + Eigen::Vector3d::Random();
    sourceNormal.normalize();
    Eigen::Vector3d probe = Eigen::Vector3d::Random();
    Eigen::Vector3d probeNormal = probe + Eigen::Vector3d::Random();
    probeNormal.normalize();
    Eigen::Array4d result = analyticUniformDielectric(
        epsilon(0), sourceNormal, source, probeNormal, probe);

    /*! \class AnisotropicLiquid
     *  \test \b AnisotropicLiquidUniformTest_numerical tests the numerical
     * evaluation of the AnisotropicLiquid Green's function against analytical result
     * for a uniform dielectric
     */
    WHEN("the derivatives are evaluated numerically") {
      AnisotropicLiquid<Stencil> gf(epsilon, euler);
      THEN("the value of the Green's function is") {
        double value = result(0);
        double gf_value = gf.kernelS(source, probe);
        REQUIRE(value == Approx(gf_value));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "probe point is") {
        double derProbe = result(1);
        double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "source point is") {
        double derSource = result(2);
        double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
        REQUIRE(derSource == Approx(gf_derSource));
      }
    }

    /*! \class AnisotropicLiquid
     *  \test \b AnisotropicLiquidUniformTest_directional_AD tests the automatic
     * evaluation (directional derivative only)
     *  of the AnisotropicLiquid Green's function against analytical result for a
     * uniform dielectric
     */
    WHEN("the derivatives are evaluated via AD") {
      AnisotropicLiquid<AD_directional> gf(epsilon, euler);
      THEN("the value of the Green's function is") {
        double value = result(0);
        double gf_value = gf.kernelS(source, probe);
        REQUIRE(value == Approx(gf_value));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "probe point is") {
        double derProbe = result(1);
        double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "source point is") {
        double derSource = result(2);
        double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
        REQUIRE(derSource == Approx(gf_derSource));
      }
    }

    /*! \class AnisotropicLiquid
     *  \test \b AnisotropicLiquidUniformTest_gradient_AD tests the automatic
     * evaluation (full gradient)
     *  of the AnisotropicLiquid Green's function against analytical result for a
     * uniform dielectric
     */
    WHEN("the derivatives are evaluated via AD using the full gradient") {
      AnisotropicLiquid<AD_gradient> gf(epsilon, euler);
      THEN("the value of the Green's function is") {
        double value = result(0);
        double gf_value = gf.kernelS(source, probe);
        REQUIRE(value == Approx(gf_value));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "probe point is") {
        double derProbe = result(1);
        double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "source point is") {
        double derSource = result(2);
        double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
        REQUIRE(derSource == Approx(gf_derSource));
      }
    }

    /*! \class AnisotropicLiquid
     *  \test \b AnisotropicLiquidUniformTest_hessian_AD tests the automatic
     * evaluation (full hessian)
     *  of the AnisotropicLiquid Green's function against analytical result for a
     * uniform dielectric
     */
    WHEN("the derivatives are evaluated via AD using the full hessian") {
      AnisotropicLiquid<AD_hessian> gf(epsilon, euler);
      THEN("the value of the Green's function is") {
        double value = result(0);
        double gf_value = gf.kernelS(source, probe);
        REQUIRE(value == Approx(gf_value));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "probe point is") {
        double derProbe = result(1);
        double gf_derProbe = gf.derivativeProbe(probeNormal, source, probe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "source point is") {
        double derSource = result(2);
        double gf_derSource = gf.derivativeSource(sourceNormal, source, probe);
        REQUIRE(derSource == Approx(gf_derSource));
      }
      /*
         AND_THEN("the value of the Green's function hessian is")
         {
         double hessian = result(4);
         double gf_hessian = gf.hessian(sourceNormal, source, probeNormal, probe);
         REQUIRE(hessian == Approx(gf_hessian));
         }
         */
    }
  }
}
