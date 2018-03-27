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

#include "catch.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#include <Eigen/Core>

#include "green/DerivativeTypes.hpp"
#include "green/SphericalDiffuse.hpp"
#include "green/dielectric_profile/OneLayerErf.hpp"
#include "green/dielectric_profile/OneLayerTanh.hpp"

using namespace pcm;
using dielectric_profile::OneLayerErf;
using dielectric_profile::OneLayerTanh;
using green::SphericalDiffuse;

SCENARIO("Evaluation of the spherical diffuse Green's function and its derivatives",
         "[green][green_spherical_diffuse]") {
  int maxL = 10;
  // High dielectric constant inside
  double eps1 = 80.0;
  // Low dielectric constant outside
  double eps2 = 2.0;
  double sphereRadius = 100.0;
  double width = 5.0;
  // Evaluation inside the sphere
  Eigen::Vector3d source1 = (Eigen::Vector3d() << 1.0, 0.0, 0.0).finished();
  Eigen::Vector3d sourceNormal1 = source1;
  sourceNormal1.normalize();
  Eigen::Vector3d probe1 = (Eigen::Vector3d() << 2.0, 0.0, 0.0).finished();
  Eigen::Vector3d probeNormal1 = probe1;
  probeNormal1.normalize();
  // Evaluation outside the sphere
  Eigen::Vector3d source2 = (Eigen::Vector3d() << 150.0, 150.0, 150.0).finished();
  Eigen::Vector3d sourceNormal2 = source2;
  sourceNormal2.normalize();
  Eigen::Vector3d probe2 = (Eigen::Vector3d() << 151.0, 150.0, 150.0).finished();
  Eigen::Vector3d probeNormal2 = probe2;
  probeNormal2.normalize();
  GIVEN("A permittivity profile modelled by the logarithmic profile") {
    WHEN("the spherical droplet is centered at the origin") {
      Eigen::Vector3d sphereCenter = Eigen::Vector3d::Zero();
      SphericalDiffuse<> gf(eps1, eps2, width, sphereRadius, sphereCenter, maxL);
      THEN("the value of the Green's function inside the droplet is") {
        double value = 0.0204749331147992819;
        double gf_value = gf.kernelS(source1, probe1);
        INFO("ref_value = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << value);
        INFO("gf_value  = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << gf_value);
        REQUIRE(value == Approx(gf_value));
      }
      AND_THEN("the value of the Green's function outside the droplet is") {
        double value = 0.487473605103865837;
        double gf_value = gf.kernelS(source2, probe2);
        INFO("ref_value = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << value);
        INFO("gf_value  = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << gf_value);
        REQUIRE(value == Approx(gf_value));
      }
      THEN("the value of the Green's function directional derivative wrt the probe "
           "point inside the droplet is") {
        double derProbe = -0.0124999765637075211;
        double gf_derProbe = gf.derivativeProbe(probeNormal1, source1, probe1);
        INFO("ref_derProbe = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derProbe);
        INFO("gf_derProbe  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derProbe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "probe point outside the droplet is") {
        double derProbe = -0.290662955233778053;
        double gf_derProbe = gf.derivativeProbe(probeNormal2, source2, probe2);
        INFO("ref_derProbe = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derProbe);
        INFO("gf_derProbe  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derProbe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      THEN("the value of the Green's function directional derivative wrt the source "
           "point inside the droplet is") {
        double derSource = 0.00937504723101320603;
        double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
        INFO("ref_derSource = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derSource);
        INFO("gf_derSource  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derSource);
        REQUIRE(derSource == Approx(gf_derSource));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "source point outside the droplet is") {
        double derSource = -0.901954472843635724;
        double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
        INFO("ref_derSource = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derSource);
        INFO("gf_derSource  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derSource);
        REQUIRE(derSource == Approx(gf_derSource));
      }
    }

    AND_WHEN("the spherical droplet is centered away from the origin") {
      Eigen::Vector3d sphereCenter =
          (Eigen::Vector3d() << 25.0, 0.0, 0.0).finished();
      SphericalDiffuse<> gf(eps1, eps2, width, sphereRadius, sphereCenter, maxL);
      THEN("the value of the Green's function inside the droplet is") {
        double value = 0.0173845118865097424;
        double gf_value = gf.kernelS(source1, probe1);
        INFO("ref_value = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << value);
        INFO("gf_value  = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << gf_value);
        REQUIRE(value == Approx(gf_value));
      }
      AND_THEN("the value of the Green's function outside the droplet is") {
        double value = 0.501385603574595162;
        double gf_value = gf.kernelS(source2, probe2);
        INFO("ref_value = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << value);
        INFO("gf_value  = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << gf_value);
        REQUIRE(value == Approx(gf_value));
      }
      THEN("the value of the Green's function directional derivative wrt the probe "
           "point inside the droplet is") {
        double derProbe = -0.0124834385275933435;
        double gf_derProbe = gf.derivativeProbe(probeNormal1, source1, probe1);
        INFO("ref_derProbe = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derProbe);
        INFO("gf_derProbe  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derProbe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "probe point outside the droplet is") {
        double derProbe = -0.291262383059764929;
        double gf_derProbe = gf.derivativeProbe(probeNormal2, source2, probe2);
        INFO("ref_derProbe = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derProbe);
        INFO("gf_derProbe  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derProbe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      THEN("the value of the Green's function directional derivative wrt the source "
           "point inside the droplet is") {
        double derSource = 0.012483564647454426;
        double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
        INFO("ref_derSource = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derSource);
        INFO("gf_derSource  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derSource);
        REQUIRE(derSource == Approx(gf_derSource));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "source point outside the droplet is") {
        double derSource = 2.720604738008503;
        double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
        INFO("ref_derSource = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derSource);
        INFO("gf_derSource  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derSource);
        REQUIRE(derSource == Approx(gf_derSource));
      }
    }
  }
  GIVEN("A permittivity profile modelled by the hyperbolic tangent function") {
    WHEN("the spherical droplet is centered at the origin") {
      Eigen::Vector3d sphereCenter = Eigen::Vector3d::Zero();
      SphericalDiffuse<OneLayerTanh> gf(
          eps1, eps2, width, sphereRadius, sphereCenter, maxL);
      THEN("the value of the Green's function inside the droplet is") {
        double value = 0.0204265162808963327;
        double gf_value = gf.kernelS(source1, probe1);
        INFO("ref_value = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << value);
        INFO("gf_value  = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << gf_value);
        REQUIRE(value == Approx(gf_value));
      }
      AND_THEN("the value of the Green's function outside the droplet is") {
        double value = 0.42831756921770181;
        double gf_value = gf.kernelS(source2, probe2);
        INFO("ref_value = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << value);
        INFO("gf_value  = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << gf_value);
        REQUIRE(value == Approx(gf_value));
      }
      THEN("the value of the Green's function directional derivative wrt the probe "
           "point inside the droplet is") {
        double derProbe = -0.0124999769448609632;
        double gf_derProbe = gf.derivativeProbe(probeNormal1, source1, probe1);
        INFO("ref_derProbe = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derProbe);
        INFO("gf_derProbe  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derProbe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "probe point outside the droplet is") {
        double derProbe = -0.288392889131738883;
        double gf_derProbe = gf.derivativeProbe(probeNormal2, source2, probe2);
        INFO("ref_derProbe = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derProbe);
        INFO("gf_derProbe  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derProbe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      THEN("the value of the Green's function directional derivative wrt the source "
           "point inside the droplet is") {
        double derSource = 0.00937504646828998811;
        double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
        INFO("ref_derSource = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derSource);
        INFO("gf_derSource  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derSource);
        REQUIRE(derSource == Approx(gf_derSource));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "source point outside the droplet is") {
        double derSource = -1.51576208806569745;
        double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
        INFO("ref_derSource = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derSource);
        INFO("gf_derSource  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derSource);
        REQUIRE(derSource == Approx(gf_derSource));
      }
    }

    AND_WHEN("the spherical droplet is centered away from the origin") {
      Eigen::Vector3d sphereCenter =
          (Eigen::Vector3d() << 25.0, 0.0, 0.0).finished();
      SphericalDiffuse<OneLayerTanh> gf(
          eps1, eps2, width, sphereRadius, sphereCenter, maxL);
      THEN("the value of the Green's function inside the droplet is") {
        double value = 0.017335802360862948;
        double gf_value = gf.kernelS(source1, probe1);
        INFO("ref_value = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << value);
        INFO("gf_value  = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << gf_value);
        REQUIRE(value == Approx(gf_value));
      }
      AND_THEN("the value of the Green's function outside the droplet is") {
        double value = 0.470725113611886292;
        double gf_value = gf.kernelS(source2, probe2);
        INFO("ref_value = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << value);
        INFO("gf_value  = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << gf_value);
        REQUIRE(value == Approx(gf_value));
      }
      THEN("the value of the Green's function directional derivative wrt the probe "
           "point inside the droplet is") {
        double derProbe = -0.0124834504460630269;
        double gf_derProbe = gf.derivativeProbe(probeNormal1, source1, probe1);
        INFO("ref_derProbe = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derProbe);
        INFO("gf_derProbe  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derProbe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "probe point outside the droplet is") {
        double derProbe = -0.290870450536662162;
        double gf_derProbe = gf.derivativeProbe(probeNormal2, source2, probe2);
        INFO("ref_derProbe = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derProbe);
        INFO("gf_derProbe  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derProbe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      THEN("the value of the Green's function directional derivative wrt the source "
           "point inside the droplet is") {
        double derSource = 0.012483552270584064;
        double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
        INFO("ref_derSource = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derSource);
        INFO("gf_derSource  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derSource);
        REQUIRE(derSource == Approx(gf_derSource));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "source point outside the droplet is") {
        double derSource = 3.75908425086779463;
        double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
        INFO("ref_derSource = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derSource);
        INFO("gf_derSource  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derSource);
        REQUIRE(derSource == Approx(gf_derSource));
      }
    }
  }

  GIVEN("A permittivity profile modelled by the error function") {
    WHEN("the spherical droplet is centered at the origin") {
      Eigen::Vector3d sphereCenter = Eigen::Vector3d::Zero();
      SphericalDiffuse<OneLayerErf> gf(
          eps1, eps2, width, sphereRadius, sphereCenter, maxL);
      THEN("the value of the Green's function inside the droplet is") {
        double value = 0.0204465601345278428;
        double gf_value = gf.kernelS(source1, probe1);
        INFO("ref_value = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << value);
        INFO("gf_value  = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << gf_value);
        REQUIRE(value == Approx(gf_value));
      }
      AND_THEN("the value of the Green's function outside the droplet is") {
        double value = 0.546440841068159378;
        double gf_value = gf.kernelS(source2, probe2);
        INFO("ref_value = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << value);
        INFO("gf_value  = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << gf_value);
        REQUIRE(value == Approx(gf_value));
      }
      THEN("the value of the Green's function directional derivative wrt the probe "
           "point inside the droplet is") {
        double derProbe = -0.0124999769314168563;
        double gf_derProbe = gf.derivativeProbe(probeNormal1, source1, probe1);
        INFO("ref_derProbe = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derProbe);
        INFO("gf_derProbe  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derProbe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "probe point outside the droplet is") {
        double derProbe = -0.292045737975143993;
        double gf_derProbe = gf.derivativeProbe(probeNormal2, source2, probe2);
        INFO("ref_derProbe = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derProbe);
        INFO("gf_derProbe  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derProbe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      THEN("the value of the Green's function directional derivative wrt the source "
           "point inside the droplet is") {
        double derSource = 0.00937504649523024369;
        double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
        INFO("ref_derSource = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derSource);
        INFO("gf_derSource  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derSource);
        REQUIRE(derSource == Approx(gf_derSource));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "source point outside the droplet is") {
        double derSource = -1.60173403890240262;
        double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
        INFO("ref_derSource = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derSource);
        INFO("gf_derSource  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derSource);
        REQUIRE(derSource == Approx(gf_derSource));
      }
    }

    AND_WHEN("the spherical droplet is centered away from the origin") {
      Eigen::Vector3d sphereCenter =
          (Eigen::Vector3d() << 25.0, 0.0, 0.0).finished();
      SphericalDiffuse<OneLayerErf> gf(
          eps1, eps2, width, sphereRadius, sphereCenter, maxL);
      THEN("the value of the Green's function inside the droplet is") {
        double value = 0.0173558563009676217;
        double gf_value = gf.kernelS(source1, probe1);
        INFO("ref_value = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << value);
        INFO("gf_value  = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << gf_value);
        REQUIRE(value == Approx(gf_value));
      }
      AND_THEN("the value of the Green's function outside the droplet is") {
        double value = 0.37238664692877399;
        double gf_value = gf.kernelS(source2, probe2);
        INFO("ref_value = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << value);
        INFO("gf_value  = " << std::setprecision(
                                   std::numeric_limits<long double>::digits10)
                            << gf_value);
        REQUIRE(value == Approx(gf_value));
      }
      THEN("the value of the Green's function directional derivative wrt the probe "
           "point inside the droplet is") {
        double derProbe = -0.0124834500456541542;
        double gf_derProbe = gf.derivativeProbe(probeNormal1, source1, probe1);
        INFO("ref_derProbe = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derProbe);
        INFO("gf_derProbe  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derProbe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "probe point outside the droplet is") {
        double derProbe = -0.279220523679291066;
        double gf_derProbe = gf.derivativeProbe(probeNormal2, source2, probe2);
        INFO("ref_derProbe = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derProbe);
        INFO("gf_derProbe  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derProbe);
        REQUIRE(derProbe == Approx(gf_derProbe));
      }
      THEN("the value of the Green's function directional derivative wrt the source "
           "point inside the droplet is") {
        double derSource = 0.0124835526869870872;
        double gf_derSource = gf.derivativeSource(sourceNormal1, source1, probe1);
        INFO("ref_derSource = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derSource);
        INFO("gf_derSource  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derSource);
        REQUIRE(derSource == Approx(gf_derSource));
      }
      AND_THEN("the value of the Green's function directional derivative wrt the "
               "source point outside the droplet is") {
        double derSource = 4.79631817304426722;
        double gf_derSource = gf.derivativeSource(sourceNormal2, source2, probe2);
        INFO("ref_derSource = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << derSource);
        INFO("gf_derSource  = "
             << std::setprecision(std::numeric_limits<long double>::digits10)
             << gf_derSource);
        REQUIRE(derSource == Approx(gf_derSource));
      }
    }
  }
}
