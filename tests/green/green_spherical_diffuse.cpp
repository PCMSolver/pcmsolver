/**
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
#include "green/dielectric_profile/MembraneTanh.hpp"
#include "green/dielectric_profile/OneLayerErf.hpp"
#include "green/dielectric_profile/OneLayerTanh.hpp"

using namespace pcm;
using green::SphericalDiffuse;
using dielectric_profile::OneLayerErf;
using dielectric_profile::OneLayerTanh;
using dielectric_profile::MembraneTanh;

SCENARIO("Evaluation of the spherical diffuse Green's function and its derivatives",
         "[green][green_spherical_diffuse]") {
  GIVEN("A permittivity profile modelled by the hyperbolic tangent function") {
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
    WHEN("the spherical droplet is centered at the origin") {
      Eigen::Vector3d sphereCenter = Eigen::Vector3d::Zero();
      SphericalDiffuse<> gf(eps1, eps2, width, sphereRadius, sphereCenter, maxL);
      THEN("the value of the Green's function inside the droplet is") {
        double value = 0.0125073117820161759;
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
        double value = 0.49987476890456245;
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
        double derProbe = -0.0125062279977827817;
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
        double derProbe = -0.289952813654903441;
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
        double derSource = 0.0124986222807023373;
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
        double derSource = 0.288676184370117994;
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
        double value = 0.0125233347952455679;
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
        double value = 0.499902162204405809;
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
        double derProbe = -0.0125024378794594271;
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
        double derProbe = -0.289953098676354326;
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
        double derSource = 0.0124899030161038094;
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
        double derSource = 0.288675899050294671;
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
    int maxL = 3;
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
    WHEN("the spherical droplet is centered at the origin") {
      Eigen::Vector3d sphereCenter = Eigen::Vector3d::Zero();
      SphericalDiffuse<OneLayerErf> gf(
          eps1, eps2, width, sphereRadius, sphereCenter, maxL);
      THEN("the value of the Green's function inside the droplet is") {
        double value = 0.0125073115838534744;
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
        double value = 0.499874677668854961;
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
        double derProbe = -0.0125062686709680876;
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
        double derProbe = -0.289951551978862021;
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
        double derSource = 0.0124986220027042266;
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
        double derSource = 0.288674990983894819;
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
        double value = 0.0125233447844991293;
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
        double value = 0.499897402784606149;
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
        double derProbe = -0.0125024373412441214;
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
        double derProbe = -0.289953249665575452;
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
        double derSource = 0.0124899033273999371;
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
        double derSource = 0.288675321960529807;
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
