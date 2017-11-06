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
#include <vector>

#include <Eigen/Core>

#include "TestingMolecules.hpp"
#include "cavity/GePolCavity.hpp"
#include "green/IonicLiquid.hpp"
#include "green/UniformDielectric.hpp"
#include "green/Vacuum.hpp"
#include "utils/ChargeDistribution.hpp"

using namespace pcm;
using cavity::GePolCavity;
using green::Vacuum;
using green::UniformDielectric;
using green::IonicLiquid;
using utils::ChargeDistribution;

SCENARIO("Calculation of the dipolar potential",
         "[utils][dipolar_potential][utils_dipolar-potential]") {
  GIVEN("A classical charge distribution of dipoles") {
    double radius = 1.181 * 1.10 / bohrToAngstrom();
    Sphere point(Eigen::Vector3d::Zero(), radius);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);
    ChargeDistribution dist;
    dist.dipoles = Eigen::Vector3d::UnitZ();
    dist.dipolesSites = Eigen::Vector3d::Zero();

    Eigen::VectorXd dipolar = Eigen::VectorXd::Zero(cavity.size());
    for (size_t i = 0; i < cavity.size(); ++i) {
      Eigen::Vector3d center = cavity.elementCenter(i);
      dipolar(i) = center(2) / std::pow(center.norm(), 3);
    }

    WHEN("the dipolar potential is calculated in vacuum") {
      Vacuum<> gf;
      Eigen::VectorXd pot = computeDipolarPotential(
          gf.exportDerivativeProbe(), cavity.elementCenter(), dist);
      Eigen::VectorXd alt = computeDipolarPotential(cavity.elementCenter(), dist);
      THEN("comparison with the analytical results in") {
        REQUIRE(pot.size() == dipolar.size());
        for (int i = 0; i < dipolar.size(); ++i) {
          REQUIRE(pot(i) == Approx(dipolar(i)));
          REQUIRE(alt(i) == Approx(dipolar(i)));
        }
      }
    }

    AND_WHEN("the dipolar potential is calculated in a uniform dielectric") {
      double permittivity = 78.39;
      UniformDielectric<> gf(permittivity);
      Eigen::VectorXd pot = computeDipolarPotential(
          gf.exportDerivativeProbe(), cavity.elementCenter(), dist);
      THEN("comparison with the analytical results in") {
        dipolar /= permittivity;
        REQUIRE(pot.size() == dipolar.size());
        for (int i = 0; i < dipolar.size(); ++i) {
          REQUIRE(pot(i) == Approx(dipolar(i)));
        }
      }
    }

    AND_WHEN(
        "the dipolar potential is calculated in an ionic liquid with kappa = 0.0") {
      double permittivity = 78.39;
      double kappa = 0.0;
      IonicLiquid<> gf(permittivity, kappa);
      Eigen::VectorXd pot = computeDipolarPotential(
          gf.exportDerivativeProbe(), cavity.elementCenter(), dist);
      THEN("comparison with the analytical results in") {
        dipolar /= permittivity;
        REQUIRE(pot.size() == dipolar.size());
        for (int i = 0; i < dipolar.size(); ++i) {
          REQUIRE(pot(i) == Approx(dipolar(i)));
        }
      }
    }
  }
}
