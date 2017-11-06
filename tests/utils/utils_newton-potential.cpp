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
using utils::nuclearChargeDistribution;

SCENARIO("Calculation of the Newton potential",
         "[utils][newton_potential][utils_newton-potential]") {
  GIVEN("A classical charge distribution of monopoles") {
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    Molecule molec = C2H4();
    GePolCavity cavity(molec, area, probeRadius, minRadius, "newton");
    ChargeDistribution dist = nuclearChargeDistribution(molec);

    WHEN("the Newton potential is calculated in vacuum") {
      Vacuum<> gf;
      Eigen::VectorXd newton =
          computeNewtonPotential(gf.exportKernelS(), cavity.elementCenter(), dist);
      THEN("comparison with the computeMEP method results in") {
        Eigen::VectorXd mep = computeMEP(molec, cavity.elements());
        REQUIRE(newton.size() == mep.size());
        for (int i = 0; i < mep.size(); ++i) {
          REQUIRE(newton(i) == Approx(mep(i)));
        }
      }
    }

    AND_WHEN("the Newton potential is calculated in a uniform dielectric") {
      double permittivity = 78.39;
      UniformDielectric<> gf(permittivity);
      Eigen::VectorXd newton =
          computeNewtonPotential(gf.exportKernelS(), cavity.elementCenter(), dist);
      THEN("comparison with the computeMEP method results in") {
        Eigen::VectorXd mep = computeMEP(molec, cavity.elements());
        mep /= permittivity;
        REQUIRE(newton.size() == mep.size());
        for (int i = 0; i < mep.size(); ++i) {
          REQUIRE(newton(i) == Approx(mep(i)));
        }
      }
    }

    AND_WHEN(
        "the Newton potential is calculated in an ionic liquid with kappa = 0.0") {
      double permittivity = 78.39;
      double kappa = 0.0;
      IonicLiquid<> gf(permittivity, kappa);
      Eigen::VectorXd newton =
          computeNewtonPotential(gf.exportKernelS(), cavity.elementCenter(), dist);
      THEN("comparison with the computeMEP method results in") {
        Eigen::VectorXd mep = computeMEP(molec, cavity.elements());
        mep /= permittivity;
        REQUIRE(newton.size() == mep.size());
        for (int i = 0; i < mep.size(); ++i) {
          REQUIRE(newton(i) == Approx(mep(i)));
        }
      }
    }
  }
}
