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

#include <iostream>

#include <Eigen/Core>

#include "TestingMolecules.hpp"
#include "bi_operators/Collocation.hpp"
#include "cavity/GePolCavity.hpp"
#include "green/DerivativeTypes.hpp"
#include "green/SphericalSharp.hpp"
#include "green/Vacuum.hpp"
#include "solver/IEFSolver.hpp"
#include "utils/Molecule.hpp"

using namespace pcm;
using bi_operators::Collocation;
using cavity::GePolCavity;
using green::SphericalSharp;
using green::Vacuum;
using solver::IEFSolver;

SCENARIO("Test solver for the IEFPCM for a point charge in a spherical sharp "
         "environment and a GePol cavity",
         "[solver][iefpcm][iefpcm_sharp-gepol-point][anisotropic]") {
  GIVEN("A spherical nanoparticle modelled as a spherical sharp permittivity") {
    double epsNP = 114.0;
    double epsSolv = 35.7;
    double sphereRadius = 100.0;
    int maxL = 200;
    Eigen::Vector3d offset;
    offset << 105.0, 106.0, 107.0;

    Vacuum<> gf_i;
    bool symm = true;
    Collocation op;

    double charge = 8.0;
    double totalASC = -7.7727501347;
    /*! \class IEFSolver
     *  \test \b pointChargeSharpGePol tests IEFSolver using a point charge with a
     * GePol cavity and a spherical sharp interface
     *  The spherical sharp interface is centered at the origin, while the point
     * charge is away from the origin.
     */
    WHEN("the spherical sharp layer is centered at the origin and the charge is "
         "away from the origin") {
      Molecule point = dummy<0>(2.0, offset);
      double area = 1.0;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius);

      SphericalSharp<> gf_o(
          epsNP, epsSolv, sphereRadius, Eigen::Vector3d::Zero(), maxL);
      IEFSolver solver(symm);
      solver.buildSystemMatrix(cavity, gf_i, gf_o, op);
      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge, offset);
      for (int i = 0; i < size; ++i) {
        INFO("fake_mep(" << i << ") = " << fake_mep(i));
      }
      THEN("the apparent surface charge is") {
        Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
        fake_asc = solver.computeCharge(fake_mep);
        double totalFakeASC = fake_asc.sum();
        for (int i = 0; i < size; ++i) {
          INFO("fake_asc(" << i << ") = " << fake_asc(i));
        }

        CAPTURE(totalASC);
        CAPTURE(totalFakeASC);
        CAPTURE(totalASC - totalFakeASC);
        REQUIRE(totalASC == Approx(totalFakeASC));
      }
    }

    AND_WHEN("the spherical sharp layer is centered away from the origin and the "
             "charge is at the origin") {
      Molecule point = dummy<0>(2.0);
      double area = 1.0;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius);

      SphericalSharp<> gf_o(epsNP, epsSolv, sphereRadius, offset, maxL);

      IEFSolver solver(symm);
      solver.buildSystemMatrix(cavity, gf_i, gf_o, op);

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);
      for (int i = 0; i < size; ++i) {
        INFO("fake_mep(" << i << ") = " << fake_mep(i));
      }
      THEN("the apparent surface charge is") {
        Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
        fake_asc = solver.computeCharge(fake_mep);
        double totalFakeASC = fake_asc.sum();
        for (int i = 0; i < size; ++i) {
          INFO("fake_asc(" << i << ") = " << fake_asc(i));
        }

        CAPTURE(totalASC);
        CAPTURE(totalFakeASC);
        CAPTURE(totalASC - totalFakeASC);
        REQUIRE(totalASC == Approx(totalFakeASC));
      }
    }
  }
}
