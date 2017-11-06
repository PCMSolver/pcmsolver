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

#include <catch.hpp>

#include <cmath>
#include <vector>

#include <Eigen/Core>

#include "TestingMolecules.hpp"
#include "cavity/GePolCavity.hpp"
#include "utils/Molecule.hpp"
#include "utils/Symmetry.hpp"

using namespace pcm;
using cavity::GePolCavity;

SCENARIO("GePol cavity for the C2H4 molecule in D2h symmetry",
         "[gepol][gepol_C2H4_D2h]") {
  GIVEN("The C2H4 molecule in D2h symmetry") {
    Molecule molec = C2H4();

    WHEN("the addition of spheres is enabled") {
      double area = 0.2 / bohr2ToAngstrom2();
      double probeRadius = 1.385 / bohrToAngstrom();
      double minRadius = 0.2 / bohrToAngstrom();
      GePolCavity cavity(molec, area, probeRadius, minRadius, "d2h");
      cavity.saveCavity("c2h4_d2h.npz");

      /*! \class GePolCavity
       *  \test \b GePolCavityD2hAddTest_size tests GePol cavity size for C2H4 in D2h
       * symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 576;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityD2hAddTest_irreducible_size tests GePol cavity
       * irreducible size for C2H4 in D2h symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 72;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityD2hAddTest_area tests GePol cavity surface area for
       * C2H4 in D2h symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 281.81993683500656;
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityD2hAddTest_volume tests GePol cavity volume for C2H4 in
       * D2h symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 406.54737252764619;
        Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
        Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
        double actualVolume = 0;
        for (int i = 0; i < cavity.size(); ++i) {
          actualVolume +=
              cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
        }
        actualVolume /= 3;
        REQUIRE(volume == Approx(actualVolume));
      }
    }

    WHEN("the addition of spheres is disabled") {
      double area = 0.2 / bohr2ToAngstrom2();
      double probeRadius = 1.385 / bohrToAngstrom();
      double minRadius = 100.0 / bohrToAngstrom();
      GePolCavity cavity(molec, area, probeRadius, minRadius, "d2h_noadd");
      cavity.saveCavity("c2h4_d2h_noadd.npz");

      /*! \class GePolCavity
       *  \test \b GePolCavityD2hTest_size tests GePol cavity size for C2H4 in D2h
       * symmetry without added spheres
       */
      THEN("the size of the cavity is") {
        int size = 576;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityD2hTest_irreducible_size tests GePol cavity irreducible
       * size for C2H4 in D2h symmetry without added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 72;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityD2hTest_area tests GePol cavity surface area for C2H4
       * in D2h symmetry without added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 281.81993683500656;
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityD2hTest_volume tests GePol cavity volume for C2H4 in
       * D2h symmetry without added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 406.54737252764619;
        Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
        Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
        double actualVolume = 0;
        for (int i = 0; i < cavity.size(); ++i) {
          actualVolume +=
              cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
        }
        actualVolume /= 3;
        REQUIRE(volume == Approx(actualVolume));
      }
    }
  }
}
