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

#include <Eigen/Core>

#include "TestingMolecules.hpp"
#include "cavity/GePolCavity.hpp"
#include "utils/Molecule.hpp"
#include "utils/Symmetry.hpp"

using namespace pcm;
using cavity::GePolCavity;

SCENARIO("GePol cavity for the H3+ molecule in C2v symmetry",
         "[gepol][gepol_H3+_C2v]") {
  GIVEN("The H3+ molecule in C2v symmetry") {
    Molecule molec = H3<5>();

    WHEN("the addition of spheres is enabled") {
      double area = 0.2 / bohr2ToAngstrom2();
      double probeRadius = 1.385 / bohrToAngstrom();
      double minRadius = 0.2 / bohrToAngstrom();
      GePolCavity cavity(molec, area, probeRadius, minRadius, "c2v");
      cavity.saveCavity("h3+_c2v.npz");

      /*! \class GePolCavity
       *  \test \b GePolCavityC2vAddTest_size tests GePol cavity size for H3+ in C2v
       * symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 312;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC2vAddTest_irreducible_size tests GePol cavity
       * irreducible size for H3+ in C2v symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 78;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC2vAddTest_area tests GePol cavity surface area for H3+
       * in C2v symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 178.74700256128352;
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC2vAddTest_volume tests GePol cavity volume for H3+ in
       * C2v symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 196.47360294559090;
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
      GePolCavity cavity(molec, area, probeRadius, minRadius, "c2v_noadd");
      cavity.saveCavity("h3+_c2v_noadd.npz");

      /*! \class GePolCavity
       *  \test \b GePolCavityC2vTest_size tests GePol cavity size for H3+ in C2v
       * symmetry without added spheres
       */
      THEN("the size of the cavity is") {
        int size = 288;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC2vTest_irreducible_size tests GePol cavity irreducible
       * size for H3+ in C2v symmetry without added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 72;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC2vTest_area tests GePol cavity surface area for H3+ in
       * C2v symmetry without added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 181.87043332808548;
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC2vTest_volume tests GePol cavity volume for H3+ in C2v
       * symmetry without added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 192.48281460140359;
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
