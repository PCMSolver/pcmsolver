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

SCENARIO(
    "GePol cavity for the linear CO2 molecule in different Abelian point groups",
    "[gepol][gepol_CO2_symmetry]") {
  GIVEN("The linear CO2 molecule") {
    double area = 0.2 / bohr2ToAngstrom2();
    double probeRadius = 1.385 / bohrToAngstrom();
    double minRadius = 0.2 / bohrToAngstrom();
    WHEN("the point group is C1") {
      Molecule molec = CO2<0>();
      GePolCavity cavity(molec, area, probeRadius, minRadius, "co2.c1");

      /*! \class GePolCavity
       *  \test \b GePolCavityCO2C1Test_size tests GePol cavity size for CO2 in C1
       * symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 448;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2C1Test_irreducible_size tests GePol cavity
       * irreducible size for CO2 in C1 symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 448;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2C1Test_area tests GePol cavity surface area for CO2
       * in C1 symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 250.68176442433020;
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2C1Test_volume tests GePol cavity volume for CO2 in
       * C1 symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 352.55869984340751;
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

  GIVEN("The linear CO2 molecule") {
    double area = 0.2 / bohr2ToAngstrom2();
    double probeRadius = 1.385 / bohrToAngstrom();
    double minRadius = 0.2 / bohrToAngstrom();
    WHEN("the point group is C2") {
      Molecule molec = CO2<1>();
      GePolCavity cavity(molec, area, probeRadius, minRadius, "co2.c2");

      /*! \class GePolCavity
       *  \test \b GePolCavityCO2C1Test_size tests GePol cavity size for CO2 in C2
       * symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 448;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2C1Test_irreducible_size tests GePol cavity
       * irreducible size for CO2 in C2 symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 224;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2C1Test_area tests GePol cavity surface area for CO2
       * in C2 symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 250.68176442433020;
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2C1Test_volume tests GePol cavity volume for CO2 in
       * C2 symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 352.55869984340751;
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

  GIVEN("The linear CO2 molecule") {
    double area = 0.2 / bohr2ToAngstrom2();
    double probeRadius = 1.385 / bohrToAngstrom();
    double minRadius = 0.2 / bohrToAngstrom();
    WHEN("the point group is Cs") {
      Molecule molec = CO2<2>();
      GePolCavity cavity(molec, area, probeRadius, minRadius, "co2.cs");

      /*! \class GePolCavity
       *  \test \b GePolCavityCO2CsTest_size tests GePol cavity size for CO2 in Cs
       * symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 448;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2CsTest_irreducible_size tests GePol cavity
       * irreducible size for CO2 in Cs symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 224;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2CsTest_area tests GePol cavity surface area for CO2
       * in Cs symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 250.68176442433020;
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2CsTest_volume tests GePol cavity volume for CO2 in
       * Cs symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 352.55869984340751;
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

  GIVEN("The linear CO2 molecule") {
    double area = 0.2 / bohr2ToAngstrom2();
    double probeRadius = 1.385 / bohrToAngstrom();
    double minRadius = 0.2 / bohrToAngstrom();
    WHEN("the point group is Ci") {
      Molecule molec = CO2<3>();
      GePolCavity cavity(molec, area, probeRadius, minRadius, "co2.ci");

      /*! \class GePolCavity
       *  \test \b GePolCavityCO2CiTest_size tests GePol cavity size for CO2 in Ci
       * symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 448;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2CiTest_irreducible_size tests GePol cavity
       * irreducible size for CO2 in Ci symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 224;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2CiTest_area tests GePol cavity surface area for CO2
       * in Ci symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 250.68176442433020;
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2CiTest_volume tests GePol cavity volume for CO2 in
       * Ci symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 352.55869984340751;
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

  GIVEN("The linear CO2 molecule") {
    double area = 0.2 / bohr2ToAngstrom2();
    double probeRadius = 1.385 / bohrToAngstrom();
    double minRadius = 0.2 / bohrToAngstrom();
    WHEN("the point group is D2") {
      Molecule molec = CO2<4>();
      GePolCavity cavity(molec, area, probeRadius, minRadius, "co2.d2");

      /*! \class GePolCavity
       *  \test \b GePolCavityCO2D2Test_size tests GePol cavity size for CO2 in D2
       * symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 448;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2D2Test_irreducible_size tests GePol cavity
       * irreducible size for CO2 in D2 symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 112;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2D2Test_area tests GePol cavity surface area for CO2
       * in D2 symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 250.68176442433020;
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2D2Test_volume tests GePol cavity volume for CO2 in
       * D2 symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 352.55869984340751;
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

  GIVEN("The linear CO2 molecule") {
    double area = 0.2 / bohr2ToAngstrom2();
    double probeRadius = 1.385 / bohrToAngstrom();
    double minRadius = 0.2 / bohrToAngstrom();
    WHEN("the point group is C2v") {
      Molecule molec = CO2<5>();
      GePolCavity cavity(molec, area, probeRadius, minRadius, "co2.c2v");

      /*! \class GePolCavity
       *  \test \b GePolCavityCO2C2vTest_size tests GePol cavity size for CO2 in C2v
       * symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 448;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2C2vTest_irreducible_size tests GePol cavity
       * irreducible size for CO2 in C2v symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 112;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2C2vTest_area tests GePol cavity surface area for CO2
       * in C2v symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 250.68176442433020;
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2C2vTest_volume tests GePol cavity volume for CO2 in
       * C2v symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 352.55869984340751;
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

  GIVEN("The linear CO2 molecule") {
    double area = 0.2 / bohr2ToAngstrom2();
    double probeRadius = 1.385 / bohrToAngstrom();
    double minRadius = 0.2 / bohrToAngstrom();
    WHEN("the point group is C2h") {
      Molecule molec = CO2<6>();
      GePolCavity cavity(molec, area, probeRadius, minRadius, "co2.c2h");

      /*! \class GePolCavity
       *  \test \b GePolCavityCO2C2hTest_size tests GePol cavity size for CO2 in C2h
       * symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 448;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2C2hTest_irreducible_size tests GePol cavity
       * irreducible size for CO2 in C2h symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 112;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2C2hTest_area tests GePol cavity surface area for CO2
       * in C2h symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 250.68176442433020;
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2C2hTest_volume tests GePol cavity volume for CO2 in
       * C2h symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 352.55869984340751;
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

  GIVEN("The linear CO2 molecule") {
    double area = 0.2 / bohr2ToAngstrom2();
    double probeRadius = 1.385 / bohrToAngstrom();
    double minRadius = 0.2 / bohrToAngstrom();
    WHEN("the point group is D2h") {
      Molecule molec = CO2<7>();
      GePolCavity cavity(molec, area, probeRadius, minRadius, "co2.d2h");

      /*! \class GePolCavity
       *  \test \b GePolCavityCO2D2hTest_size tests GePol cavity size for CO2 in D2h
       * symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 448;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2D2hTest_irreducible_size tests GePol cavity
       * irreducible size for CO2 in D2h symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 56;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2D2hTest_area tests GePol cavity surface area for CO2
       * in D2h symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 250.68176442433020;
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCO2D2hTest_volume tests GePol cavity volume for CO2 in
       * D2h symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 352.55869984340751;
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
