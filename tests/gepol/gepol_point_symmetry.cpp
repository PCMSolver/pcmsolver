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

using namespace pcm;
using cavity::GePolCavity;

SCENARIO("GePol cavity for a single sphere in different Abelian point groups",
         "[gepol][gepol_point_symmetry]") {
  GIVEN("A single sphere") {
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    WHEN("the point group is C1") {
      Molecule point = dummy<0>();
      GePolCavity cavity(point, area, probeRadius, minRadius, "c1");

      /*! \class GePolCavity
       *  \test \b GePolCavityC1Test_size tests GePol cavity size for a point charge
       * in C1 symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 32;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC1Test_irreducible_size tests GePol cavity irreducible
       * size for a point charge in C1 symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 32;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC1Test_area tests GePol cavity surface area for a point
       * charge in C1 symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 4.0 * M_PI * pow(1.0, 2);
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC1Test_volume tests GePol cavity volume for a point
       * charge in C1 symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
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

  GIVEN("A single sphere") {
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    WHEN("the point group is C2") {
      Molecule point = dummy<1>();
      GePolCavity cavity(point, area, probeRadius, minRadius, "c2");

      /*! \class GePolCavity
       *  \test \b GePolCavityC1Test_size tests GePol cavity size for a point charge
       * in C2 symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 32;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC1Test_irreducible_size tests GePol cavity irreducible
       * size for a point charge in C2 symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 16;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC1Test_area tests GePol cavity surface area for a point
       * charge in C2 symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 4.0 * M_PI * pow(1.0, 2);
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC1Test_volume tests GePol cavity volume for a point
       * charge in C2 symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
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

  GIVEN("A single sphere") {
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    WHEN("the point group is Cs") {
      Molecule point = dummy<2>();
      GePolCavity cavity(point, area, probeRadius, minRadius, "cs");

      /*! \class GePolCavity
       *  \test \b GePolCavityCsTest_size tests GePol cavity size for a point charge
       * in Cs symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 32;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCsTest_irreducible_size tests GePol cavity irreducible
       * size for a point charge in Cs symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 16;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCsTest_area tests GePol cavity surface area for a point
       * charge in Cs symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 4.0 * M_PI * pow(1.0, 2);
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCsTest_volume tests GePol cavity volume for a point
       * charge in Cs symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
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

  GIVEN("A single sphere") {
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    WHEN("the point group is Ci") {
      Molecule point = dummy<3>();
      GePolCavity cavity(point, area, probeRadius, minRadius, "ci");

      /*! \class GePolCavity
       *  \test \b GePolCavityCiTest_size tests GePol cavity size for a point charge
       * in Ci symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 32;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCiTest_irreducible_size tests GePol cavity irreducible
       * size for a point charge in Ci symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 16;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCiTest_area tests GePol cavity surface area for a point
       * charge in Ci symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 4.0 * M_PI * pow(1.0, 2);
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityCiTest_volume tests GePol cavity volume for a point
       * charge in Ci symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
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

  GIVEN("A single sphere") {
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    WHEN("the point group is D2") {
      Molecule point = dummy<4>();
      GePolCavity cavity(point, area, probeRadius, minRadius, "d2");

      /*! \class GePolCavity
       *  \test \b GePolCavityD2Test_size tests GePol cavity size for a point charge
       * in D2 symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 32;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityD2Test_irreducible_size tests GePol cavity irreducible
       * size for a point charge in D2 symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 8;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityD2Test_area tests GePol cavity surface area for a point
       * charge in D2 symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 4.0 * M_PI * pow(1.0, 2);
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityD2Test_volume tests GePol cavity volume for a point
       * charge in D2 symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
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

  GIVEN("A single sphere") {
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    WHEN("the point group is C2v") {
      Molecule point = dummy<5>();
      GePolCavity cavity(point, area, probeRadius, minRadius, "c2v");

      /*! \class GePolCavity
       *  \test \b GePolCavityC2vTest_size tests GePol cavity size for a point charge
       * in C2v symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 32;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC2vTest_irreducible_size tests GePol cavity irreducible
       * size for a point charge in C2v symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 8;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC2vTest_area tests GePol cavity surface area for a
       * point charge in C2v symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 4.0 * M_PI * pow(1.0, 2);
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC2vTest_volume tests GePol cavity volume for a point
       * charge in C2v symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
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

  GIVEN("A single sphere") {
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    WHEN("the point group is C2h") {
      Molecule point = dummy<6>();
      GePolCavity cavity(point, area, probeRadius, minRadius, "c2h");

      /*! \class GePolCavity
       *  \test \b GePolCavityC2hTest_size tests GePol cavity size for a point charge
       * in C2h symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 32;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC2hTest_irreducible_size tests GePol cavity irreducible
       * size for a point charge in C2h symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 8;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC2hTest_area tests GePol cavity surface area for a
       * point charge in C2h symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 4.0 * M_PI * pow(1.0, 2);
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityC2hTest_volume tests GePol cavity volume for a point
       * charge in C2h symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
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

  GIVEN("A single sphere") {
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    WHEN("the point group is D2h") {
      Molecule point = dummy<7>();
      GePolCavity cavity(point, area, probeRadius, minRadius, "d2h");

      /*! \class GePolCavity
       *  \test \b GePolCavityD2hTest_size tests GePol cavity size for a point charge
       * in D2h symmetry with added spheres
       */
      THEN("the size of the cavity is") {
        int size = 32;
        int actualSize = cavity.size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityD2hTest_irreducible_size tests GePol cavity irreducible
       * size for a point charge in D2h symmetry with added spheres
       */
      AND_THEN("the irreducible size of the cavity is") {
        int size = 4;
        int actualSize = cavity.irreducible_size();
        REQUIRE(size == actualSize);
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityD2hTest_area tests GePol cavity surface area for a
       * point charge in D2h symmetry with added spheres
       */
      AND_THEN("the surface area of the cavity is") {
        double area = 4.0 * M_PI * pow(1.0, 2);
        double actualArea = cavity.elementArea().sum();
        REQUIRE(area == Approx(actualArea));
      }
      /*! \class GePolCavity
       *  \test \b GePolCavityD2hTest_volume tests GePol cavity volume for a point
       * charge in D2h symmetry with added spheres
       */
      AND_THEN("the volume of the cavity is") {
        double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
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
