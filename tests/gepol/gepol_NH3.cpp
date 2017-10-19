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

#include "LoggerInterface.hpp"
#include "TestingMolecules.hpp"
#include "cavity/GePolCavity.hpp"
#include "utils/Molecule.hpp"
#include "utils/Symmetry.hpp"

using namespace pcm;
using cavity::GePolCavity;

TEST_CASE("GePol cavity for an ammonia molecule", "[gepol][gepol_NH3]") {
  Molecule molec = NH3();
  double area = 0.3 / bohr2ToAngstrom2();
  double probeRadius = 1.385 / bohrToAngstrom();
  double minRadius = 0.2 / bohrToAngstrom();
  GePolCavity cavity(molec, area, probeRadius, minRadius, "nh3");
  cavity.saveCavity("nh3.npz");

  /*! \class GePolCavity
   *  \test \b GePolCavityNH3Test_size tests GePol cavity size for ammonia
   */
  SECTION("Test size") {
    int size = 230;
    int actualSize = cavity.size();
    REQUIRE(size == actualSize);
  }

  /*! \class GePolCavity
   *  \test \b GePolCavityNH3Test_area tests GePol cavity surface area for ammonia
   */
  SECTION("Test surface area") {
    double area = 147.13247859942391;
    double actualArea = cavity.elementArea().sum();
    REQUIRE(area == Approx(actualArea));
  }

  /*! \class GePolCavity
   *  \test \b GePolCavityNH3Test_volume tests GePol cavity volume for ammonia
   */
  SECTION("Test volume") {
    double volume = 153.12929788519045;
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
