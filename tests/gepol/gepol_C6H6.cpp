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

TEST_CASE("GePol cavity for the benzene molecule", "[gepol][gepol_C6H6]") {
  double area = 0.3 / bohr2ToAngstrom2();
  double probeRadius = 1.385 / bohrToAngstrom();
  // Addition of spheres is enabled, but will not happen in this particular case
  double minRadius = 10.0 / bohrToAngstrom();
  Molecule molec = C6H6();
  GePolCavity cavity(molec, area, probeRadius, minRadius, "c6h6");
  cavity.saveCavity("c6h6.npz");

  /*! \class GePolCavity
   *  \test \b GePolCavityC6H6AddTest_size tests GePol cavity size for C6H6 in C1
   * symmetry with added spheres
   */
  SECTION("Test size") {
    int size = 644;
    int actualSize = cavity.size();
    REQUIRE(size == actualSize);
  }

  /*! \class GePolCavity
   *  \test \b GePolCavityC6H6AddTest_area tests GePol cavity surface area for C6H6
   * in C1 symmetry with added spheres
   */
  SECTION("Test surface area") {
    double area = 391.06094362589062;
    double actualArea = cavity.elementArea().sum();
    REQUIRE(area == Approx(actualArea));
  }

  /*! \class GePolCavity
   *  \test \b GePolCavityC6H6AddTest_volume tests GePol cavity volume for C6H6 in C1
   * symmetry with added spheres
   */
  SECTION("Test volume") {
    double volume = 567.67444339622284;
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
