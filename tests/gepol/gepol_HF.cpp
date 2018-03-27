/*
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

#include <Eigen/Core>

#include "TestingMolecules.hpp"
#include "cavity/GePolCavity.hpp"
#include "utils/Molecule.hpp"
#include "utils/Symmetry.hpp"

using namespace pcm;
using cavity::GePolCavity;

TEST_CASE("GePol cavity for an hydrogen fluoride molecule", "[gepol][gepol_HF]") {
  Molecule molec = HF();
  double area = 0.02 / bohr2ToAngstrom2();
  double probeRadius = 0.0 / bohrToAngstrom();
  double minRadius = 0.2 / bohrToAngstrom();
  GePolCavity cavity(molec, area, probeRadius, minRadius);

  /*! \class GePolCavity
   *  \test \b GePolCavityHFTest_size tests GePol cavity size for ammonia
   */
  SECTION("Test size") {
    size_t ref_size = 1688;
    size_t size = cavity.size();
    REQUIRE(size == ref_size);
  }

  /*! \class GePolCavity
   *  \test \b GePolCavityHFTest_area tests GePol cavity surface area for ammonia
   */
  SECTION("Test surface area") {
    double ref_area = 110.64517236179323;
    double area = cavity.elementArea().sum();
    REQUIRE(area == Approx(ref_area));
  }

  /*! \class GePolCavity
   *  \test \b GePolCavityHFTest_volume tests GePol cavity volume for ammonia
   */
  SECTION("Test volume") {
    double ref_volume = 105.98002389543836;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double volume = 0;
    for (size_t i = 0; i < cavity.size(); ++i) {
      volume +=
          cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
    }
    volume /= 3;
    REQUIRE(volume == Approx(ref_volume));
  }
}
