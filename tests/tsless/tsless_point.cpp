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
#include "cavity/TsLessCavity.hpp"

using namespace pcm;
using cavity::TsLessCavity;

TEST_CASE("TsLess cavity for a single sphere", "[tsless][tsless_point]") {
  Molecule point = dummy<0>();
  double area = 0.4;
  double minDistance = 0.1;
  double probeRadius = 0.0;
  int derOrder = 4;
  TsLessCavity cavity(point, area, probeRadius, 100.0, minDistance, derOrder);

  /*! \class TsLessCavity
   *  \test \b TsLessCavityTest_size tests TsLess cavity size for a point charge
   */
  SECTION("Test size") {
    size_t ref_size = 32;
    size_t size = cavity.size();
    REQUIRE(ref_size == size);
  }

  /*! \class TsLessCavity
   *  \test \b TsLessCavityTest_area tests TsLess cavity surface area for a point
   * charge
   */
  SECTION("Test surface area") {
    double ref_area = 4.0 * M_PI * pow(1.0, 2);
    double area = cavity.elementArea().sum();
    REQUIRE(ref_area == Approx(area));
  }

  /*! \class TsLessCavity
   *  \test \b TsLessCavityTest_volume tests TsLess cavity volume for a point charge
   */
  SECTION("Test volume") {
    double ref_volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double volume = 0;
    for (size_t i = 0; i < cavity.size(); ++i) {
      volume +=
          cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
    }
    volume /= 3;
    REQUIRE(ref_volume == Approx(volume));
  }
}
