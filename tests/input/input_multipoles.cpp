/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2019 Roberto Di Remigio, Luca Frediani and contributors.
 *
 *     This file is part of PCMSolver.
 *
 *     PCMSolver is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     PCMSolver is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
 *
 *     For information on the complete list of contributors to the
 *     PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
 */
/* pcmsolver_copyright_end */

#include "catch.hpp"

#include <iostream>
#include <string>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include "bi_operators/BIOperatorData.hpp"
#include "cavity/CavityData.hpp"
#include "green/GreenData.hpp"
#include "interface/Input.hpp"
#include "solver/SolverData.hpp"
#include "utils/Sphere.hpp"

using pcm::Input;
using pcm::utils::Sphere;

/*! \class Input
 *  \test \b InputMultipolesTest_Multipoles tests input reading on an input file
 * parsed
 * by go_pcm.py
 */
TEST_CASE("Input reading using GetKw for an input file for a classical point "
          "multipoles distribution",
          "[input][multipoles][input_multipoles]") {
  Input parsedInput("@multipoles.inp");
  std::string units = "ANGSTROM";
  int CODATAyear = 2002;
  std::string cavityType = "GEPOL";
  std::string radiiSet = "BONDI";
  std::string mode = "IMPLICIT";
  std::string solverType = "IEFPCM";
  std::string greenInsideType = "VACUUM_DERIVATIVE";
  std::string greenOutsideType = "UNIFORMDIELECTRIC_DERIVATIVE";
  double area = 10.0 * angstrom2ToBohr2();
  std::string solvent = "Acetonitrile"; // Name in the Solvent object
  REQUIRE(units == parsedInput.units());
  REQUIRE(CODATAyear == parsedInput.CODATAyear());
  REQUIRE(cavityType == parsedInput.cavityParams().cavityType);
  REQUIRE(radiiSet == parsedInput.radiiSet());
  REQUIRE(mode == parsedInput.mode());
  REQUIRE(solverType == parsedInput.solverParams().solverType);
  REQUIRE(greenInsideType == parsedInput.insideGreenParams().greensFunctionType);
  REQUIRE(greenOutsideType ==
          parsedInput.outsideStaticGreenParams().greensFunctionType);
  REQUIRE(area == Approx(parsedInput.cavityParams().area));
  REQUIRE(solvent == parsedInput.solvent().name);
  for (size_t j = 0; j < 3; ++j) {
    REQUIRE(0.0 == Approx(parsedInput.multipoles().monopolesSites(j, 0)));
  }
  REQUIRE(1.0 == Approx(parsedInput.multipoles().monopoles(0)));
  for (size_t j = 0; j < 3; ++j) {
    REQUIRE(0.0 == Approx(parsedInput.multipoles().dipolesSites(j, 0)));
    REQUIRE(Eigen::Vector3d::UnitZ()(j) ==
            Approx(parsedInput.multipoles().dipoles(j, 0)));
  }
}
