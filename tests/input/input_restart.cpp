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
#include <string>
#include <vector>

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
 *  \test \b InputRestartTest_Restart tests input reading on an input file parsed by
 * go_pcm.py
 */
TEST_CASE("Input reading using GetKw for an input file for a restart cavity",
          "[input][input_restart]") {
  std::string filename = "@restart.inp";
  Input parsedInput = Input(filename);
  std::string units = "AU";
  int CODATAyear = 2010;
  std::string cavityType = "RESTART";
  std::string cavFilename = "cavity.npz";
  double area = 0.3;
  bool scaling = true;
  std::string radiiSet = "BONDI";
  double minimalRadius = 100.0;
  double diagonalScaling = 1.07;
  std::string mode = "IMPLICIT";
  std::string solvent = "Water"; // Name in the Solvent object
  std::string solverType = "IEFPCM";
  double correction = 0.0;
  bool hermitivitize = true;
  double probeRadius = 1.385 * angstromToBohr(); // The value for water
  std::string greenInsideType = "VACUUM_DERIVATIVE";
  std::string greenOutsideType = "UNIFORMDIELECTRIC_DERIVATIVE";

  REQUIRE(units == parsedInput.units());
  REQUIRE(CODATAyear == parsedInput.CODATAyear());
  REQUIRE(cavityType == parsedInput.cavityParams().cavityType);
  REQUIRE(cavFilename == parsedInput.cavityParams().filename);
  REQUIRE(area == Approx(parsedInput.cavityParams().area));
  REQUIRE(diagonalScaling == Approx(parsedInput.integratorScaling()));
  REQUIRE(scaling == parsedInput.scaling());
  REQUIRE(radiiSet == parsedInput.radiiSet());
  REQUIRE(minimalRadius == Approx(parsedInput.cavityParams().minimalRadius));
  REQUIRE(mode == parsedInput.mode());
  REQUIRE(solvent == parsedInput.solvent().name);
  REQUIRE(solverType == parsedInput.solverParams().solverType);
  REQUIRE(correction == Approx(parsedInput.correction()));
  REQUIRE(hermitivitize == parsedInput.hermitivitize());
  REQUIRE(probeRadius == Approx(parsedInput.cavityParams().probeRadius));
  REQUIRE(greenInsideType == parsedInput.insideGreenParams().greensFunctionType);
  REQUIRE(greenOutsideType ==
          parsedInput.outsideStaticGreenParams().greensFunctionType);
}
