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
 * pcmsolver.py
 */
TEST_CASE("Input reading using GetKw for an input file for a restart cavity",
          "[input][input_restart]") {
  std::string filename = "@restart.inp";
  Input parsedInput = Input(filename);
  std::string units = "AU";
  int CODATAyear = 2010;
  std::string type = "RESTART";
  std::string cavFilename = "cavity.npz";
  int patchLevel = 2;
  double coarsity = 0.5;
  double area = 0.3;
  double minDistance = 0.1;
  int derOrder = 4;
  bool scaling = true;
  std::string radiiSet = "BONDI";
  double minimalRadius = 100.0;
  double diagonalScaling = 1.07;
  std::string mode = "IMPLICIT";
  std::string solvent = "Water"; // Name in the Solvent object
  std::string solverType = "IEFPCM";
  int equationType = 1;
  double correction = 0.0;
  bool hermitivitize = true;
  double probeRadius = 1.385 * angstromToBohr(); // The value for water
  std::string greenInsideType = "VACUUM";
  std::string greenOutsideType = "UNIFORMDIELECTRIC";
  int derivativeInsideType = 1;
  int derivativeOutsideType = 1;

  REQUIRE(units == parsedInput.units());
  REQUIRE(CODATAyear == parsedInput.CODATAyear());
  REQUIRE(type == parsedInput.cavityType());
  REQUIRE(cavFilename == parsedInput.cavityParams().filename);
  REQUIRE(patchLevel == parsedInput.cavityParams().patchLevel);
  REQUIRE(coarsity == Approx(parsedInput.cavityParams().coarsity));
  REQUIRE(area == Approx(parsedInput.cavityParams().area));
  REQUIRE(minDistance == Approx(parsedInput.cavityParams().minDistance));
  REQUIRE(diagonalScaling == Approx(parsedInput.integratorScaling()));
  REQUIRE(derOrder == parsedInput.cavityParams().derOrder);
  REQUIRE(scaling == parsedInput.scaling());
  REQUIRE(radiiSet == parsedInput.radiiSet());
  REQUIRE(minimalRadius == Approx(parsedInput.cavityParams().minimalRadius));
  REQUIRE(mode == parsedInput.mode());
  REQUIRE(solvent == parsedInput.solvent().name);
  REQUIRE(solverType == parsedInput.solverType());
  REQUIRE(equationType == parsedInput.equationType());
  REQUIRE(correction == Approx(parsedInput.correction()));
  REQUIRE(hermitivitize == parsedInput.hermitivitize());
  REQUIRE(probeRadius == Approx(parsedInput.cavityParams().probeRadius));
  REQUIRE(greenInsideType == parsedInput.greenInsideType());
  REQUIRE(greenOutsideType == parsedInput.greenOutsideType());
  REQUIRE(derivativeInsideType == parsedInput.insideGreenParams().howDerivative);
  REQUIRE(derivativeOutsideType ==
          parsedInput.outsideStaticGreenParams().howDerivative);
}
