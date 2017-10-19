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
 *  \test \b Input_TsLess tests input reading on an input file parsed by pcmsolver.py
 */
TEST_CASE("Input reading using GetKw for an input file for a TsLess cavity",
          "[input][input_tsless]") {
  std::string filename = "@tsless.inp";
  Input parsedInput = Input(filename);
  std::string units = "ANGSTROM";
  int CODATAyear = 2002;
  std::string type = "TSLESS";
  double area = 0.6 * angstrom2ToBohr2();
  double minDistance = 5.0 * angstromToBohr();
  int derOrder = 25;
  bool scaling = false;
  double diagonalScaling = 1.0694;
  std::string radiiSet = "UFF";
  double minimalRadius = 0.1 * angstromToBohr();
  std::string mode = "ATOMS";
  std::vector<int> atoms;
  std::vector<double> radii;
  atoms.push_back(1);
  atoms.push_back(2);
  atoms.push_back(4);
  radii.push_back(1.4 * angstromToBohr());
  radii.push_back(4.3 * angstromToBohr());
  radii.push_back(1.2 * angstromToBohr());
  std::string solvent = "Cyclohexane"; // Name in the Solvent object
  std::string solverType = "CPCM";
  double correction = 0.5;
  bool hermitivitize = false;
  double probeRadius = 2.815 * angstromToBohr(); // The value for water
  std::string greenInsideType = "VACUUM";
  std::string greenOutsideType = "UNIFORMDIELECTRIC";
  int derivativeInsideType = 1;
  int derivativeOutsideType = 1;

  REQUIRE(units == parsedInput.units());
  REQUIRE(CODATAyear == parsedInput.CODATAyear());
  REQUIRE(type == parsedInput.cavityType());
  REQUIRE(area == Approx(parsedInput.cavityParams().area));
  REQUIRE(minDistance == Approx(parsedInput.cavityParams().minDistance));
  REQUIRE(derOrder == parsedInput.cavityParams().derOrder);
  REQUIRE(scaling == parsedInput.scaling());
  REQUIRE(diagonalScaling == Approx(parsedInput.integratorScaling()));
  REQUIRE(radiiSet == parsedInput.radiiSet());
  REQUIRE(minimalRadius == Approx(parsedInput.cavityParams().minimalRadius));
  REQUIRE(mode == parsedInput.mode());
  for (size_t i = 0; i < atoms.size(); ++i) {
    REQUIRE(atoms[i] == parsedInput.atoms(i));
    REQUIRE(radii[i] == Approx(parsedInput.radii(i)));
  }
  REQUIRE(solvent == parsedInput.solvent().name);
  REQUIRE(solverType == parsedInput.solverType());
  REQUIRE(correction == Approx(parsedInput.correction()));
  REQUIRE(hermitivitize == parsedInput.hermitivitize());
  REQUIRE(probeRadius == Approx(parsedInput.cavityParams().probeRadius));
  REQUIRE(greenInsideType == parsedInput.greenInsideType());
  REQUIRE(greenOutsideType == parsedInput.greenOutsideType());
  REQUIRE(derivativeInsideType == parsedInput.insideGreenParams().howDerivative);
  REQUIRE(derivativeOutsideType ==
          parsedInput.outsideStaticGreenParams().howDerivative);
}
