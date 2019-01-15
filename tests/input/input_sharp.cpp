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
 *  \test \b InputSharpTest_Sharp tests input reading on an input file parsed by
 * pcmsolver.py
 */
TEST_CASE("Input reading using GetKw for an input file for a sharp environment",
          "[input][input_sharp]") {
  std::string filename = "@sharp.inp";
  Input parsedInput = Input(filename);
  std::string units = "ANGSTROM";
  std::string cavityType = "GEPOL";
  double area = 0.3 * angstrom2ToBohr2();
  bool scaling = true;
  double diagonalScaling = 1.07;
  std::string radiiSet = "BONDI";
  std::string mode = "IMPLICIT";
  std::string solverType = "IEFPCM";
  double probeRadius = 1.5 * angstromToBohr(); // The value for water
  std::string greenInsideType = "VACUUM_DERIVATIVE";
  std::string greenOutsideType = "SPHERICALSHARP_DERIVATIVE";
  double epsilonInside = 1.0;
  double epsilonStatic1 = 114.0;
  double epsilonDynamic1 = 7.3;
  double epsilonStatic2 = 35.7;
  double epsilonDynamic2 = 1.81;
  double center = 100.0 * angstromToBohr();
  Eigen::Vector3d origin = (Eigen::Vector3d() << 70.0, 1.0, 23.0).finished();
  origin *= angstromToBohr();

  REQUIRE(units == parsedInput.units());
  REQUIRE(cavityType == parsedInput.cavityParams().cavityType);
  REQUIRE(area == Approx(parsedInput.cavityParams().area));
  REQUIRE(scaling == parsedInput.scaling());
  REQUIRE(radiiSet == parsedInput.radiiSet());
  REQUIRE(mode == parsedInput.mode());
  REQUIRE(diagonalScaling == Approx(parsedInput.integratorScaling()));
  REQUIRE(solverType == parsedInput.solverParams().solverType);
  REQUIRE(probeRadius == Approx(parsedInput.cavityParams().probeRadius));
  REQUIRE(greenInsideType == parsedInput.insideGreenParams().greensFunctionType);
  REQUIRE(greenOutsideType ==
          parsedInput.outsideStaticGreenParams().greensFunctionType);
  REQUIRE(epsilonInside == Approx(parsedInput.insideGreenParams().epsilon));
  REQUIRE(epsilonStatic1 == Approx(parsedInput.outsideStaticGreenParams().epsilon1));
  REQUIRE(epsilonStatic2 == Approx(parsedInput.outsideStaticGreenParams().epsilon2));
  REQUIRE(epsilonDynamic1 ==
          Approx(parsedInput.outsideDynamicGreenParams().epsilon1));
  REQUIRE(epsilonDynamic2 ==
          Approx(parsedInput.outsideDynamicGreenParams().epsilon2));
  REQUIRE(center == Approx(parsedInput.outsideDynamicGreenParams().center));
  for (int i = 0; i < 3; ++i) {
    REQUIRE(origin(i) == Approx(parsedInput.outsideStaticGreenParams().origin(i)));
    REQUIRE(origin(i) == Approx(parsedInput.outsideDynamicGreenParams().origin(i)));
  }
}
