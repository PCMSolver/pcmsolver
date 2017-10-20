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
 *  \test \b InputDiffuseTest_Diffuse tests input reading on an input file parsed by
 * pcmsolver.py
 */
TEST_CASE("Input reading using GetKw for an input file for a diffuse environment",
          "[input][input_diffuse]") {
  std::string filename = "@diffuse.inp";
  Input parsedInput = Input(filename);
  std::string units = "ANGSTROM";
  int CODATAyear = 1998;
  std::string type = "WAVELET";
  int patchLevel = 1;
  double coarsity = 0.3;
  bool scaling = true;
  double diagonalScaling = 1.07;
  std::string radiiSet = "BONDI";
  std::string mode = "EXPLICIT";
  Eigen::Vector3d c1, c2, c3;
  c1 << 0.00, 0.00, 0.00;
  c2 << 0.00, -0.96, 0.00;
  c3 << -0.905, 0.32, 0.00;
  c2 *= angstromToBohr();
  c3 *= angstromToBohr();
  Sphere sph1(c1, 1.80 * angstromToBohr());
  Sphere sph2(c2, 1.44 * angstromToBohr());
  Sphere sph3(c3, 1.44 * angstromToBohr());
  std::vector<Sphere> spheres;
  spheres.push_back(sph1);
  spheres.push_back(sph2);
  spheres.push_back(sph3);
  std::string solverType = "WAVELET";
  int equationType = 0;
  double probeRadius = 1.385 * angstromToBohr(); // The value for water
  std::string greenInsideType = "VACUUM";
  std::string greenOutsideType = "SPHERICALDIFFUSE";
  int derivativeInsideType = 0;
  int derivativeOutsideType = 0;
  double epsilonInside = 1.0;
  double epsilonStatic1 = 78.39;
  double epsilonDynamic1 = 10.423;
  double epsilonStatic2 = 20.0;
  double epsilonDynamic2 = 4.0;
  double center = 100.0 * angstromToBohr();
  double width = 5.0 * angstromToBohr();
  int profile = 1;
  Eigen::Vector3d origin = (Eigen::Vector3d() << 70.0, 1.0, 23.0).finished();
  origin *= angstromToBohr();

  REQUIRE(units == parsedInput.units());
  REQUIRE(CODATAyear == parsedInput.CODATAyear());
  REQUIRE(type == parsedInput.cavityType());
  REQUIRE(patchLevel == parsedInput.cavityParams().patchLevel);
  REQUIRE(coarsity == Approx(parsedInput.cavityParams().coarsity));
  REQUIRE(scaling == parsedInput.scaling());
  REQUIRE(radiiSet == parsedInput.radiiSet());
  REQUIRE(mode == parsedInput.mode());
  REQUIRE(diagonalScaling == Approx(parsedInput.integratorScaling()));
  for (size_t i = 0; i < spheres.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      REQUIRE(spheres[i].center(j) == Approx(parsedInput.spheres(i).center(j)));
    }
    REQUIRE(spheres[i].radius == Approx(parsedInput.spheres(i).radius));
  }
  REQUIRE(solverType == parsedInput.solverType());
  REQUIRE(equationType == parsedInput.equationType());
  REQUIRE(probeRadius == Approx(parsedInput.cavityParams().probeRadius));
  REQUIRE(greenInsideType == parsedInput.greenInsideType());
  REQUIRE(greenOutsideType == parsedInput.greenOutsideType());
  REQUIRE(derivativeInsideType == parsedInput.insideGreenParams().howDerivative);
  REQUIRE(derivativeOutsideType ==
          parsedInput.outsideStaticGreenParams().howDerivative);
  REQUIRE(epsilonInside == Approx(parsedInput.insideGreenParams().epsilon));
  REQUIRE(epsilonStatic1 == Approx(parsedInput.outsideStaticGreenParams().epsilon1));
  REQUIRE(epsilonStatic2 == Approx(parsedInput.outsideStaticGreenParams().epsilon2));
  REQUIRE(epsilonDynamic1 ==
          Approx(parsedInput.outsideDynamicGreenParams().epsilon1));
  REQUIRE(epsilonDynamic2 ==
          Approx(parsedInput.outsideDynamicGreenParams().epsilon2));
  REQUIRE(center == Approx(parsedInput.outsideDynamicGreenParams().center));
  REQUIRE(width == Approx(parsedInput.outsideDynamicGreenParams().width));
  REQUIRE(profile == parsedInput.outsideDynamicGreenParams().howProfile);
  for (int i = 0; i < 3; ++i) {
    REQUIRE(origin(i) == Approx(parsedInput.outsideStaticGreenParams().origin(i)));
    REQUIRE(origin(i) == Approx(parsedInput.outsideDynamicGreenParams().origin(i)));
  }
}
