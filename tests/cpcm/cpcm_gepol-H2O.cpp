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

#include <Eigen/Core>

#include "TestingMolecules.hpp"
#include "bi_operators/Collocation.hpp"
#include "cavity/GePolCavity.hpp"
#include "green/DerivativeTypes.hpp"
#include "green/UniformDielectric.hpp"
#include "green/Vacuum.hpp"
#include "solver/CPCMSolver.hpp"
#include "utils/Molecule.hpp"
#include "utils/Symmetry.hpp"

using namespace pcm;
using bi_operators::Collocation;
using cavity::GePolCavity;
using green::Vacuum;
using green::UniformDielectric;
using solver::CPCMSolver;

/*! \class CPCMSolver
 *  \test \b NH3GePol tests CPCMSolver using water and a GePol cavity
 */
TEST_CASE("Test solver for the C-PCM with H2O molecule and a GePol cavity",
          "[solver][cpcm][cpcm_gepol-H2O]") {
  Molecule molecule = H2O();
  double area = 0.2 / bohr2ToAngstrom2();
  double probeRadius = 0.0;
  double minRadius = 100.0;
  GePolCavity cavity = GePolCavity(molecule, area, probeRadius, minRadius);

  double permittivity = 78.39;
  Vacuum<> gf_i;
  UniformDielectric<> gf_o(permittivity);
  bool symm = true;
  double correction = 0.0;

  CPCMSolver solver(symm, correction);
  solver.buildSystemMatrix(cavity, gf_i, gf_o, Collocation());

  double Ocharge = 8.0;
  double Hcharge = 1.0;
  int size = cavity.size();
  Eigen::VectorXd fake_mep = computeMEP(molecule, cavity.elements());
  // The total ASC for a conductor is -Q
  // for CPCM it will be -Q*(epsilon-1)/(epsilon + correction)
  Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
  fake_asc = solver.computeCharge(fake_mep);
  double totalASC =
      -(Ocharge + 2.0 * Hcharge) * (permittivity - 1) / (permittivity + correction);
  double totalFakeASC = fake_asc.sum();
  CAPTURE(totalASC - totalFakeASC);
  REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));

  Eigen::VectorXd reference =
      cnpy::custom::npy_load<double>("ASC-cpcm_gepol-H2O.npy");
  for (int i = 0; i < cavity.size(); ++i) {
    REQUIRE(reference(i) == Approx(fake_asc(i)));
  }
}
