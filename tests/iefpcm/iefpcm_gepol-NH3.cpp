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

#include "bi_operators/Collocation.hpp"
#include "green/DerivativeTypes.hpp"
#include "cavity/GePolCavity.hpp"
#include "utils/Molecule.hpp"
#include "green/Vacuum.hpp"
#include "TestingMolecules.hpp"
#include "green/UniformDielectric.hpp"
#include "solver/IEFSolver.hpp"
#include "utils/Symmetry.hpp"

using namespace pcm;
using bi_operators::Collocation;
using cavity::GePolCavity;
using green::Vacuum;
using green::UniformDielectric;
using solver::IEFSolver;

/*! \class IEFSolver
 *  \test \b NH3GePol tests IEFSolver using ammonia and a GePol cavity
 */
TEST_CASE("Test solver for the IEFPCM with NH3 molecule and a GePol cavity",
          "[solver][iefpcm][iefpcm_gepol-NH3]") {
  Molecule molec = NH3();

  double area = 0.4;
  double probeRadius = 0.0;
  double minRadius = 100.0;
  GePolCavity cavity = GePolCavity(molec, area, probeRadius, minRadius);
  cavity.saveCavity("nh3.npz");

  double permittivity = 78.39;
  Vacuum<> gf_i;
  UniformDielectric<> gf_o(permittivity);

  Collocation op;

  bool symm = true;
  IEFSolver solver(symm);
  solver.buildSystemMatrix(cavity, gf_i, gf_o, op);

  double Ncharge = 7.0;
  double Hcharge = 1.0;
  int size = cavity.size();
  Eigen::VectorXd fake_mep = computeMEP(molec, cavity.elements());
  // The total ASC for a dielectric is -Q*(epsilon-1)/epsilon
  Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
  fake_asc = solver.computeCharge(fake_mep);
  double totalASC = -(Ncharge + 3.0 * Hcharge) * (permittivity - 1) / permittivity;
  double totalFakeASC = fake_asc.sum();
  CAPTURE(totalASC - totalFakeASC);
  REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
}
