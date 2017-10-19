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

#include <catch.hpp>

#include <cmath>
#include <vector>

#include <Eigen/Core>

#include "TestingMolecules.hpp"
#include "bi_operators/Collocation.hpp"
#include "cavity/GePolCavity.hpp"
#include "green/DerivativeTypes.hpp"
#include "green/UniformDielectric.hpp"
#include "green/Vacuum.hpp"
#include "solver/CPCMSolver.hpp"
#include "utils/Molecule.hpp"

using namespace pcm;
using bi_operators::Collocation;
using cavity::GePolCavity;
using green::Vacuum;
using green::UniformDielectric;
using solver::CPCMSolver;

/*! \class CPCMSolver
 *  \test \b C2H4GePolD2h tests CPCMSolver using C2H4 with a GePol cavity in D2h
 * symmetry
 */
TEST_CASE("Test solver for the CPCM and the C2H4 molecule in D2h symmetry",
          "[cpcm][cpcm_symmetry][cpcm_gepol-C2H4_D2h]") {
  Molecule molec = C2H4();
  double area = 0.2 / bohr2ToAngstrom2();
  double probeRadius = 1.385 / bohrToAngstrom();
  double minRadius = 100.0 / bohrToAngstrom();
  GePolCavity cavity(molec, area, probeRadius, minRadius, "cpcm_d2h_noadd");

  double permittivity = 78.39;
  Vacuum<> gf_i;
  UniformDielectric<> gf_o(permittivity);
  bool symm = true;
  double correction = 0.0;

  Collocation S;

  CPCMSolver solver(symm, correction);
  solver.buildSystemMatrix(cavity, gf_i, gf_o, S);

  double Ccharge = 6.0;
  double Hcharge = 1.0;
  int size = cavity.size();
  Eigen::VectorXd fake_mep = computeMEP(molec, cavity.elements());
  // The total ASC for a dielectric is -Q*(epsilon-1)/epsilon
  int irr_size = cavity.irreducible_size();
  Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
  fake_asc = solver.computeCharge(fake_mep);

  for (int i = 0; i < size; ++i) {
    INFO("fake_mep(" << i << ") = " << fake_mep(i));
  }
  for (int i = 0; i < size; ++i) {
    INFO("fake_asc(" << i << ") = " << fake_asc(i));
  }

  // The total ASC for a conductor is -Q
  // for CPCM it will be -Q*(epsilon-1)/epsilon
  double totalASC =
      -(2.0 * Ccharge + 4.0 * Hcharge) * (permittivity - 1) / permittivity;
  // Renormalize
  int nr_irrep = cavity.pointGroup().nrIrrep();
  double totalFakeASC = fake_asc.sum() * nr_irrep;
  CAPTURE(totalASC);
  CAPTURE(totalFakeASC);
  CAPTURE(totalASC - totalFakeASC);
  REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));

  Eigen::VectorXd reference =
      cnpy::custom::npy_load<double>("ASC-cpcm_gepol-C2H4_D2h.npy");
  for (int i = 0; i < cavity.size(); ++i) {
    REQUIRE(reference(i) == Approx(fake_asc(i)));
  }
}
