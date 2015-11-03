/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013-2015 Roberto Di Remigio, Luca Frediani and contributors
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
 *     PCMSolver API, see: <http://pcmsolver.github.io/pcmsolver-doc>
 */
/* pcmsolver_copyright_end */

#include <catch.hpp>

#include <cmath>
#include <vector>

#include "Config.hpp"

#include <Eigen/Core>

#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "GePolCavity.hpp"
#include "Molecule.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "IEFSolver.hpp"
#include "TestingMolecules.hpp"

/*! \class IEFSolver
 *  \test \b C2H4GePolD2h tests IEFSolver using C2H4 with a GePol cavity in D2h symmetry
 */
TEST_CASE("Test solver for the IEFPCM and the C2H4 molecule in D2h symmetry", "[iefpcm][iefpcm_symmetry][iefpcm_gepol-C2H4_D2h]")
{
  Molecule molec = C2H4();
  double area = 0.2 / convertBohr2ToAngstrom2;
  double probeRadius = 1.385 / convertBohrToAngstrom;
  double minRadius = 100.0 / convertBohrToAngstrom;
  GePolCavity cavity = GePolCavity(molec, area, probeRadius, minRadius, "ief_d2h_noadd");

  double permittivity = 78.39;
  Vacuum<AD_directional, CollocationIntegrator> gfInside = Vacuum<AD_directional, CollocationIntegrator>();
  UniformDielectric<AD_directional, CollocationIntegrator> gfOutside =
    UniformDielectric<AD_directional, CollocationIntegrator>(permittivity);
  bool symm = true;
  IEFSolver solver(symm);
  solver.buildSystemMatrix(cavity, gfInside, gfOutside);

  double Ccharge = 6.0;
  double Hcharge = 1.0;
  size_t size = cavity.size();
  Eigen::VectorXd fake_mep = computeMEP(molec, cavity.elements());
  // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
  size_t irr_size = cavity.irreducible_size();
  Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(irr_size);
  fake_asc = solver.computeCharge(fake_mep);

  for (size_t i = 0; i < size; ++i) {
    INFO("fake_mep(" << i << ") = " << fake_mep(i));
  }
  for (size_t i = 0; i < size; ++i) {
    INFO("fake_asc(" << i << ") = " << fake_asc(i));
  }

  double totalASC = - (2.0 * Ccharge + 4.0 * Hcharge) * (permittivity - 1) / permittivity;
  int nr_irrep = cavity.pointGroup().nrIrrep();
  double totalFakeASC = fake_asc.sum() * nr_irrep;
  CAPTURE(totalASC);
  CAPTURE(totalFakeASC);
  CAPTURE(totalASC - totalFakeASC);
  REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
}
