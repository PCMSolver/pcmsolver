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
#include "cavity/GePolCavity.hpp"
#include "green/DerivativeTypes.hpp"
#include "green/UniformDielectric.hpp"
#include "green/Vacuum.hpp"
#include "solver/IEFSolver.hpp"

using namespace pcm;
using bi_operators::Collocation;
using cavity::GePolCavity;
using green::Vacuum;
using green::UniformDielectric;
using solver::IEFSolver;

/*! \class IEFSolver
 *  \test \b NH3GePolRestart tests IEFSolver using ammonia with a GePol cavity read
 * from .npz file
 */
TEST_CASE("Test solver for the IEFPCM for NH3 and a restarted GePol cavity",
          "[solver][iefpcm][iefpcm_gepol-NH3_from-file]") {
  Eigen::Vector3d N(-0.000000000, -0.104038047, 0.000000000);
  Eigen::Vector3d H1(-0.901584415, 0.481847022, -1.561590016);
  Eigen::Vector3d H2(-0.901584415, 0.481847022, 1.561590016);
  Eigen::Vector3d H3(1.803168833, 0.481847022, 0.000000000);
  // Load  cavity
  GePolCavity cavity;
  cavity.loadCavity("nh3.npz");

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
  Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
  for (int i = 0; i < size; ++i) {
    Eigen::Vector3d center = cavity.elementCenter(i);
    double Ndistance = (center - N).norm();
    double H1distance = (center - H1).norm();
    double H2distance = (center - H2).norm();
    double H3distance = (center - H3).norm();
    fake_mep(i) = Ncharge / Ndistance + Hcharge / H1distance + Hcharge / H2distance +
                  Hcharge / H3distance;
  }
  // The total ASC for a dielectric is -Q*(epsilon-1)/epsilon
  Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
  fake_asc = solver.computeCharge(fake_mep);
  double totalASC = -(Ncharge + 3.0 * Hcharge) * (permittivity - 1) / permittivity;
  double totalFakeASC = fake_asc.sum();
  CAPTURE(totalASC - totalFakeASC);
  REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
}
