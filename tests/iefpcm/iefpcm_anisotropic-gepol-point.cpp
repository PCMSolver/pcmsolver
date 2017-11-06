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
#include "solver/IEFSolver.hpp"
#include "utils/Molecule.hpp"

using namespace pcm;
using bi_operators::Collocation;
using cavity::GePolCavity;
using green::Vacuum;
using green::UniformDielectric;
using solver::IEFSolver;

SCENARIO(
    "Test solver for the anisotropic IEFPCM for a point charge and a GePol cavity",
    "[solver][iefpcm][iefpcm_anisotropic-gepol-point][anisotropic]") {
  GIVEN("An isotropic environment modelled and a point charge forcing the use of an "
        "anisotropic solver") {
    double permittivity = 78.39;
    Vacuum<> gf_i;
    UniformDielectric<> gf_o(permittivity);
    bool symm = true;

    Collocation op;

    double charge = 8.0;
    double totalASC = -charge * (permittivity - 1) / permittivity;

    /*! \class IEFSolver
     *  \test \b anisotropicPointChargeGePol tests IEFSolver using a point charge
     * with a GePol cavity
     *  We are forcing the usage of the buildAnisotropicMatrix method.
     *  The results are compared with Gauss' theorem and the results from the
     * buildIsotropicMatrix method
     *  The point charge is at the origin.
     */
    WHEN("the point charge is located at the origin") {
      Molecule point = dummy<0>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity = GePolCavity(point, area, probeRadius, minRadius);

      IEFSolver aniso_solver(symm);
      aniso_solver.buildAnisotropicMatrix(cavity, gf_i, gf_o, op);

      IEFSolver iso_solver(symm);
      iso_solver.buildIsotropicMatrix(cavity, gf_i, gf_o, op);

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);
      for (int i = 0; i < size; ++i) {
        INFO("fake_mep(" << i << ") = " << fake_mep(i));
      }

      THEN("the apparent surface charge is") {
        Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
        Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);
        double totalAnisoASC = aniso_fake_asc.sum();
        double totalIsoASC = iso_fake_asc.sum();

        for (int i = 0; i < size; ++i) {
          INFO("aniso_fake_asc(" << i << ") = " << aniso_fake_asc(i));
        }
        for (int i = 0; i < size; ++i) {
          INFO("iso_fake_asc(" << i << ") = " << iso_fake_asc(i));
        }

        CAPTURE(totalASC);
        CAPTURE(totalAnisoASC);
        CAPTURE(totalASC - totalAnisoASC);
        REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
        REQUIRE(totalIsoASC == Approx(totalAnisoASC));
      }
    }

    /*! \class IEFSolver
     *  \test \b anisotropicPointChargeShiftedGePol tests IEFSolver using a point
     * charge with a GePol cavity
     *  We are forcing the usage of the buildAnisotropicMatrix method.
     *  The results are compared with Gauss' theorem and the results from the
     * buildIsotropicMatrix method
     *  The point charge is away from the origin.
     */
    AND_WHEN("the point charge is located away from the origin") {
      Eigen::Vector3d origin = 100 * Eigen::Vector3d::Random();
      Molecule point = dummy<0>(2.929075493, origin);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity = GePolCavity(point, area, probeRadius, minRadius);

      IEFSolver aniso_solver(symm);
      aniso_solver.buildAnisotropicMatrix(cavity, gf_i, gf_o, op);

      IEFSolver iso_solver(symm);
      iso_solver.buildIsotropicMatrix(cavity, gf_i, gf_o, op);

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge, origin);
      for (int i = 0; i < size; ++i) {
        INFO("fake_mep(" << i << ") = " << fake_mep(i));
      }
      THEN("the apparent surface charge is") {
        Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
        Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

        for (int i = 0; i < size; ++i) {
          INFO("aniso_fake_asc(" << i << ") = " << aniso_fake_asc(i));
        }
        for (int i = 0; i < size; ++i) {
          INFO("iso_fake_asc(" << i << ") = " << iso_fake_asc(i));
        }

        double totalAnisoASC = aniso_fake_asc.sum();
        double totalIsoASC = iso_fake_asc.sum();
        CAPTURE(totalASC);
        CAPTURE(totalAnisoASC);
        CAPTURE(totalASC - totalAnisoASC);
        REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
        REQUIRE(totalIsoASC == Approx(totalAnisoASC));
      }
    }
  }
}
