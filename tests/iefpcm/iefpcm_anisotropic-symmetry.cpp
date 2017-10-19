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

using namespace pcm;
using bi_operators::Collocation;
using cavity::GePolCavity;
using green::Vacuum;
using green::UniformDielectric;
using solver::IEFSolver;

SCENARIO("Test solver for the anisotropic IEFPCM for a point charge in different "
         "Abelian point groups",
         "[solver][iefpcm][iefpcm_anisotropic-symmetry][anisotropic]") {
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
     *  \test \b pointChargeGePolC1 tests IEFSolver using a point charge with a GePol
     * cavity in C1 symmetry
     *  We are forcing the usage of the buildAnisotropicMatrix method.
     *  The results are compared with Gauss' theorem and the results from the
     * buildIsotropicMatrix method
     *  The point charge is at the origin.
     */
    WHEN("the point group is C1") {
      Molecule point = dummy<0>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius, "C1");

      IEFSolver aniso_solver(symm);
      aniso_solver.buildAnisotropicMatrix(cavity, gf_i, gf_o, op);

      IEFSolver iso_solver(symm);
      iso_solver.buildIsotropicMatrix(cavity, gf_i, gf_o, op);

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

      THEN("the total apparent surface charge is") {
        Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
        Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

        int nr_irrep = cavity.pointGroup().nrIrrep();
        double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
        double totalIsoASC = iso_fake_asc.sum() * nr_irrep;

        CAPTURE(totalASC);
        CAPTURE(totalAnisoASC);
        CAPTURE(totalASC - totalAnisoASC);
        REQUIRE(totalIsoASC == Approx(totalAnisoASC));
        REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
      }
    }

    /*! \class IEFSolver
     *  \test \b pointChargeGePolC2 tests IEFSolver using a point charge with a GePol
     * cavity in C2 symmetry
     *  We are forcing the usage of the buildAnisotropicMatrix method.
     *  The results are compared with Gauss' theorem and the results from the
     * buildIsotropicMatrix method
     *  The point charge is at the origin.
     */
    WHEN("the point group is C2") {
      Molecule point = dummy<1>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius, "C2");

      IEFSolver aniso_solver(symm);
      aniso_solver.buildAnisotropicMatrix(cavity, gf_i, gf_o, op);

      IEFSolver iso_solver(symm);
      iso_solver.buildIsotropicMatrix(cavity, gf_i, gf_o, op);

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

      THEN("the total apparent surface charge is") {
        Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
        Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

        int nr_irrep = cavity.pointGroup().nrIrrep();
        double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
        double totalIsoASC = iso_fake_asc.sum() * nr_irrep;

        CAPTURE(totalASC);
        CAPTURE(totalAnisoASC);
        CAPTURE(totalASC - totalAnisoASC);
        REQUIRE(totalIsoASC == Approx(totalAnisoASC));
        REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
      }
    }

    /*! \class IEFSolver
     *  \test \b pointChargeGePolCs tests IEFSolver using a point charge with a GePol
     * cavity in Cs symmetry
     *  We are forcing the usage of the buildAnisotropicMatrix method.
     *  The results are compared with Gauss' theorem and the results from the
     * buildIsotropicMatrix method
     *  The point charge is at the origin.
     */
    WHEN("the point group is Cs") {
      Molecule point = dummy<2>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius, "Cs");

      IEFSolver aniso_solver(symm);
      aniso_solver.buildAnisotropicMatrix(cavity, gf_i, gf_o, op);

      IEFSolver iso_solver(symm);
      iso_solver.buildIsotropicMatrix(cavity, gf_i, gf_o, op);

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

      THEN("the total apparent surface charge is") {
        Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
        Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

        int nr_irrep = cavity.pointGroup().nrIrrep();
        double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
        double totalIsoASC = iso_fake_asc.sum() * nr_irrep;

        CAPTURE(totalASC);
        CAPTURE(totalAnisoASC);
        CAPTURE(totalASC - totalAnisoASC);
        REQUIRE(totalIsoASC == Approx(totalAnisoASC));
        REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
      }
    }

    /*! \class IEFSolver
     *  \test \b pointChargeGePolCi tests IEFSolver using a point charge with a GePol
     * cavity in Ci symmetry
     *  We are forcing the usage of the buildAnisotropicMatrix method.
     *  The results are compared with Gauss' theorem and the results from the
     * buildIsotropicMatrix method
     *  The point charge is at the origin.
     */
    WHEN("the point group is Ci") {
      Molecule point = dummy<3>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius, "Ci");

      IEFSolver aniso_solver(symm);
      aniso_solver.buildAnisotropicMatrix(cavity, gf_i, gf_o, op);

      IEFSolver iso_solver(symm);
      iso_solver.buildIsotropicMatrix(cavity, gf_i, gf_o, op);

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

      THEN("the total apparent surface charge is") {
        Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
        Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

        int nr_irrep = cavity.pointGroup().nrIrrep();
        double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
        double totalIsoASC = iso_fake_asc.sum() * nr_irrep;

        CAPTURE(totalASC);
        CAPTURE(totalAnisoASC);
        CAPTURE(totalASC - totalAnisoASC);
        REQUIRE(totalIsoASC == Approx(totalAnisoASC));
        REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
      }
    }

    /*! \class IEFSolver
     *  \test \b pointChargeGePolD2 tests IEFSolver using a point charge with a GePol
     * cavity in D2 symmetry
     */
    WHEN("the point group is D2") {
      Molecule point = dummy<4>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius, "D2");

      IEFSolver aniso_solver(symm);
      aniso_solver.buildAnisotropicMatrix(cavity, gf_i, gf_o, op);

      IEFSolver iso_solver(symm);
      iso_solver.buildIsotropicMatrix(cavity, gf_i, gf_o, op);

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

      THEN("the total apparent surface charge is") {
        Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
        Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

        int nr_irrep = cavity.pointGroup().nrIrrep();
        double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
        double totalIsoASC = iso_fake_asc.sum() * nr_irrep;

        CAPTURE(totalASC);
        CAPTURE(totalAnisoASC);
        CAPTURE(totalASC - totalAnisoASC);
        REQUIRE(totalIsoASC == Approx(totalAnisoASC));
        REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
      }
    }

    /*! \class IEFSolver
     *  \test \b pointChargeGePolC2v tests IEFSolver using a point charge with a
     * GePol cavity in C2v symmetry
     *  We are forcing the usage of the buildAnisotropicMatrix method.
     *  The results are compared with Gauss' theorem and the results from the
     * buildIsotropicMatrix method
     *  The point charge is at the origin.
     */
    WHEN("the point group is C2v") {
      Molecule point = dummy<5>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius, "C2v");

      IEFSolver aniso_solver(symm);
      aniso_solver.buildAnisotropicMatrix(cavity, gf_i, gf_o, op);

      IEFSolver iso_solver(symm);
      iso_solver.buildIsotropicMatrix(cavity, gf_i, gf_o, op);

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

      THEN("the total apparent surface charge is") {
        Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
        Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

        int nr_irrep = cavity.pointGroup().nrIrrep();
        double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
        double totalIsoASC = iso_fake_asc.sum() * nr_irrep;

        CAPTURE(totalASC);
        CAPTURE(totalAnisoASC);
        CAPTURE(totalASC - totalAnisoASC);
        REQUIRE(totalIsoASC == Approx(totalAnisoASC));
        REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
      }
    }

    /*! \class IEFSolver
     *  \test \b pointChargeGePolC2h tests IEFSolver using a point charge with a
     * GePol cavity in C2h symmetry
     *  We are forcing the usage of the buildAnisotropicMatrix method.
     *  The results are compared with Gauss' theorem and the results from the
     * buildIsotropicMatrix method
     *  The point charge is at the origin.
     */
    WHEN("the point group is C2h") {
      Molecule point = dummy<6>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius, "C2h");

      IEFSolver aniso_solver(symm);
      aniso_solver.buildAnisotropicMatrix(cavity, gf_i, gf_o, op);

      IEFSolver iso_solver(symm);
      iso_solver.buildIsotropicMatrix(cavity, gf_i, gf_o, op);

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

      THEN("the total apparent surface charge is") {
        Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
        Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

        int nr_irrep = cavity.pointGroup().nrIrrep();
        double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
        double totalIsoASC = iso_fake_asc.sum() * nr_irrep;

        CAPTURE(totalASC);
        CAPTURE(totalAnisoASC);
        CAPTURE(totalASC - totalAnisoASC);
        REQUIRE(totalIsoASC == Approx(totalAnisoASC));
        REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
      }
    }

    /*! \class IEFSolver
     *  \test \b pointChargeGePolD2h tests IEFSolver using a point charge with a
     * GePol cavity in D2h symmetry
     *  We are forcing the usage of the buildAnisotropicMatrix method.
     *  The results are compared with Gauss' theorem and the results from the
     * buildIsotropicMatrix method
     *  The point charge is at the origin.
     */
    WHEN("the point group is D2h") {
      Molecule point = dummy<7>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius, "D2h");

      IEFSolver aniso_solver(symm);
      aniso_solver.buildAnisotropicMatrix(cavity, gf_i, gf_o, op);

      IEFSolver iso_solver(symm);
      iso_solver.buildIsotropicMatrix(cavity, gf_i, gf_o, op);

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

      THEN("the total apparent surface charge is") {
        Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);
        Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

        int nr_irrep = cavity.pointGroup().nrIrrep();
        double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
        double totalIsoASC = iso_fake_asc.sum() * nr_irrep;

        CAPTURE(totalASC);
        CAPTURE(totalAnisoASC);
        CAPTURE(totalASC - totalAnisoASC);
        REQUIRE(totalIsoASC == Approx(totalAnisoASC));
        REQUIRE(totalASC == Approx(totalAnisoASC).epsilon(1.0e-03));
      }
    }
  }
}
