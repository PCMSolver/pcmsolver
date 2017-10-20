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

using namespace pcm;
using bi_operators::Collocation;
using cavity::GePolCavity;
using green::Vacuum;
using green::UniformDielectric;
using solver::CPCMSolver;

SCENARIO(
    "Test solver for the C-PCM for a point charge in different Abelian point groups",
    "[solver][cpcm][cpcm_symmetry]") {
  GIVEN("An isotropic environment modelled as a perfect conductor and a point "
        "charge") {
    double permittivity = 78.39;
    Vacuum<> gf_i;
    UniformDielectric<> gf_o(permittivity);
    bool symm = true;
    double correction = 0.0;

    Collocation S;

    double charge = 8.0;
    // The total ASC for a conductor is -Q
    // for CPCM it will be -Q*(epsilon-1)/epsilon
    double totalASC = -charge * (permittivity - 1) / permittivity;

    /*! \class CPCMSolver
     *  \test \b pointChargeGePolC1 tests CPCMSolver using a point charge with a
     * GePol cavity in C1 symmetry
     */
    WHEN("the point group is C1") {
      Molecule point = dummy<0>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius, "C1");

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

      CPCMSolver solver(symm, correction);
      solver.buildSystemMatrix(cavity, gf_i, gf_o, S);
      THEN("the total apparent surface charge is") {
        int irr_size = cavity.irreducible_size();
        Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
        fake_asc = solver.computeCharge(fake_mep);
        int nr_irrep = cavity.pointGroup().nrIrrep();
        double totalFakeASC = fake_asc.sum() * nr_irrep;
        CAPTURE(totalASC - totalFakeASC);
        REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
      }
    }

    /*! \class CPCMSolver
     *  \test \b pointChargeGePolC2 tests CPCMSolver using a point charge with a
     * GePol cavity in C2 symmetry
     */
    WHEN("the point group is C2") {
      Molecule point = dummy<1>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius, "C2");

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);

      CPCMSolver solver(symm, correction);
      solver.buildSystemMatrix(cavity, gf_i, gf_o, S);
      THEN("the total apparent surface charge is") {
        int irr_size = cavity.irreducible_size();
        Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
        fake_asc = solver.computeCharge(fake_mep);
        int nr_irrep = cavity.pointGroup().nrIrrep();
        double totalFakeASC = fake_asc.sum() * nr_irrep;
        CAPTURE(totalASC - totalFakeASC);
        REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
      }
    }

    /*! \class CPCMSolver
     *  \test \b pointChargeGePolCs tests CPCMSolver using a point charge with a
     * GePol cavity in Cs symmetry
     */
    WHEN("the point group is Cs") {
      Molecule point = dummy<2>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius, "Cs");

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);
      for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
      }

      CPCMSolver solver(symm, correction);
      solver.buildSystemMatrix(cavity, gf_i, gf_o, S);
      THEN("the total apparent surface charge is") {
        int irr_size = cavity.irreducible_size();
        Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
        fake_asc = solver.computeCharge(fake_mep);
        int nr_irrep = cavity.pointGroup().nrIrrep();
        double totalFakeASC = fake_asc.sum() * nr_irrep;
        CAPTURE(totalASC - totalFakeASC);
        REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
      }
    }

    /*! \class CPCMSolver
     *  \test \b pointChargeGePolCi tests CPCMSolver using a point charge with a
     * GePol cavity in Ci symmetry
     */
    WHEN("the point group is Ci") {
      Molecule point = dummy<3>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius, "Ci");

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);
      for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
      }

      CPCMSolver solver(symm, correction);
      solver.buildSystemMatrix(cavity, gf_i, gf_o, S);
      THEN("the total apparent surface charge is") {
        int irr_size = cavity.irreducible_size();
        Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
        fake_asc = solver.computeCharge(fake_mep);
        int nr_irrep = cavity.pointGroup().nrIrrep();
        double totalFakeASC = fake_asc.sum() * nr_irrep;
        CAPTURE(totalASC - totalFakeASC);
        REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
      }
    }

    /*! \class CPCMSolver
     *  \test \b pointChargeGePolD2 tests CPCMSolver using a point charge with a
     * GePol cavity in D2 symmetry
     */
    WHEN("the point group is D2") {
      Molecule point = dummy<4>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius, "D2");

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);
      for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
      }

      CPCMSolver solver(symm, correction);
      solver.buildSystemMatrix(cavity, gf_i, gf_o, S);
      THEN("the total apparent surface charge is") {
        int irr_size = cavity.irreducible_size();
        Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
        fake_asc = solver.computeCharge(fake_mep);
        int nr_irrep = cavity.pointGroup().nrIrrep();
        double totalFakeASC = fake_asc.sum() * nr_irrep;
        CAPTURE(totalASC - totalFakeASC);
        REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
      }
    }

    /*! \class CPCMSolver
     *  \test \b pointChargeGePolC2v tests CPCMSolver using a point charge with a
     * GePol cavity in C2v symmetry
     */
    WHEN("the point group is C2v") {
      Molecule point = dummy<5>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius, "C2v");

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);
      for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
      }

      CPCMSolver solver(symm, correction);
      solver.buildSystemMatrix(cavity, gf_i, gf_o, S);
      THEN("the total apparent surface charge is") {
        int irr_size = cavity.irreducible_size();
        Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
        fake_asc = solver.computeCharge(fake_mep);
        int nr_irrep = cavity.pointGroup().nrIrrep();
        double totalFakeASC = fake_asc.sum() * nr_irrep;
        CAPTURE(totalASC - totalFakeASC);
        REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
      }
    }

    /*! \class CPCMSolver
     *  \test \b pointChargeGePolC2h tests CPCMSolver using a point charge with a
     * GePol cavity in C2h symmetry
     */
    WHEN("the point group is C2h") {
      Molecule point = dummy<6>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius, "C2h");

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);
      for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
      }

      CPCMSolver solver(symm, correction);
      solver.buildSystemMatrix(cavity, gf_i, gf_o, S);
      THEN("the total apparent surface charge is") {
        int irr_size = cavity.irreducible_size();
        Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
        fake_asc = solver.computeCharge(fake_mep);
        int nr_irrep = cavity.pointGroup().nrIrrep();
        double totalFakeASC = fake_asc.sum() * nr_irrep;
        CAPTURE(totalASC - totalFakeASC);
        REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
      }
    }

    /*! \class CPCMSolver
     *  \test \b pointChargeGePolD2h tests CPCMSolver using a point charge with a
     * GePol cavity in D2h symmetry
     */
    WHEN("the point group is D2h") {
      Molecule point = dummy<7>(2.929075493);
      double area = 0.4;
      double probeRadius = 0.0;
      double minRadius = 100.0;
      GePolCavity cavity(point, area, probeRadius, minRadius, "D2h");

      int size = cavity.size();
      Eigen::VectorXd fake_mep = computeMEP(cavity.elements(), charge);
      for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
      }

      CPCMSolver solver(symm, correction);
      solver.buildSystemMatrix(cavity, gf_i, gf_o, S);
      THEN("the total apparent surface charge is") {
        int irr_size = cavity.irreducible_size();
        Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
        fake_asc = solver.computeCharge(fake_mep);
        int nr_irrep = cavity.pointGroup().nrIrrep();
        double totalFakeASC = fake_asc.sum() * nr_irrep;
        CAPTURE(totalASC - totalFakeASC);
        REQUIRE(totalASC == Approx(totalFakeASC).epsilon(1.0e-03));
      }
    }
  }
}
