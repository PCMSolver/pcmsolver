/* pcmsolver_copyright_start */
/*
 *     PCMSolver, an API for the Polarizable Continuum Model
 *     Copyright (C) 2013 Roberto Di Remigio, Luca Frediani and contributors
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

#define BOOST_TEST_MODULE CPCMSolverpointChargeGePolSymmetry

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <iostream>

#include "Config.hpp"

#include <Eigen/Dense>
#include <boost/filesystem.hpp>

#include "CPCMSolver.hpp"
#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "GePolCavity.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"
#include "TestingMolecules.hpp"

namespace fs = boost::filesystem;

/*! \class CPCMSolver
 *  \test \b pointChargeGePolC1 tests CPCMSolver using a point charge with a GePol cavity in C1 symmetry
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolC1)
{
    // Set up cavity
    Molecule point = dummy<0>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.c1");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    double correction = 0.0;
    CPCMSolver solver(gfInside, gfOutside, symm, correction);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int size = cavity.size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    // The total ASC for a conductor is -Q
    // for CPCM it will be -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
    solver.computeCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum();
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}

/*! \class CPCMSolver
 *  \test \b pointChargeGePolC2 tests CPCMSolver using a point charge with a GePol cavity in C2 symmetry
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolC2)
{
    // Set up cavity
    Molecule point = dummy<1>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.c2");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    double correction = 0.0;
    CPCMSolver solver(gfInside, gfOutside, symm, correction);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int irr_size = cavity.irreducible_size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(irr_size);
    // Calculate it only on the irreducible portion of the cavity
    // then replicate it according to the point group
    for (int i = 0; i < irr_size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    int nr_irrep = cavity.pointGroup().nrIrrep();
    // The total ASC for a conductor is -Q
    // for CPCM it will be -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(irr_size);
    solver.computeCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum() * nr_irrep;
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}

/*! \class CPCMSolver
 *  \test \b pointChargeGePolCs tests CPCMSolver using a point charge with a GePol cavity in Cs symmetry
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolCs)
{
    // Set up cavity
    Molecule point = dummy<2>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.cs");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    double correction = 0.0;
    CPCMSolver solver(gfInside, gfOutside, symm, correction);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int irr_size = cavity.irreducible_size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(irr_size);
    // Calculate it only on the irreducible portion of the cavity
    // then replicate it according to the point group
    for (int i = 0; i < irr_size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    int nr_irrep = cavity.pointGroup().nrIrrep();
    // The total ASC for a conductor is -Q
    // for CPCM it will be -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(irr_size);
    solver.computeCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum() * nr_irrep;
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}

/*! \class CPCMSolver
 *  \test \b pointChargeGePolCi tests CPCMSolver using a point charge with a GePol cavity in Ci symmetry
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolCi)
{
    // Set up cavity
    Molecule point = dummy<3>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.ci");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    double correction = 0.0;
    CPCMSolver solver(gfInside, gfOutside, symm, correction);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int irr_size = cavity.irreducible_size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(irr_size);
    // Calculate it only on the irreducible portion of the cavity
    // then replicate it according to the point group
    for (int i = 0; i < irr_size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    int nr_irrep = cavity.pointGroup().nrIrrep();
    // The total ASC for a conductor is -Q
    // for CPCM it will be -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(irr_size);
    solver.computeCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum() * nr_irrep;
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}

/*! \class CPCMSolver
 *  \test \b pointChargeGePolD2 tests CPCMSolver using a point charge with a GePol cavity in D2 symmetry
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolD2)
{
    // Set up cavity
    Molecule point = dummy<4>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.d2");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    double correction = 0.0;
    CPCMSolver solver(gfInside, gfOutside, symm, correction);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int irr_size = cavity.irreducible_size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(irr_size);
    // Calculate it only on the irreducible portion of the cavity
    // then replicate it according to the point group
    for (int i = 0; i < irr_size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    int nr_irrep = cavity.pointGroup().nrIrrep();
    // The total ASC for a conductor is -Q
    // for CPCM it will be -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(irr_size);
    solver.computeCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum() * nr_irrep;
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}

/*! \class CPCMSolver
 *  \test \b pointChargeGePolC2v tests CPCMSolver using a point charge with a GePol cavity in C2v symmetry
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolC2v)
{
    // Set up cavity
    Molecule point = dummy<5>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.c2v");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    double correction = 0.0;
    CPCMSolver solver(gfInside, gfOutside, symm, correction);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int irr_size = cavity.irreducible_size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(irr_size);
    // Calculate it only on the irreducible portion of the cavity
    // then replicate it according to the point group
    for (int i = 0; i < irr_size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    int nr_irrep = cavity.pointGroup().nrIrrep();
    // The total ASC for a conductor is -Q
    // for CPCM it will be -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(irr_size);
    solver.computeCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum() * nr_irrep;
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}

/*! \class CPCMSolver
 *  \test \b pointChargeGePolC2h tests CPCMSolver using a point charge with a GePol cavity in C2h symmetry
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolC2h)
{
    // Set up cavity
    Molecule point = dummy<6>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
        GePolCavity cavity(point, area, probeRadius, minRadius);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.c2h");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    double correction = 0.0;
    CPCMSolver solver(gfInside, gfOutside, symm, correction);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int irr_size = cavity.irreducible_size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(irr_size);
    // Calculate it only on the irreducible portion of the cavity
    // then replicate it according to the point group
    for (int i = 0; i < irr_size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    int nr_irrep = cavity.pointGroup().nrIrrep();
    // The total ASC for a conductor is -Q
    // for CPCM it will be -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(irr_size);
    solver.computeCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum() * nr_irrep;
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}

/*! \class CPCMSolver
 *  \test \b pointChargeGePolD2h tests CPCMSolver using a point charge with a GePol cavity in D2h symmetry
 */
BOOST_AUTO_TEST_CASE(pointChargeGePolD2h)
{
    // Set up cavity
    Molecule point = dummy<7>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.d2h");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity, diag);
    bool symm = true;
    double correction = 0.0;
    CPCMSolver solver(gfInside, gfOutside, symm, correction);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int irr_size = cavity.irreducible_size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(irr_size);
    // Calculate it only on the irreducible portion of the cavity
    // then replicate it according to the point group
    for (int i = 0; i < irr_size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    int nr_irrep = cavity.pointGroup().nrIrrep();
    // The total ASC for a conductor is -Q
    // for CPCM it will be -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(irr_size);
    solver.computeCharge(fake_mep, fake_asc);
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum() * nr_irrep;
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}
