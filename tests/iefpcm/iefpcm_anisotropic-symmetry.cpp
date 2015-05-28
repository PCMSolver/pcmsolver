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

#define BOOST_TEST_MODULE IEFSolveranisotropicPointChargeGePolSymmetry

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <iostream>

#include "Config.hpp"

#include <Eigen/Dense>

#include <boost/filesystem.hpp>

#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "GePolCavity.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "IEFSolver.hpp"
#include "TestingMolecules.hpp"

namespace fs = boost::filesystem;

/*! \class IEFSolver
 *  \test \b anisotropicPointChargeGePolC1 tests IEFSolver using a point charge with a GePol cavity in C1 symmetry
 *  We are forcing the usage of the buildAnisotropicMatrix method.
 *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
 *  The point charge is at the origin.
 */
BOOST_AUTO_TEST_CASE(anisotropicPointChargeGePolC1)
{
    // Set up cavity
    Molecule point = dummy<0>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.c1");

    double permittivity = 78.39;
    Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> > * gfInside = new Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> >();
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> > * gfOutside = new
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> >(permittivity);
    bool symm = true;
    IEFSolver aniso_solver(gfInside, gfOutside, symm);
    aniso_solver.buildAnisotropicMatrix(cavity);

    IEFSolver iso_solver(gfInside, gfOutside, symm);
    iso_solver.buildIsotropicMatrix(cavity);

    double charge = 8.0;
    int size = cavity.size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    Eigen::VectorXd aniso_fake_asc = Eigen::VectorXd::Zero(size);
    aniso_solver.computeCharge(fake_mep, aniso_fake_asc);

    Eigen::VectorXd iso_fake_asc = Eigen::VectorXd::Zero(size);
    iso_solver.computeCharge(fake_mep, iso_fake_asc);

    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalAnisoASC = aniso_fake_asc.sum();
    double totalIsoASC = iso_fake_asc.sum();
    BOOST_TEST_MESSAGE("totalASC = " << totalASC);
    BOOST_TEST_MESSAGE("totalAnisoASC = " << totalAnisoASC);
    BOOST_TEST_MESSAGE("totalASC - totalAnisoASC = " << totalASC - totalAnisoASC);
    BOOST_REQUIRE_CLOSE(totalASC, totalAnisoASC, 4e-02);
    BOOST_REQUIRE_CLOSE(totalIsoASC, totalAnisoASC, 1.0e-09);
}

/*! \class IEFSolver
 *  \test \b anisotropicPointChargeGePolC2 tests IEFSolver using a point charge with a GePol cavity in C2 symmetry
 *  We are forcing the usage of the buildAnisotropicMatrix method.
 *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
 *  The point charge is at the origin.
 */
BOOST_AUTO_TEST_CASE(anisotropicPointChargeGePolC2)
{
    // Set up cavity
    Molecule point = dummy<1>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.c2");

    double permittivity = 78.39;
    Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> > * gfInside = new Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> >();
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> > * gfOutside = new
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> >(permittivity);
    bool symm = true;
    IEFSolver aniso_solver(gfInside, gfOutside, symm);
    aniso_solver.buildAnisotropicMatrix(cavity);

    IEFSolver iso_solver(gfInside, gfOutside, symm);
    iso_solver.buildIsotropicMatrix(cavity);

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

    Eigen::VectorXd aniso_fake_asc = Eigen::VectorXd::Zero(irr_size);
    aniso_solver.computeCharge(fake_mep, aniso_fake_asc);

    Eigen::VectorXd iso_fake_asc = Eigen::VectorXd::Zero(irr_size);
    iso_solver.computeCharge(fake_mep, iso_fake_asc);

    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
    double totalIsoASC = iso_fake_asc.sum() * nr_irrep;
    BOOST_TEST_MESSAGE("totalASC = " << totalASC);
    BOOST_TEST_MESSAGE("totalAnisoASC = " << totalAnisoASC);
    BOOST_TEST_MESSAGE("totalASC - totalAnisoASC = " << totalASC - totalAnisoASC);
    BOOST_REQUIRE_CLOSE(totalASC, totalAnisoASC, 4e-02);
    BOOST_REQUIRE_CLOSE(totalIsoASC, totalAnisoASC, 1.0e-09);
}

/*! \class IEFSolver
 *  \test \b anisotropicPointChargeGePolCs tests IEFSolver using a point charge with a GePol cavity in Cs symmetry
 *  We are forcing the usage of the buildAnisotropicMatrix method.
 *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
 *  The point charge is at the origin.
 */
BOOST_AUTO_TEST_CASE(anisotropicPointChargeGePolCs)
{
    // Set up cavity
    Molecule point = dummy<2>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.cs");

    double permittivity = 78.39;
    Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> > * gfInside = new Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> >();
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> > * gfOutside = new
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> >(permittivity);
    bool symm = true;
    IEFSolver aniso_solver(gfInside, gfOutside, symm);
    aniso_solver.buildAnisotropicMatrix(cavity);

    IEFSolver iso_solver(gfInside, gfOutside, symm);
    iso_solver.buildIsotropicMatrix(cavity);

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

    Eigen::VectorXd aniso_fake_asc = Eigen::VectorXd::Zero(irr_size);
    aniso_solver.computeCharge(fake_mep, aniso_fake_asc);

    Eigen::VectorXd iso_fake_asc = Eigen::VectorXd::Zero(irr_size);
    iso_solver.computeCharge(fake_mep, iso_fake_asc);

    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
    double totalIsoASC = iso_fake_asc.sum() * nr_irrep;
    BOOST_TEST_MESSAGE("totalASC = " << totalASC);
    BOOST_TEST_MESSAGE("totalAnisoASC = " << totalAnisoASC);
    BOOST_TEST_MESSAGE("totalASC - totalAnisoASC = " << totalASC - totalAnisoASC);
    BOOST_REQUIRE_CLOSE(totalASC, totalAnisoASC, 4e-02);
    BOOST_REQUIRE_CLOSE(totalIsoASC, totalAnisoASC, 1.0e-09);
}

/*! \class IEFSolver
 *  \test \b anisotropicPointChargeGePolCi tests IEFSolver using a point charge with a GePol cavity in Ci symmetry
 *  We are forcing the usage of the buildAnisotropicMatrix method.
 *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
 *  The point charge is at the origin.
 */
BOOST_AUTO_TEST_CASE(anisotropicPointChargeGePolCi)
{
    // Set up cavity
    Molecule point = dummy<3>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.ci");

    double permittivity = 78.39;
    Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> > * gfInside = new Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> >();
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> > * gfOutside = new
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> >(permittivity);
    bool symm = true;
    IEFSolver aniso_solver(gfInside, gfOutside, symm);
    aniso_solver.buildAnisotropicMatrix(cavity);

    IEFSolver iso_solver(gfInside, gfOutside, symm);
    iso_solver.buildIsotropicMatrix(cavity);

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

    Eigen::VectorXd aniso_fake_asc = Eigen::VectorXd::Zero(irr_size);
    aniso_solver.computeCharge(fake_mep, aniso_fake_asc);

    Eigen::VectorXd iso_fake_asc = Eigen::VectorXd::Zero(irr_size);
    iso_solver.computeCharge(fake_mep, iso_fake_asc);

    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
    double totalIsoASC = iso_fake_asc.sum() * nr_irrep;
    BOOST_TEST_MESSAGE("totalASC = " << totalASC);
    BOOST_TEST_MESSAGE("totalAnisoASC = " << totalAnisoASC);
    BOOST_TEST_MESSAGE("totalASC - totalAnisoASC = " << totalASC - totalAnisoASC);
    BOOST_REQUIRE_CLOSE(totalASC, totalAnisoASC, 4e-02);
    BOOST_REQUIRE_CLOSE(totalIsoASC, totalAnisoASC, 1.0e-09);
}

/*! \class IEFSolver
 *  \test \b anisotropicPointChargeGePolD2 tests IEFSolver using a point charge with a GePol cavity in D2 symmetry
 *  We are forcing the usage of the buildAnisotropicMatrix method.
 *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
 *  The point charge is at the origin.
 */
BOOST_AUTO_TEST_CASE(anisotropicPointChargeGePolD2)
{
    // Set up cavity
    Molecule point = dummy<4>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.d2");

    double permittivity = 78.39;
    Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> > * gfInside = new Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> >();
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> > * gfOutside = new
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> >(permittivity);
    bool symm = true;
    IEFSolver aniso_solver(gfInside, gfOutside, symm);
    aniso_solver.buildAnisotropicMatrix(cavity);

    IEFSolver iso_solver(gfInside, gfOutside, symm);
    iso_solver.buildIsotropicMatrix(cavity);

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

    Eigen::VectorXd aniso_fake_asc = Eigen::VectorXd::Zero(irr_size);
    aniso_solver.computeCharge(fake_mep, aniso_fake_asc);

    Eigen::VectorXd iso_fake_asc = Eigen::VectorXd::Zero(irr_size);
    iso_solver.computeCharge(fake_mep, iso_fake_asc);

    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
    double totalIsoASC = iso_fake_asc.sum() * nr_irrep;
    BOOST_TEST_MESSAGE("totalASC = " << totalASC);
    BOOST_TEST_MESSAGE("totalAnisoASC = " << totalAnisoASC);
    BOOST_TEST_MESSAGE("totalASC - totalAnisoASC = " << totalASC - totalAnisoASC);
    BOOST_REQUIRE_CLOSE(totalASC, totalAnisoASC, 4e-02);
    BOOST_REQUIRE_CLOSE(totalIsoASC, totalAnisoASC, 1.0e-09);
}

/*! \class IEFSolver
 *  \test \b anisotropicPointChargeGePolC2v tests IEFSolver using a point charge with a GePol cavity in C2v symmetry
 *  We are forcing the usage of the buildAnisotropicMatrix method.
 *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
 *  The point charge is at the origin.
 */
BOOST_AUTO_TEST_CASE(anisotropicPointChargeGePolC2v)
{
    // Set up cavity
    Molecule point = dummy<5>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.c2v");

    double permittivity = 78.39;
    Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> > * gfInside = new Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> >();
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> > * gfOutside = new
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> >(permittivity);
    bool symm = true;
    IEFSolver aniso_solver(gfInside, gfOutside, symm);
    aniso_solver.buildAnisotropicMatrix(cavity);

    IEFSolver iso_solver(gfInside, gfOutside, symm);
    iso_solver.buildIsotropicMatrix(cavity);

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

    Eigen::VectorXd aniso_fake_asc = Eigen::VectorXd::Zero(irr_size);
    aniso_solver.computeCharge(fake_mep, aniso_fake_asc);

    Eigen::VectorXd iso_fake_asc = Eigen::VectorXd::Zero(irr_size);
    iso_solver.computeCharge(fake_mep, iso_fake_asc);

    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
    double totalIsoASC = iso_fake_asc.sum() * nr_irrep;
    BOOST_TEST_MESSAGE("totalASC = " << totalASC);
    BOOST_TEST_MESSAGE("totalAnisoASC = " << totalAnisoASC);
    BOOST_TEST_MESSAGE("totalASC - totalAnisoASC = " << totalASC - totalAnisoASC);
    BOOST_REQUIRE_CLOSE(totalASC, totalAnisoASC, 4e-02);
    BOOST_REQUIRE_CLOSE(totalIsoASC, totalAnisoASC, 1.0e-09);
}

/*! \class IEFSolver
 *  \test \b anisotropicPointChargeGePolC2h tests IEFSolver using a point charge with a GePol cavity in C2h symmetry
 *  We are forcing the usage of the buildAnisotropicMatrix method.
 *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
 *  The point charge is at the origin.
 */
BOOST_AUTO_TEST_CASE(anisotropicPointChargeGePolC2h)
{
    // Set up cavity
    Molecule point = dummy<6>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.c2h");

    double permittivity = 78.39;
    Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> > * gfInside = new Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> >();
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> > * gfOutside = new
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> >(permittivity);
    bool symm = true;
    IEFSolver aniso_solver(gfInside, gfOutside, symm);
    aniso_solver.buildAnisotropicMatrix(cavity);

    IEFSolver iso_solver(gfInside, gfOutside, symm);
    iso_solver.buildIsotropicMatrix(cavity);

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

    Eigen::VectorXd aniso_fake_asc = Eigen::VectorXd::Zero(irr_size);
    aniso_solver.computeCharge(fake_mep, aniso_fake_asc);

    Eigen::VectorXd iso_fake_asc = Eigen::VectorXd::Zero(irr_size);
    iso_solver.computeCharge(fake_mep, iso_fake_asc);

    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
    double totalIsoASC = iso_fake_asc.sum() * nr_irrep;
    BOOST_TEST_MESSAGE("totalASC = " << totalASC);
    BOOST_TEST_MESSAGE("totalAnisoASC = " << totalAnisoASC);
    BOOST_TEST_MESSAGE("totalASC - totalAnisoASC = " << totalASC - totalAnisoASC);
    BOOST_REQUIRE_CLOSE(totalASC, totalAnisoASC, 4e-02);
    BOOST_REQUIRE_CLOSE(totalIsoASC, totalAnisoASC, 1.0e-09);
}

/*! \class IEFSolver
 *  \test \b anisotropicPointChargeGePolD2h tests IEFSolver using a point charge with a GePol cavity in D2h symmetry
 *  We are forcing the usage of the buildAnisotropicMatrix method.
 *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
 *  The point charge is at the origin.
 */
BOOST_AUTO_TEST_CASE(anisotropicPointChargeGePolD2h)
{
    // Set up cavity
    Molecule point = dummy<7>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity(point, area, probeRadius, minRadius);
    fs::rename("PEDRA.OUT", "PEDRA.OUT.d2h");

    double permittivity = 78.39;
    Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> > * gfInside = new Vacuum<AD_directional, CollocationIntegrator<AD_directional, Uniform> >();
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> > * gfOutside = new
    UniformDielectric<AD_directional, CollocationIntegrator<AD_directional, Uniform> >(permittivity);
    bool symm = true;
    IEFSolver aniso_solver(gfInside, gfOutside, symm);
    aniso_solver.buildAnisotropicMatrix(cavity);

    IEFSolver iso_solver(gfInside, gfOutside, symm);
    iso_solver.buildIsotropicMatrix(cavity);

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

    Eigen::VectorXd aniso_fake_asc = Eigen::VectorXd::Zero(irr_size);
    aniso_solver.computeCharge(fake_mep, aniso_fake_asc);

    Eigen::VectorXd iso_fake_asc = Eigen::VectorXd::Zero(irr_size);
    iso_solver.computeCharge(fake_mep, iso_fake_asc);

    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalAnisoASC = aniso_fake_asc.sum() * nr_irrep;
    double totalIsoASC = iso_fake_asc.sum() * nr_irrep;
    BOOST_TEST_MESSAGE("totalASC = " << totalASC);
    BOOST_TEST_MESSAGE("totalAnisoASC = " << totalAnisoASC);
    BOOST_TEST_MESSAGE("totalASC - totalAnisoASC = " << totalASC - totalAnisoASC);
    BOOST_REQUIRE_CLOSE(totalASC, totalAnisoASC, 4e-02);
    BOOST_REQUIRE_CLOSE(totalIsoASC, totalAnisoASC, 1.0e-09);
}
