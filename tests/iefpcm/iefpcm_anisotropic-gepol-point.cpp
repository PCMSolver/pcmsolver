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

#define BOOST_TEST_MODULE IEFSolveranisotropicPointChargeGePol

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <iostream>

#include "Config.hpp"

#include <Eigen/Core>

#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "GePolCavity.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "IEFSolver.hpp"
#include "TestingMolecules.hpp"

/*! \class IEFSolver
 *  \test \b anisotropicPointChargeGePol tests IEFSolver using a point charge with a GePol cavity
 *  We are forcing the usage of the buildAnisotropicMatrix method.
 *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
 *  The point charge is at the origin.
 */
BOOST_AUTO_TEST_CASE(anisotropicPointChargeGePol)
{
    // Set up cavity
    Molecule point = dummy<0>(2.929075493);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity = GePolCavity(point, area, probeRadius, minRadius);
    cavity.saveCavity("point.npz");

    double permittivity = 78.39;
    pcm::shared_ptr<Vacuum<AD_directional, CollocationIntegrator> > gfInside = 
	    pcm::make_shared<Vacuum<AD_directional, CollocationIntegrator> >(Vacuum<AD_directional, CollocationIntegrator>());
    pcm::shared_ptr<UniformDielectric<AD_directional, CollocationIntegrator> > gfOutside 
	    = pcm::make_shared<UniformDielectric<AD_directional, CollocationIntegrator> >(UniformDielectric<AD_directional, CollocationIntegrator>(permittivity));
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
    Eigen::VectorXd aniso_fake_asc = aniso_solver.computeCharge(fake_mep);

    Eigen::VectorXd iso_fake_asc = iso_solver.computeCharge(fake_mep);

    for (int i = 0; i < size; ++i) {
        BOOST_TEST_MESSAGE("fake_mep(" << i << ") = " << fake_mep(i));
    }
    for (int i = 0; i < size; ++i) {
        BOOST_TEST_MESSAGE("aniso_fake_asc(" << i << ") = " << aniso_fake_asc(i));
    }
    for (int i = 0; i < size; ++i) {
        BOOST_TEST_MESSAGE("iso_fake_asc(" << i << ") = " << iso_fake_asc(i));
    }

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
 *  \test \b anisotropicPointChargeShiftedGePol tests IEFSolver using a point charge with a GePol cavity
 *  We are forcing the usage of the buildAnisotropicMatrix method.
 *  The results are compared with Gauss' theorem and the results from the buildIsotropicMatrix method
 *  The point charge is away from the origin.
 */
BOOST_AUTO_TEST_CASE(anisotropicPointChargeShiftedGePol)
{
    // Set up cavity
    Eigen::Vector3d origin = 100 * Eigen::Vector3d::Random();
    Molecule point = dummy<0>(2.929075493, origin);
    double area = 0.4;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity = GePolCavity(point, area, probeRadius, minRadius);
    cavity.saveCavity("point.npz");

    double permittivity = 78.39;
    pcm::shared_ptr<Vacuum<AD_directional, CollocationIntegrator> > gfInside = 
	    pcm::make_shared<Vacuum<AD_directional, CollocationIntegrator> >(Vacuum<AD_directional, CollocationIntegrator>());
    pcm::shared_ptr<UniformDielectric<AD_directional, CollocationIntegrator> > gfOutside 
	    = pcm::make_shared<UniformDielectric<AD_directional, CollocationIntegrator> >(UniformDielectric<AD_directional, CollocationIntegrator>(permittivity));
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
        double distance = (center - origin).norm();
        fake_mep(i) = charge / distance;
    }
    Eigen::VectorXd aniso_fake_asc = Eigen::VectorXd::Zero(size);
    aniso_fake_asc = aniso_solver.computeCharge(fake_mep);

    Eigen::VectorXd iso_fake_asc = Eigen::VectorXd::Zero(size);
    iso_fake_asc = iso_solver.computeCharge(fake_mep);

    for (int i = 0; i < size; ++i) {
        BOOST_TEST_MESSAGE("fake_mep(" << i << ") = " << fake_mep(i));
    }
    for (int i = 0; i < size; ++i) {
        BOOST_TEST_MESSAGE("aniso_fake_asc(" << i << ") = " << aniso_fake_asc(i));
    }
    for (int i = 0; i < size; ++i) {
        BOOST_TEST_MESSAGE("iso_fake_asc(" << i << ") = " << iso_fake_asc(i));
    }

    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    double totalASC = - charge * (permittivity - 1) / permittivity;
    double totalAnisoASC = aniso_fake_asc.sum();
    double totalIsoASC = iso_fake_asc.sum();
    BOOST_TEST_MESSAGE("totalASC = " << totalASC);
    BOOST_TEST_MESSAGE("totalAnisoFakeASC = " << totalAnisoASC);
    BOOST_TEST_MESSAGE("totalASC - totalAnisoASC = " << totalASC - totalAnisoASC);
    BOOST_REQUIRE_CLOSE(totalASC, totalAnisoASC, 4e-02);
    BOOST_REQUIRE_CLOSE(totalIsoASC, totalAnisoASC, 1.0e-09);
}
