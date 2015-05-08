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

#define BOOST_TEST_MODULE IEFSolverpointChargeGePol

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <iostream>

#include "Config.hpp"

#include <Eigen/Dense>

#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "GePolCavity.hpp"
#include "Vacuum.hpp"
#include "TanhSphericalDiffuse.hpp"
#include "IEFSolver.hpp"
#include "TestingMolecules.hpp"

/*! \class IEFSolver
 *  \test \b pointChargeDiffuseGePol tests IEFSolver using a point charge with a GePol cavity and a spherical diffuse interface
 */
BOOST_AUTO_TEST_CASE(pointChargeGePol)
{
    // Set up cavity
    Molecule point = dummy<0>(2.929075493);
    double area = 100.0;
    double probeRadius = 0.0;
    double minRadius = 100.0;
    GePolCavity cavity = GePolCavity(point, area, probeRadius, minRadius);
    cavity.saveCavity("point.npz");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double eps1 = 78.39;
    double eps2 = 78.39;
    double center = 100.0;
    double width = 5.0;
    Eigen::Vector3d origin = 100 * Eigen::Vector3d::Random();
    std::cout << "origin " << origin.transpose() << std::endl;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(diag);
    TanhSphericalDiffuse * gfOutside = new TanhSphericalDiffuse(eps1, eps2, width, center, origin, diag);
    bool symm = true;
    IEFSolver solver(gfInside, gfOutside, symm);
    solver.buildSystemMatrix(cavity);

    double charge = 8.0;
    int size = cavity.size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double distance = center.norm();
        fake_mep(i) = charge / distance;
    }
    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
    solver.computeCharge(fake_mep, fake_asc);
    double totalASC = - charge * (eps1 - 1) / eps1;
    double totalFakeASC = fake_asc.sum();
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 4e-02);
}
