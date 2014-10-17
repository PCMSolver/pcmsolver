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

#define BOOST_TEST_MODULE WEMSolverNH3

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <iostream>
#include <vector>

#include "Config.hpp"

#include "CollocationIntegrator.hpp"
#include "DerivativeTypes.hpp"
#include "WEMSolver.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"
#include "WaveletCavity.hpp"

/*! \class WEMSolver
 *  \test \b NH3 tests WEMSolver with linear ansatz functions using ammonia and a wavelet cavity
 */
BOOST_AUTO_TEST_CASE(NH3)
{
    // Set up cavity

    Eigen::Vector3d N( -0.000000000,   -0.104038047,    0.000000000);
    Eigen::Vector3d H1(-0.901584415,    0.481847022,   -1.561590016);
    Eigen::Vector3d H2(-0.901584415,    0.481847022,    1.561590016);
    Eigen::Vector3d H3( 1.803168833,    0.481847022,    0.000000000);
    std::vector<Sphere> spheres;
    Sphere sph1(N,  2.929075493);
    Sphere sph2(H1, 2.267671349);
    Sphere sph3(H2, 2.267671349);
    Sphere sph4(H3, 2.267671349);
    spheres.push_back(sph1);
    spheres.push_back(sph2);
    spheres.push_back(sph3);
    spheres.push_back(sph4);
/*
    Eigen::Vector3d C1(-0.694303272975, -0.000000000000, -1.202568545351);
    Eigen::Vector3d C2( 0.694303272975,  0.000000000000, -1.202568545351);
    Eigen::Vector3d C3( 1.388606546154,  0.000000000000,  0.000000000000);
    Eigen::Vector3d C4( 0.694303272975,  0.000000000000,  1.202568545351);
    Eigen::Vector3d C5(-0.694303272975, -0.000000000000,  1.202568545351);
    Eigen::Vector3d C6(-1.388606546154, -0.000000000000,  0.000000000000);
    Eigen::Vector3d H1(-1.235418032354, -0.000000000000, -2.139806800843);
    Eigen::Vector3d H2( 1.235418032354,  0.000000000000, -2.139806800843);
    Eigen::Vector3d H3( 2.470836065313,  0.000000000000,  0.000000000000);
    Eigen::Vector3d H4( 1.235418032354,  0.000000000000,  2.139806800843);
    Eigen::Vector3d H5(-1.235418032354, -0.000000000000,  2.139806800843);
    Eigen::Vector3d H6(-2.470836065313, -0.000000000000,  0.000000000000);
    
    std::vector<Sphere> spheres;
    Sphere sph1(C1, 3.212534412); 
    Sphere sph2(C2, 3.212534412);
    Sphere sph3(C3, 3.212534412);
    Sphere sph4(C4, 3.212534412);
    Sphere sph5(C5, 3.212534412);
    Sphere sph6(C6, 3.212534412);
    
    Sphere sph7(H1, 2.267671349); 
    Sphere sph8(H2, 2.267671349);
    Sphere sph9(H3, 2.267671349);
    Sphere sph10(H4, 2.267671349);
    Sphere sph11(H5, 2.267671349);
    Sphere sph12(H6, 2.267671349);
    
    spheres.push_back(sph1);
    spheres.push_back(sph2);
    spheres.push_back(sph3);
    spheres.push_back(sph4);
    spheres.push_back(sph5);
    spheres.push_back(sph6);
    spheres.push_back(sph7);
    spheres.push_back(sph8);
    spheres.push_back(sph9);
    spheres.push_back(sph10);
    spheres.push_back(sph11);
    spheres.push_back(sph12);
 */
    double probeRadius = 1.385; // Probe Radius for water
    int patchLevel = 2;
    double coarsity = 0.5;
    WaveletCavity cavity(spheres, probeRadius, patchLevel, coarsity);
    cavity.readCavity("molec_dyadic.dat");

    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>();
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity);
    int firstKind = 0;
    FILE* debugFile = fopen("debug.out","w");
    fclose(debugFile);
    WEMSolver solver(gfInside, gfOutside, "Linear", firstKind);
    solver.buildSystemMatrix(cavity);
    cavity.uploadPoints(solver.getQuadratureLevel(), solver.getT_());

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
    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
    solver.compCharge(fake_mep, fake_asc);
    double totalASC = - (Ncharge + 3.0 * Hcharge) * (permittivity - 1) / permittivity;
    double totalFakeASC = fake_asc.sum();
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 3e-2);
}
