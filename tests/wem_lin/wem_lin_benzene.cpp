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

#define BOOST_TEST_MODULE WEMSolverbenzene
#define DEBUG

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
 *  \test \b  benzene tests WEMSolver with linear ansatz functions using ammonia and a wavelet cavity
 */
BOOST_AUTO_TEST_CASE(benzene)
{
  printf("STARTING TEST\n");
    // Set up cavity
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
    */
    Eigen::Vector3d C1( 5.274,  1.999, -8.568);
    Eigen::Vector3d C2( 6.627,  2.018, -8.209);
    Eigen::Vector3d C3( 7.366,  0.829, -8.202);
    Eigen::Vector3d C4( 6.752, -0.379, -8.554);
    Eigen::Vector3d C5( 5.399, -0.398, -8.912);
    Eigen::Vector3d C6( 4.660,  0.791, -8.919);
    Eigen::Vector3d H1( 4.704,  2.916, -8.573);
    Eigen::Vector3d H2( 7.101,  2.950, -7.938);
    Eigen::Vector3d H3( 8.410,  0.844, -7.926);
    Eigen::Vector3d H4( 7.322, -1.296, -8.548);
    Eigen::Vector3d H5( 4.925, -1.330, -9.183);
    Eigen::Vector3d H6( 3.616,  0.776, -9.196);

    std::vector<Sphere> spheres;
    /*
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
    */
    Sphere sph1(C1, 1.7*1.2); 
    Sphere sph2(C2, 1.7*1.2);
    Sphere sph3(C3, 1.7*1.2);
    Sphere sph4(C4, 1.7*1.2);
    Sphere sph5(C5, 1.7*1.2);
    Sphere sph6(C6, 1.7*1.2);
    
    Sphere sph7(H1, 1.2*1.2); 
    Sphere sph8(H2, 1.2*1.2);
    Sphere sph9(H3, 1.2*1.2);
    Sphere sph10(H4, 1.2*1.2);
    Sphere sph11(H5, 1.2*1.2);
    Sphere sph12(H6, 1.2*1.2);

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
    
    double probeRadius = 1.385*1.2; // Probe Radius for water
    int patchLevel = 3;
    double coarsity = 0.5;
  printf("TEST 1\n");
    WaveletCavity cavity(spheres, probeRadius, patchLevel, coarsity);
  printf("TEST 1a\n");
    cavity.readCavity("molec_dyadic.dat");
    cavity.scaleNPoints(1./0.52917721092);

  printf("TEST 2\n");
    CollocationIntegrator * diag = new CollocationIntegrator();
    double permittivity = 78.39;
    Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>();
    UniformDielectric<AD_directional> * gfOutside = new
    UniformDielectric<AD_directional>(permittivity);
    int firstKind = 0;
#ifdef DEBUG
    FILE* debugFile = fopen("debug.out","w");
    fclose(debugFile);
#endif
    WEMSolver solver(gfInside, gfOutside, "Linear", firstKind);
    solver.buildSystemMatrix(cavity);
    cavity.uploadPoints(solver.getQuadratureLevel(), solver.getT_());

    double Hcharge = 1.0;
    double Ccharge = 6.0;
    int size = cavity.size();
    Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; ++i) {
        Eigen::Vector3d center = cavity.elementCenter(i);
        double C1mep = Ccharge/(center - C1).norm();
        double C2mep = Ccharge/(center - C2).norm();
        double C3mep = Ccharge/(center - C3).norm();
        double C4mep = Ccharge/(center - C4).norm();
        double C5mep = Ccharge/(center - C5).norm();
        double C6mep = Ccharge/(center - C6).norm();

        double H1mep = Hcharge/(center - H1).norm();
        double H2mep = Hcharge/(center - H2).norm();
        double H3mep = Hcharge/(center - H3).norm();
        double H4mep = Hcharge/(center - H4).norm();
        double H5mep = Hcharge/(center - H5).norm();
        double H6mep = Hcharge/(center - H6).norm();
        fake_mep(i) = C1mep + C2mep + C3mep + C4mep + C5mep + C6mep +
          H1mep + H2mep + H3mep + H4mep + H5mep + H6mep;
    }
    // The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
    Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
    solver.compCharge(fake_mep, fake_asc);
    double totalASC = - (6 * Ccharge + 6 * Hcharge) * ( permittivity - 1) / permittivity; 
    double totalFakeASC = fake_asc.sum();
    std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
    BOOST_REQUIRE_CLOSE(totalASC, totalFakeASC, 3e-2);
}
