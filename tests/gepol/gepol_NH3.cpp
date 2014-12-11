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

#define BOOST_TEST_MODULE GePolCavityNH3

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <cmath>

#include "Config.hpp"

#include <Eigen/Dense>

#include "GePolCavity.hpp"
#include "Molecule.hpp"
#include "PhysicalConstants.hpp"
#include "Symmetry.hpp"

struct GePolCavityNH3Test {
    GePolCavity cavity;
    GePolCavityNH3Test() { SetUp(); }
    void SetUp() {
        Eigen::Vector3d N( -0.000000000,   -0.104038047,    0.000000000);
        Eigen::Vector3d H1(-0.901584415,    0.481847022,   -1.561590016);
        Eigen::Vector3d H2(-0.901584415,    0.481847022,    1.561590016);
        Eigen::Vector3d H3( 1.803168833,    0.481847022,    0.000000000);
	
	int nAtoms = 4;
	Eigen::MatrixXd geom(3, nAtoms);
	geom.col(0) = N.transpose();
	geom.col(1) = H1.transpose();
	geom.col(2) = H2.transpose();
	geom.col(3) = H3.transpose();
	Eigen::Vector4d charges, masses;
	charges << 7.0, 1.0, 1.0, 1.0;
	masses  << 14.0030740, 1.0078250, 1.0078250, 1.0078250;
	std::vector<Atom> atoms;
	atoms.push_back( Atom("Nitrogen", "N", charges(0), masses(0), 2.929075493, N,  1.0) );
	atoms.push_back( Atom("Hydrogen", "H", charges(1), masses(1), 2.267671349, H1, 1.0) );
	atoms.push_back( Atom("Hydrogen", "H", charges(2), masses(2), 2.267671349, H2, 1.0) );
	atoms.push_back( Atom("Hydrogen", "H", charges(3), masses(3), 2.267671349, H3, 1.0) );

        std::vector<Sphere> spheres;
        Sphere sph1(N,  2.929075493);
        Sphere sph2(H1, 2.267671349);
        Sphere sph3(H2, 2.267671349);
        Sphere sph4(H3, 2.267671349);
        spheres.push_back(sph1);
        spheres.push_back(sph2);
        spheres.push_back(sph3);
        spheres.push_back(sph4);

	Molecule molec(nAtoms, charges, masses, geom, atoms, spheres);
        double area = 0.4; // Bohr^2
        double probeRadius = 1.385 / convertBohrToAngstrom;
        double minRadius = 0.2 / convertBohrToAngstrom;
        // C1
        Symmetry pGroup = buildGroup(0, 0, 0, 0);
        cavity = GePolCavity(molec, area, probeRadius, minRadius, pGroup);
        cavity.saveCavity("nh3.npz");
    }
};

/*! \class GePolCavity
 *  \test \b GePolCavityNH3Test_size tests GePol cavity size for ammonia
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityNH3Test)
{
    int size = 544;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityNH3Test_area tests GePol cavity surface area for ammonia
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityNH3Test)
{
    double area = 147.18581691164593;
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-10);
}

/*! \class GePolCavity
 *  \test \b GePolCavityNH3Test_volume tests GePol cavity volume for ammonia
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityNH3Test)
{
    double volume = 152.81441857040116;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-10);
}
