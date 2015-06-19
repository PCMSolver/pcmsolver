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

#define BOOST_TEST_MODULE GePolCavityC6H6AddTest

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>

#include "Config.hpp"

#include <Eigen/Core>

#include <boost/filesystem.hpp>

#include "GePolCavity.hpp"
#include "Molecule.hpp"
#include "PhysicalConstants.hpp"
#include "Symmetry.hpp"
#include "TestingMolecules.hpp"

namespace fs = boost::filesystem;

struct GePolCavityC6H6AddTest {
    GePolCavity cavity;
    GePolCavityC6H6AddTest() { SetUp(); }
    void SetUp() {
        double area = 0.3 / convertBohr2ToAngstrom2;
        double probeRadius = 1.385 / convertBohrToAngstrom;
        // Addition of spheres is enabled, but will not happen in this particular case
        double minRadius = 10.0 / convertBohrToAngstrom;
	Molecule molec = C6H6();
        cavity = GePolCavity(molec, area, probeRadius, minRadius);
        cavity.saveCavity("c6h6.npz");
        fs::rename("PEDRA.OUT", "PEDRA.OUT.c6h6");
        fs::rename("cavity.off", "cavity.off.c6h6");
    }
};

/*! \class GePolCavity
 *  \test \b GePolCavityC6H6AddTest_size tests GePol cavity size for C6H6 in C1 symmetry with added spheres
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityC6H6AddTest)
{
    int size = 644;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC6H6AddTest_area tests GePol cavity surface area for C6H6 in C1 symmetry with added spheres
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityC6H6AddTest)
{
    double area = 391.06094362589062;
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-10);
}

/*! \class GePolCavity
 *  \test \b GePolCavityC6H6AddTest_volume tests GePol cavity volume for C6H6 in C1 symmetry with added spheres
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityC6H6AddTest)
{
    double volume = 567.67444339622284;
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
