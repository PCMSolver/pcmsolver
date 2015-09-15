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

#define BOOST_TEST_MODULE GePolCavity

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <cmath>

#include "Config.hpp"

#include <Eigen/Dense> 

#include "GePolCavity.hpp"
#include "Symmetry.hpp"

struct GePolCavityTest {
    GePolCavity cavity;
    GePolCavityTest() { SetUp(); }
    void SetUp() {
        Eigen::Vector3d origin(0.0, 0.0, 0.0);
        std::vector<Sphere> spheres;
        Sphere sph1(origin,  1.0);
        spheres.push_back(sph1);
        double area = 0.4;
        // C1
        Symmetry pGroup = buildGroup(0, 0, 0, 0);
        cavity = GePolCavity(spheres, area, 0.0, 100.0, pGroup);
        cavity.saveCavity("point.npz");
    }
};

/* \class GePolCavity
 *  \test \b GePolCavityTest_size tests GePol cavity size for a point charge
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityTest)
{
    int size = 32;
    size_t actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/* \class GePolCavity
 *  \test \b GePolCavityTest_area tests GePol cavity surface area for a point charge
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityTest)
{
    double area = 4.0 * M_PI * pow(1.0, 2);
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-12);
}

/* \class GePolCavity
 *  \test \b GePolCavityTest_volume tests GePol cavity volume for a point charge
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityTest)
{
    double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( size_t i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-12);
}
