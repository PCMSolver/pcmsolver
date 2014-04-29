#define BOOST_TEST_MODULE GePolCavityNH3

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "GePolCavity.hpp"
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
        std::vector<Sphere> spheres;
        Sphere sph1(N,  2.929075493);
        Sphere sph2(H1, 2.267671349);
        Sphere sph3(H2, 2.267671349);
        Sphere sph4(H3, 2.267671349);
        spheres.push_back(sph1);
        spheres.push_back(sph2);
        spheres.push_back(sph3);
        spheres.push_back(sph4);
        double area = 0.4; // Bohr^2
        double probeRadius = 1.385 / convertBohrToAngstrom;
        double minRadius = 0.2 / convertBohrToAngstrom;
        // C1
        Symmetry pGroup = buildGroup(0, 0, 0, 0);
        cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
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
