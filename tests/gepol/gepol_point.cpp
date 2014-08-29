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

/*! \class GePolCavity
 *  \test \b GePolCavityTest_size tests GePol cavity size for a point charge
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityTest)
{
    int size = 32;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityTest_area tests GePol cavity surface area for a point charge
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityTest)
{
    double area = 4.0 * M_PI * pow(1.0, 2);
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-12);
}

/*! \class GePolCavity
 *  \test \b GePolCavityTest_volume tests GePol cavity volume for a point charge
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityTest)
{
    double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
    Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
    Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cavity.size(); ++i ) {
        actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(
                            i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-12);
}
