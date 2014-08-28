#define BOOST_TEST_MODULE GePolCavityNH3RestartTest

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <cmath>

#include "Config.hpp"

#include <Eigen/Dense>

#include "GePolCavity.hpp"

struct GePolCavityNH3RestartTest {
    GePolCavity cavity;
    GePolCavityNH3RestartTest() { SetUp(); }
    void SetUp() {
        cavity.loadCavity("nh3.npz");
    }
};

/*! \class GePolCavity
 *  \test \b GePolCavityNH3RestartTest_size tests GePol cavity size for ammonia loading the cavity from a .npz file
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityNH3RestartTest)
{
    int size = 544;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityNH3RestartTest_area tests GePol cavity surface area for ammonia loading the cavity from from a .npz file
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityNH3RestartTest)
{
    double area = 147.18581691164593;
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-10);
}

/*! \class GePolCavity
 *  \test \b GePolCavityNH3RestartTest_volume tests GePol cavity volume for ammonia loading the cavity from from a .npz file
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityNH3RestartTest)
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
