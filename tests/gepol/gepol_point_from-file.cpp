#define BOOST_TEST_MODULE GePolCavityRestartTest

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "GePolCavity.hpp"

struct GePolCavityRestartTest {
    GePolCavity cavity;
    GePolCavityRestartTest() { SetUp(); }
    void SetUp() {
        cavity.loadCavity("point.npz");
    }
};

/*! \class GePolCavity
 *  \test \b GePolCavityRestartTest_size tests GePol cavity size for a point charge loading the cavity from a .npz file
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityRestartTest)
{
    int size = 32;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityRestartTest_area tests GePol cavity surface area for a point charge loading the cavity from from a .npz file
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityRestartTest)
{
    double area = 4.0 * M_PI * pow(1.0, 2);
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-12);
}

/*! \class GePolCavity
 *  \test \b GePolCavityRestartTest_volume tests GePol cavity volume for a point charge loading the cavity from from a .npz file
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityRestartTest)
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
