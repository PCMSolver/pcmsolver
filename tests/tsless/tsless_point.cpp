#define BOOST_TEST_MODULE TsLessCavity

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <cmath>

#include "Config.hpp"

#include <Eigen/Dense>

#include "TsLessCavity.hpp"

struct TsLessCavityTest {
    TsLessCavity cav;
    TsLessCavityTest() { SetUp(); }
    void SetUp() {
        Eigen::Vector3d origin(0.0, 0.0, 0.0);
        std::vector<Sphere> spheres;
        Sphere sph1(origin,  1.0);
        spheres.push_back(sph1);
        double area = 0.4;
        double minDistance = 0.1;
        double probeRadius = 0.0;
        int derOrder = 4;
        cav = TsLessCavity(spheres, area, probeRadius, minDistance, derOrder);
    }
};

BOOST_FIXTURE_TEST_CASE(size, TsLessCavityTest)
{
    int size = 32;
    int actualSize = cav.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

BOOST_FIXTURE_TEST_CASE(area, TsLessCavityTest)
{
    double area = 4.0 * M_PI * pow(1.0, 2);
    double actualArea = cav.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-12);
}

BOOST_FIXTURE_TEST_CASE(volume, TsLessCavityTest)
{
    double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
    Eigen::Matrix3Xd elementCenter = cav.elementCenter();
    Eigen::Matrix3Xd elementNormal = cav.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cav.size(); ++i ) {
        actualVolume += cav.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-12);
}
