#define BOOST_TEST_MODULE TsLessCavityNH3

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "TsLessCavity.hpp"

struct TsLessCavityNH3Test {
protected:
    TsLessCavity cav;
    TsLessCavityNH3Test() { SetUp(); }
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
        double area = 0.4;
        double minDistance = 0.1;
        double probeRadius = 0.0;
        int derOrder = 4;
        cav = TsLessCavity(spheres, area, probeRadius, minDistance, derOrder);
    }
};

BOOST_FIXTURE_TEST_CASE(size, TsLessCavityNH3Test)
{
    int size = 544;
    int actualSize = cav.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

BOOST_FIXTURE_TEST_CASE(area, TsLessCavityNH3Test)
{
    double area = 147.18581691164593;
    double actualArea = cav.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-12);
}

BOOST_FIXTURE_TEST_CASE(volume, TsLessCavityNH3Test)
{
    double volume = 152.81441857040116;
    Eigen::Matrix3Xd elementCenter = cav.elementCenter();
    Eigen::Matrix3Xd elementNormal = cav.elementNormal();
    double actualVolume = 0;
    for ( int i = 0; i < cav.size(); ++i ) {
        actualVolume += cav.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
    }
    actualVolume /= 3;
    BOOST_REQUIRE_CLOSE(volume, actualVolume, 1.0e-12);
}
