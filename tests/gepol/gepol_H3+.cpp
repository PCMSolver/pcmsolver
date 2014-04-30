#define BOOST_TEST_MODULE GePolCavityH3+

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "GePolCavity.hpp"
#include "PhysicalConstants.hpp"
#include "Symmetry.hpp"

struct GePolCavityH3Test {
    GePolCavity cavity;
    GePolCavityH3Test() { SetUp(); }
    void SetUp() {
        Eigen::Vector3d H1( 0.735000, 0.000000, -1.333333);
        Eigen::Vector3d H2(-0.735000, 0.000000, -1.333333);
        Eigen::Vector3d H3( 0.000000, 0.000000,  2.666667);
        std::vector<Sphere> spheres;
        double radiusH = (1.20 * 1.20) / convertBohrToAngstrom;
        Sphere sph2(H1, radiusH);
        Sphere sph3(H2, radiusH);
        Sphere sph4(H3, radiusH);
        spheres.push_back(sph2);
        spheres.push_back(sph3);
        spheres.push_back(sph4);
        double area = 0.2 / convertBohr2ToAngstrom2;
        double probeRadius = 1.385 / convertBohrToAngstrom;
        double minRadius = 0.2 / convertBohrToAngstrom;
        // C1
        Symmetry pGroup = buildGroup(0, 0, 0, 0);
        cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
        cavity.saveCavity("h3+.npz");
    }
};

/*! \class GePolCavity
 *  \test \b GePolCavityH3Test_size tests GePol cavity size for H3+
 */
BOOST_FIXTURE_TEST_CASE(size, GePolCavityH3Test)
{
    int size = 312;
    int actualSize = cavity.size();
    BOOST_REQUIRE_EQUAL(size, actualSize);
}

/*! \class GePolCavity
 *  \test \b GePolCavityH3Test_area tests GePol cavity surface area for H3+
 */
BOOST_FIXTURE_TEST_CASE(area, GePolCavityH3Test)
{
    double area = 178.74700256125493;
    double actualArea = cavity.elementArea().sum();
    BOOST_REQUIRE_CLOSE(area, actualArea, 1.0e-10);
}

/*! \class GePolCavity
 *  \test \b GePolCavityH3Test_volume tests GePol cavity volume for H3+
 */
BOOST_FIXTURE_TEST_CASE(volume, GePolCavityH3Test)
{
    double volume = 196.4736029455637;
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
