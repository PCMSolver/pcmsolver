#include "catch.hpp"

#include <Eigen/Core>

#include "cavity/TsLessCavity.hpp"
#include "TestingMolecules.hpp"

using namespace pcm;
using cavity::TsLessCavity;

TEST_CASE("TsLess cavity for a single sphere", "[tsless][tsless_point]")
{
    Molecule point = dummy<0>();
    double area = 0.4;
    double minDistance = 0.1;
    double probeRadius = 0.0;
    int derOrder = 4;
    TsLessCavity cavity(point, area, probeRadius, 100.0, minDistance, derOrder);

    /*! \class TsLessCavity
     *  \test \b TsLessCavityTest_size tests TsLess cavity size for a point charge
     */
    SECTION("Test size")
    {
        size_t ref_size = 32;
        size_t size = cavity.size();
        REQUIRE(ref_size == size);
    }

    /*! \class TsLessCavity
     *  \test \b TsLessCavityTest_area tests TsLess cavity surface area for a point charge
     */
    SECTION("Test surface area")
    {
        double ref_area = 4.0 * M_PI * pow(1.0, 2);
        double area = cavity.elementArea().sum();
        REQUIRE(ref_area == Approx(area));
    }

    /*! \class TsLessCavity
     *  \test \b TsLessCavityTest_volume tests TsLess cavity volume for a point charge
     */
    SECTION("Test volume")
    {
        double ref_volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
        Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
        Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
        double volume = 0;
        for ( size_t i = 0; i < cavity.size(); ++i ) {
            volume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
        }
        volume /= 3;
        REQUIRE(ref_volume == Approx(volume));
    }
}
