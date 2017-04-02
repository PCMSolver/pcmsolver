#include "catch.hpp"

#include <Eigen/Core>

#include "cavity/TsLessCavity.hpp"
#include "TestingMolecules.hpp"

using namespace pcm;
using cavity::TsLessCavity;

TEST_CASE("TsLess cavity for an hydrogen fluoride molecule", "[tsless][tsless_HF]")
{
    Molecule molec = HF();
    double area = 0.02 / bohr2ToAngstrom2();
    double probeRadius = 0.0 / bohrToAngstrom();
    double minRadius = 0.2 / bohrToAngstrom();
    double minDistance = 0.01;
    int derOrder = 4;
    TsLessCavity cavity(molec, area, probeRadius, minRadius, minDistance, derOrder);

    /*! \class TsLessCavity
     *  \test \b TsLessCavityHFTest_size tests TsLess cavity size for ammonia
     */
    SECTION("Test size")
    {
        size_t ref_size = 1498;
        size_t size = cavity.size();
        REQUIRE(size == ref_size);
    }

    /*! \class TsLessCavity
     *  \test \b TsLessCavityHFTest_area tests TsLess cavity surface area for ammonia
     */
    SECTION("Test surface area")
    {
        double ref_area = 111.07794350824591;
        double area = cavity.elementArea().sum();
        REQUIRE(area == Approx(ref_area));
    }

    /*! \class TsLessCavity
     *  \test \b TsLessCavityHFTest_volume tests TsLess cavity volume for ammonia
     */
    SECTION("Test volume")
    {
        double ref_volume = 106.59560252994531;
        Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
        Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
        double volume = 0;
        for ( size_t i = 0; i < cavity.size(); ++i ) {
            volume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
        }
        volume /= 3;
        REQUIRE(volume == Approx(ref_volume));
    }
}
