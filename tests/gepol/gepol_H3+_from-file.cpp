#include <vector>
#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "GePolCavity.hpp"

#include "gtestPimpl.hpp"

class GePolCavityH3RestartTest : public ::testing::Test
{
protected:
    GePolCavity cavity;
    virtual void SetUp() {
        cavity.loadCavity("h3+.npz");
    }
};

/*! \class GePolCavity 
 *  \test \b GePolCavityH3RestartTest_size tests GePol cavity size for H3+ loading the cavity from a .npz file
 */
TEST_F(GePolCavityH3RestartTest, size)
{
    int size = 312;
    int actualSize = cavity.size();
    EXPECT_EQ(size, actualSize);
}

/*! \class GePolCavity 
 *  \test \b GePolCavityH3RestartTest_area tests GePol cavity surface area for H3+ loading the cavity from from a .npz file
 */
TEST_F(GePolCavityH3RestartTest, area)
{
    double area = 178.74700256125493;
    double actualArea = cavity.elementArea().sum();
    EXPECT_NEAR(area, actualArea, 1.0e-10);
}

/*! \class GePolCavity 
 *  \test \b GePolCavityH3RestartTest_volume tests GePol cavity volume for H3+ loading the cavity from from a .npz file
 */
TEST_F(GePolCavityH3RestartTest, volume)
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
    EXPECT_NEAR(volume, actualVolume, 1.0e-10);
}
