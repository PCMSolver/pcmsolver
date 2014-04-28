#include <vector>
#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "GePolCavity.hpp"

#include "gtestPimpl.hpp"

class GePolCavityNH3RestartTest : public ::testing::Test
{
protected:
    GePolCavity cavity;
    virtual void SetUp() {
        cavity.loadCavity("nh3.npz");
    }
};

/*! \class GePolCavity 
 *  \test \b GePolCavityNH3RestartTest_size tests GePol cavity size for ammonia loading the cavity from a .npz file
 */
TEST_F(GePolCavityNH3RestartTest, size)
{
    int size = 544;
    int actualSize = cavity.size();
    EXPECT_EQ(size, actualSize);
}

/*! \class GePolCavity 
 *  \test \b GePolCavityNH3RestartTest_area tests GePol cavity surface area for ammonia loading the cavity from from a .npz file
 */
TEST_F(GePolCavityNH3RestartTest, area)
{
    double area = 147.18581691164593;
    double actualArea = cavity.elementArea().sum();
    EXPECT_NEAR(area, actualArea, 1.0e-10);
}

/*! \class GePolCavity 
 *  \test \b GePolCavityNH3RestartTest_volume tests GePol cavity volume for ammonia loading the cavity from from a .npz file
 */
TEST_F(GePolCavityNH3RestartTest, volume)
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
    EXPECT_NEAR(volume, actualVolume, 1.0e-10);
}
