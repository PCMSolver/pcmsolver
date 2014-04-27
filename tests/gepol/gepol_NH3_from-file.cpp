#include <vector>
#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "GePolCavity.hpp"

#include "gtestPimpl.hpp"

class GePolCavityNH3Test : public ::testing::Test
{
protected:
    GePolCavity cavity;
    virtual void SetUp() {
        cavity.loadCavity("nh3.npz");
    }
};

TEST_F(GePolCavityNH3Test, size)
{
    int size = 544;
    int actualSize = cavity.size();
    EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityNH3Test, area)
{
    double area = 147.18581691164593;
    double actualArea = cavity.elementArea().sum();
    EXPECT_NEAR(area, actualArea, 1.0e-10);
}

TEST_F(GePolCavityNH3Test, volume)
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
