#include <vector>
#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "GePolCavity.hpp"

#include "gtestPimpl.hpp"

class GePolCavityH3Test : public ::testing::Test
{
protected:
    GePolCavity cavity;
    virtual void SetUp() {
        cavity.loadCavity("h3+.npz");
    }
};

TEST_F(GePolCavityH3Test, size)
{
    int size = 312;
    int actualSize = cavity.size();
    EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityH3Test, area)
{
    double area = 178.74700256125493;
    double actualArea = cavity.elementArea().sum();
    EXPECT_NEAR(area, actualArea, 1.0e-10);
}

TEST_F(GePolCavityH3Test, volume)
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
