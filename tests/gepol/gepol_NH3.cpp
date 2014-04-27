#include <vector>
#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "GePolCavity.hpp"
#include "PhysicalConstants.hpp"
#include "Symmetry.hpp"

#include "gtestPimpl.hpp"

class GePolCavityNH3Test : public ::testing::Test
{
protected:
    GePolCavity cavity;
    virtual void SetUp() {
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
        double area = 0.4; // Bohr^2
        double probeRadius = 1.385 / convertBohrToAngstrom;
        double minRadius = 0.2 / convertBohrToAngstrom;
        // C1
        Symmetry pGroup = buildGroup(0, 0, 0, 0);
        cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
        cavity.saveCavity("nh3.npz");
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
