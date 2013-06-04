#include <Eigen/Dense>

#include <vector>
#include <cmath>
#include "GePolCavity.h"
#include "gtest/gtest.h"

class GePolCavityTest : public ::testing::Test
{
	protected:
		GePolCavity cav;
		virtual void SetUp()
		{
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
			cav = GePolCavity(spheres, area);
		}
};

TEST_F(GePolCavityTest, size)
{
	int size = 544;
	int actualSize = cav.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityTest, area)
{
	double area = 147.18581691164593;
 	double actualArea = cav.getElementArea().sum();
	EXPECT_DOUBLE_EQ(area, actualArea);
}

TEST_F(GePolCavityTest, volume)
{
	double volume = 152.81441857040116;
	Eigen::Matrix3Xd elementCenter = cav.getElementCenter();
	Eigen::Matrix3Xd elementNormal = cav.getElementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cav.size(); ++i )
	{
		actualVolume += cav.getElementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_DOUBLE_EQ(volume, actualVolume);
}
