#include <vector>
#include <cmath>

#include "Config.hpp"

// Disable obnoxious warnings from Eigen headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#elif (__INTEL_COMPILER)
#pragma warning push
#pragma warning disable "-Wall"
#include <Eigen/Dense>
#pragma warning pop
#endif

#include "GePolCavity.hpp"

// Disable obnoxious warnings from Google Test headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include "gtest/gtest.h"
#pragma GCC diagnostic pop
#endif

class GePolCavityNH3Test : public ::testing::Test
{
	protected:
		GePolCavity cavity;
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
			cavity = GePolCavity(spheres, area);
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
 	double actualArea = cavity.getElementArea().sum();
	EXPECT_NEAR(area, actualArea, 1.0e-10);
}

TEST_F(GePolCavityNH3Test, volume)
{
	double volume = 152.81441857040116;
	Eigen::Matrix3Xd elementCenter = cavity.getElementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.getElementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.getElementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_NEAR(volume, actualVolume, 1.0e-10);
}
