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

#include "TsLessCavity.hpp"

// Disable obnoxious warnings from Google Test headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include "gtest/gtest.h"
#pragma GCC diagnostic pop
#endif


class TsLessCavityTest : public ::testing::Test
{
	protected:
		TsLessCavity cav;
		virtual void SetUp()
		{
			Eigen::Vector3d origin(0.0, 0.0, 0.0);
			std::vector<Sphere> spheres;
			Sphere sph1(origin,  1.0);
			spheres.push_back(sph1);
			double area = 0.4;
			cav = TsLessCavity(spheres, area);
		}
};

TEST_F(TsLessCavityTest, size)
{
	int size = 32;
	int actualSize = cav.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(TsLessCavityTest, area)
{
	double area = 4.0 * M_PI * pow(1.0, 2);
 	double actualArea = cav.getElementArea().sum();
	EXPECT_DOUBLE_EQ(area, actualArea);
//	EXPECT_NEAR(area, actualArea, 1.0e-12);
}

TEST_F(TsLessCavityTest, volume)
{
	double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
	Eigen::Matrix3Xd elementCenter = cav.getElementCenter();
	Eigen::Matrix3Xd elementNormal = cav.getElementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cav.size(); ++i )
	{
		actualVolume += cav.getElementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_DOUBLE_EQ(volume, actualVolume);
//	EXPECT_NEAR(volume, actualVolume, 1.0e-12);
}
