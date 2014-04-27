#include <vector>
#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "TsLessCavity.hpp"

#include "gtestPimpl.hpp"

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
			double minDistance = 0.1;
			double probeRadius = 0.0;
			int derOrder = 4;
			cav = TsLessCavity(spheres, area, probeRadius, minDistance, derOrder);
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
 	double actualArea = cav.elementArea().sum();
	EXPECT_DOUBLE_EQ(area, actualArea);
//	EXPECT_NEAR(area, actualArea, 1.0e-12);
}

TEST_F(TsLessCavityTest, volume)
{
	double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
	Eigen::Matrix3Xd elementCenter = cav.elementCenter();
	Eigen::Matrix3Xd elementNormal = cav.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cav.size(); ++i )
	{
		actualVolume += cav.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_DOUBLE_EQ(volume, actualVolume);
//	EXPECT_NEAR(volume, actualVolume, 1.0e-12);
}
