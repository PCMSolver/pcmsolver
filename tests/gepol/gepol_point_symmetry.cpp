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

#include <boost/filesystem.hpp>

#include "GePolCavity.hpp"
#include "PhysicalConstants.hpp"

// Disable obnoxious warnings from Google Test headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include "gtest/gtest.h"
#pragma GCC diagnostic pop
#endif

namespace fs = boost::filesystem;

class GePolCavityC1Test : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
			Eigen::Vector3d origin(0.0, 0.0, 0.0);
			std::vector<Sphere> spheres;
			Sphere sph1(origin,  1.0);
			spheres.push_back(sph1);
			double area = 0.4;
			double probeRadius = 0.0;
			double minRadius = 100.0;
			int group = 0; // C1
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, group);
			fs::rename("PEDRA.OUT", "PEDRA.OUT.c1");
			fs::rename("cavity.off", "cavity.off.c1");
		}
};

TEST_F(GePolCavityC1Test, size)
{
	int size = 32;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC1Test, irreducible_size)
{
	int size = 32;
	int actualSize = cavity.irreducible_size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC1Test, area)
{
	double area = 4.0 * M_PI * pow(1.0, 2);
 	double actualArea = cavity.elementArea().sum();
	EXPECT_DOUBLE_EQ(area, actualArea);
//	EXPECT_NEAR(area, actualArea, 1.0e-12);
}

TEST_F(GePolCavityC1Test, volume)
{
	double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_DOUBLE_EQ(volume, actualVolume);
//	EXPECT_NEAR(volume, actualVolume, 1.0e-12);
}

class GePolCavityCsTest : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
			Eigen::Vector3d origin(0.0, 0.0, 0.0);
			std::vector<Sphere> spheres;
			Sphere sph1(origin,  1.0);
			spheres.push_back(sph1);
			double area = 0.4;
			double probeRadius = 0.0;
			double minRadius = 100.0;
			int group = 1; // Cs
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, group);
			fs::rename("PEDRA.OUT", "PEDRA.OUT.cs");
			fs::rename("cavity.off", "cavity.off.cs");
		}
};

TEST_F(GePolCavityCsTest, size)
{
	int size = 32;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityCsTest, irreducible_size)
{
	int size = 16;
	int actualSize = cavity.irreducible_size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityCsTest, area)
{
	double area = 4.0 * M_PI * pow(1.0, 2);
 	double actualArea = cavity.elementArea().sum();
	EXPECT_DOUBLE_EQ(area, actualArea);
//	EXPECT_NEAR(area, actualArea, 1.0e-12);
}

TEST_F(GePolCavityCsTest, volume)
{
	double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_DOUBLE_EQ(volume, actualVolume);
//	EXPECT_NEAR(volume, actualVolume, 1.0e-12);
}

class GePolCavityC2Test : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
			Eigen::Vector3d origin(0.0, 0.0, 0.0);
			std::vector<Sphere> spheres;
			Sphere sph1(origin,  1.0);
			spheres.push_back(sph1);
			double area = 0.4;
			double probeRadius = 0.0;
			double minRadius = 100.0;
			int group = 2; // C2
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, group);
			fs::rename("PEDRA.OUT", "PEDRA.OUT.c2");
			fs::rename("cavity.off", "cavity.off.c2");
		}
};

TEST_F(GePolCavityC2Test, size)
{
	int size = 32;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC2Test, irreducible_size)
{
	int size = 16;
	int actualSize = cavity.irreducible_size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC2Test, area)
{
	double area = 4.0 * M_PI * pow(1.0, 2);
 	double actualArea = cavity.elementArea().sum();
	EXPECT_DOUBLE_EQ(area, actualArea);
//	EXPECT_NEAR(area, actualArea, 1.0e-12);
}

TEST_F(GePolCavityC2Test, volume)
{
	double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_DOUBLE_EQ(volume, actualVolume);
//	EXPECT_NEAR(volume, actualVolume, 1.0e-12);
}

class GePolCavityCiTest : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
			Eigen::Vector3d origin(0.0, 0.0, 0.0);
			std::vector<Sphere> spheres;
			Sphere sph1(origin,  1.0);
			spheres.push_back(sph1);
			double area = 0.4;
			double probeRadius = 0.0;
			double minRadius = 100.0;
			int group = 3; // Ci
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, group);
			fs::rename("PEDRA.OUT", "PEDRA.OUT.ci");
			fs::rename("cavity.off", "cavity.off.ci");
		}
};

TEST_F(GePolCavityCiTest, size)
{
	int size = 32;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityCiTest, irreducible_size)
{
	int size = 16;
	int actualSize = cavity.irreducible_size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityCiTest, area)
{
	double area = 4.0 * M_PI * pow(1.0, 2);
 	double actualArea = cavity.elementArea().sum();
	EXPECT_DOUBLE_EQ(area, actualArea);
//	EXPECT_NEAR(area, actualArea, 1.0e-12);
}

TEST_F(GePolCavityCiTest, volume)
{
	double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_DOUBLE_EQ(volume, actualVolume);
//	EXPECT_NEAR(volume, actualVolume, 1.0e-12);
}

class GePolCavityC2hTest : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
			Eigen::Vector3d origin(0.0, 0.0, 0.0);
			std::vector<Sphere> spheres;
			Sphere sph1(origin,  1.0);
			spheres.push_back(sph1);
			double area = 0.4;
			double probeRadius = 0.0;
			double minRadius = 100.0;
			int group = 4; // C2h
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, group);
			fs::rename("PEDRA.OUT", "PEDRA.OUT.c2h");
			fs::rename("cavity.off", "cavity.off.c2h");
		}
};

TEST_F(GePolCavityC2hTest, size)
{
	int size = 32;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC2hTest, irreducible_size)
{
	int size = 8;
	int actualSize = cavity.irreducible_size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC2hTest, area)
{
	double area = 4.0 * M_PI * pow(1.0, 2);
 	double actualArea = cavity.elementArea().sum();
	EXPECT_DOUBLE_EQ(area, actualArea);
//	EXPECT_NEAR(area, actualArea, 1.0e-12);
}

TEST_F(GePolCavityC2hTest, volume)
{
	double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_DOUBLE_EQ(volume, actualVolume);
//	EXPECT_NEAR(volume, actualVolume, 1.0e-12);
}

class GePolCavityD2Test : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
			Eigen::Vector3d origin(0.0, 0.0, 0.0);
			std::vector<Sphere> spheres;
			Sphere sph1(origin,  1.0);
			spheres.push_back(sph1);
			double area = 0.4;
			double probeRadius = 0.0;
			double minRadius = 100.0;
			int group = 5; // D2
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, group);
			fs::rename("PEDRA.OUT", "PEDRA.OUT.d2");
			fs::rename("cavity.off", "cavity.off.d2");
		}
};

TEST_F(GePolCavityD2Test, size)
{
	int size = 32;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityD2Test, irreducible_size)
{
	int size = 8;
	int actualSize = cavity.irreducible_size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityD2Test, area)
{
	double area = 4.0 * M_PI * pow(1.0, 2);
 	double actualArea = cavity.elementArea().sum();
	EXPECT_DOUBLE_EQ(area, actualArea);
//	EXPECT_NEAR(area, actualArea, 1.0e-12);
}

TEST_F(GePolCavityD2Test, volume)
{
	double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_DOUBLE_EQ(volume, actualVolume);
//	EXPECT_NEAR(volume, actualVolume, 1.0e-12);
}

class GePolCavityC2vTest : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
			Eigen::Vector3d origin(0.0, 0.0, 0.0);
			std::vector<Sphere> spheres;
			Sphere sph1(origin,  1.0);
			spheres.push_back(sph1);
			double area = 0.4;
			double probeRadius = 0.0;
			double minRadius = 100.0;
			int group = 6; // C2v
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, group);
			fs::rename("PEDRA.OUT", "PEDRA.OUT.c2v");
			fs::rename("cavity.off", "cavity.off.c2v");
		}
};

TEST_F(GePolCavityC2vTest, size)
{
	int size = 32;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC2vTest, irreducible_size)
{
	int size = 8;
	int actualSize = cavity.irreducible_size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC2vTest, area)
{
	double area = 4.0 * M_PI * pow(1.0, 2);
 	double actualArea = cavity.elementArea().sum();
	EXPECT_DOUBLE_EQ(area, actualArea);
//	EXPECT_NEAR(area, actualArea, 1.0e-12);
}

TEST_F(GePolCavityC2vTest, volume)
{
	double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_DOUBLE_EQ(volume, actualVolume);
//	EXPECT_NEAR(volume, actualVolume, 1.0e-12);
}

class GePolCavityD2hTest : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
			Eigen::Vector3d origin(0.0, 0.0, 0.0);
			std::vector<Sphere> spheres;
			Sphere sph1(origin,  1.0);
			spheres.push_back(sph1);
			double area = 0.4;
			double probeRadius = 0.0;
			double minRadius = 100.0;
			int group = 7; // D2h
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, group);
			fs::rename("PEDRA.OUT", "PEDRA.OUT.d2h");
			fs::rename("cavity.off", "cavity.off.d2h");
		}
};

TEST_F(GePolCavityD2hTest, size)
{
	int size = 32;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityD2hTest, irreducible_size)
{
	int size = 4;
	int actualSize = cavity.irreducible_size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityD2hTest, area)
{
	double area = 4.0 * M_PI * pow(1.0, 2);
 	double actualArea = cavity.elementArea().sum();
	EXPECT_DOUBLE_EQ(area, actualArea);
//	EXPECT_NEAR(area, actualArea, 1.0e-12);
}

TEST_F(GePolCavityD2hTest, volume)
{
	double volume = 4.0 * M_PI * pow(1.0, 3) / 3.0;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_DOUBLE_EQ(volume, actualVolume);
//	EXPECT_NEAR(volume, actualVolume, 1.0e-12);
}
