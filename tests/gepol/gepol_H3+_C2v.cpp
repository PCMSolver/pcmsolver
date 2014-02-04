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

// Test C2v symmetry with addition of extra spheres enabled
class GePolCavityC2vAddTest : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
			Eigen::Vector3d H1( 0.735000, 0.000000, -1.333333);
			Eigen::Vector3d H2( 0.000000, 0.000000,  2.666667);
			Eigen::Vector3d H3(-0.735000, 0.000000, -1.333333);
			std::vector<Sphere> spheres;
			double radiusH = (1.20 * 1.20) / convertBohrToAngstrom;
			Sphere sph1(H1, radiusH);
			Sphere sph2(H2, radiusH);
			Sphere sph3(H3, radiusH);
			spheres.push_back(sph1);
			spheres.push_back(sph2);
			spheres.push_back(sph3);
			double area = 0.2 / convertBohr2ToAngstrom2;
			double probeRadius = 1.385 / convertBohrToAngstrom;
			double minRadius = 0.2 / convertBohrToAngstrom;
			int pGroup = 6; // C2v
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
			cavity.saveCavity("h3+_c2v.npz");
			fs::rename("PEDRA.OUT", "PEDRA.OUT.c2v");
			fs::rename("cavity.off", "cavity.off.c2v");
		}
};

TEST_F(GePolCavityC2vAddTest, size)
{
	int size = 312;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC2vAddTest, area)
{
	double area = 178.74700256128352;
 	double actualArea = cavity.getElementArea().sum();
	EXPECT_NEAR(area, actualArea, 1.0e-10);
}

TEST_F(GePolCavityC2vAddTest, volume)
{
	double volume = 196.47360294559090;
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

// Test C2v symmetry without addition of extra spheres enabled
class GePolCavityC2vTest : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
			Eigen::Vector3d H1( 0.735000, 0.000000, -1.333333);
			Eigen::Vector3d H2( 0.000000, 0.000000,  2.666667);
			Eigen::Vector3d H3(-0.735000, 0.000000, -1.333333);
			std::vector<Sphere> spheres;
			double radiusH = (1.20 * 1.20) / convertBohrToAngstrom;
			Sphere sph1(H1, radiusH);
			Sphere sph2(H2, radiusH);
			Sphere sph3(H3, radiusH);
			spheres.push_back(sph1);
			spheres.push_back(sph2);
			spheres.push_back(sph3);
			double area = 0.2 / convertBohr2ToAngstrom2;
			double probeRadius = 1.385 / convertBohrToAngstrom;
			double minRadius = 100.0 / convertBohrToAngstrom;
			int pGroup = 6; // C2v
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
			cavity.saveCavity("h3+_c2v_noadd.npz");
			fs::rename("PEDRA.OUT", "PEDRA.OUT.c2v_noadd");
			fs::rename("cavity.off", "cavity.off.c2v_noadd");
		}
};

TEST_F(GePolCavityC2vTest, size)
{
	int size = 288;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC2vTest, area)
{
	double area = 181.87043332808548;
 	double actualArea = cavity.getElementArea().sum();
	EXPECT_NEAR(area, actualArea, 1.0e-10);
}

TEST_F(GePolCavityC2vTest, volume)
{
	double volume = 192.48281460140359;
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
