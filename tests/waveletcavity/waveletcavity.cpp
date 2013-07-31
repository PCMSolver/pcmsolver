#include <Eigen/Dense>

#include <vector>
#include <cmath>
#include "PWCSolver.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "WaveletCavity.hpp"
#include "gtest/gtest.h"

class WaveletCavityTest : public ::testing::Test
{
	protected:
		WaveletCavity cav;
		virtual void SetUp()
		{
			Eigen::Vector3d N(0.0, 0.0, 0.0);
			std::vector<Sphere> spheres;
			Sphere sph1(N, 1.0);
			spheres.push_back(sph1);
			double probeRadius = 1.385; // Probe Radius for water
			cav = WaveletCavity(spheres, probeRadius);
			cav.readCavity("molec_dyadic.dat");
			double permittivity = 78.39;
			Vacuum * gfInside = new Vacuum(2); // Automatic directional derivative
			UniformDielectric * gfOutside = new UniformDielectric(2, permittivity);
			int firstKind = 0;
			PWCSolver solver(gfInside, gfOutside, firstKind);
			solver.buildSystemMatrix(cav);
			cav.uploadPoints(solver.getQuadratureLevel(), solver.getT_(), false);
		}
};

TEST_F(WaveletCavityTest, size)
{
	int size = 4864;
	int actualSize = cav.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(WaveletCavityTest, area)
{
	double area = 4.0 * M_PI * pow(1.0, 2);
 	double actualArea = cav.getElementArea().sum();
	EXPECT_DOUBLE_EQ(area, actualArea);
//	EXPECT_NEAR(area, actualArea, 1.0e-12);
}

TEST_F(WaveletCavityTest, volume)
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

class WaveletCavityNH3Test : public ::testing::Test
{
	protected:
		WaveletCavity cav;
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
			double probeRadius = 1.385; // Probe Radius for water
			cav = WaveletCavity(spheres, probeRadius);
			cav.readCavity("molec_dyadic.dat");
			double permittivity = 78.39;
			Vacuum * gfInside = new Vacuum(2); // Automatic directional derivative
			UniformDielectric * gfOutside = new UniformDielectric(2, permittivity);
			int firstKind = 0;
			PWCSolver solver(gfInside, gfOutside, firstKind);
			solver.buildSystemMatrix(cav);
			cav.uploadPoints(solver.getQuadratureLevel(), solver.getT_(), false);
		}
};

TEST_F(WaveletCavityNH3Test, size)
{
	int size = 4864;
	int actualSize = cav.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(WaveletCavityNH3Test, area)
{
	double area = 146.41602862089889;
 	double actualArea = cav.getElementArea().sum();
	EXPECT_NEAR(area, actualArea, 1.0e-12);
}

TEST_F(WaveletCavityNH3Test, volume)
{
	double volume = 153.3130201033002;
	Eigen::Matrix3Xd elementCenter = cav.getElementCenter();
	Eigen::Matrix3Xd elementNormal = cav.getElementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cav.size(); ++i )
	{
		actualVolume += cav.getElementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_NEAR(volume, actualVolume, 1.0e-12);
}
