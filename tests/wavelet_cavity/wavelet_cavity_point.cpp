#include <vector>
#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "PWCSolver.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "WaveletCavity.hpp"

// Disable obnoxious warnings from Google Test headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include "gtest/gtest.h"
#pragma GCC diagnostic pop
#endif

class WaveletCavityTest : public ::testing::Test
{
	protected:
		WaveletCavity cavity;
		virtual void SetUp()
		{
			Eigen::Vector3d N(0.0, 0.0, 0.0);
			std::vector<Sphere> spheres;
			Sphere sph1(N, 1.0);
			spheres.push_back(sph1);
			double probeRadius = 1.385; // Probe Radius for water
	        int patchLevel = 2;
        	double coarsity = 0.5;
			cavity = WaveletCavity(spheres, probeRadius, patchLevel, coarsity);
			cavity.readCavity("molec_dyadic.dat");
			double permittivity = 78.39;
			Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(); 
			UniformDielectric<AD_directional> * gfOutside = new UniformDielectric<AD_directional>(permittivity);
			int firstKind = 0;
			PWCSolver solver(gfInside, gfOutside, firstKind);
			solver.buildSystemMatrix(cavity);
			cavity.uploadPoints(solver.getQuadratureLevel(), solver.getT_(), false);
		}
};

TEST_F(WaveletCavityTest, size)
{
	int size = 4864;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(WaveletCavityTest, area)
{
	double area = 4.0 * M_PI * pow(1.0, 2);
 	double actualArea = cavity.elementArea().sum();
	EXPECT_DOUBLE_EQ(area, actualArea);
//	EXPECT_NEAR(area, actualArea, 1.0e-12);
}

TEST_F(WaveletCavityTest, volume)
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
