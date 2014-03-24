#include <vector>
#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"

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

class WaveletCavityNH3Test : public ::testing::Test
{
	protected:
		WaveletCavity cavity;
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
	                int patchLevel = 2;
        	        double coarsity = 0.5;
			cavity = WaveletCavity(spheres, probeRadius, patchLevel, coarsity);
			cavity.readCavity("molec_dyadic.dat");
			double permittivity = 78.39;
			Vacuum * gfInside = new Vacuum(2); // Automatic directional derivative
			UniformDielectric * gfOutside = new UniformDielectric(2, permittivity);
			int firstKind = 0;
			PWCSolver solver(gfInside, gfOutside, firstKind);
			solver.buildSystemMatrix(cavity);
			cavity.uploadPoints(solver.getQuadratureLevel(), solver.getT_(), false);
		}
};

TEST_F(WaveletCavityNH3Test, size)
{
	int size = 4288;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(WaveletCavityNH3Test, area)
{
	double area = 146.41490284471513;
 	double actualArea = cavity.elementArea().sum();
	EXPECT_NEAR(area, actualArea, 1.0e-10);
}

TEST_F(WaveletCavityNH3Test, volume)
{
	double volume = 153.3346491517182;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_NEAR(volume, actualVolume, 1.0e-10);
}
