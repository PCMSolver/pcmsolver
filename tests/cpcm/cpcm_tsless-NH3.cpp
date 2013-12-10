#include <iostream>

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
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "CPCMSolver.hpp"

// Disable obnoxious warnings from Google Test headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include "gtest/gtest.h"
#pragma GCC diagnostic pop
#endif

TEST(CPCMSolver, NH3TsLess) 
{
	// Set up cavity
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
	TsLessCavity cavity(spheres, area);		
	
	double permittivity = 78.39;
	Vacuum * gfInside = new Vacuum(2); // Automatic directional derivative
	UniformDielectric * gfOutside = new UniformDielectric(2, permittivity);
	double correction = 0.0;
	CPCMSolver solver(gfInside, gfOutside, correction);
	solver.buildSystemMatrix(cavity);

	double Ncharge = 7.0;
	double Hcharge = 1.0;
	int size = cavity.size();
	Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
	for (int i = 0; i < size; ++i)
	{
		Eigen::Vector3d center = cavity.getElementCenter(i);
		double Ndistance = (center - N).norm();
		double H1distance = (center - H1).norm();
		double H2distance = (center - H2).norm();
		double H3distance = (center - H3).norm();
		fake_mep(i) = Ncharge / Ndistance + Hcharge / H1distance + Hcharge / H2distance + Hcharge / H3distance; 
	}
	// The total ASC for a conductor is -Q
	// for CPCM it will be -Q*[(epsilon-1)/epsilon}
	Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
	solver.compCharge(fake_mep, fake_asc);
	double totalASC = - (Ncharge + 3.0 * Hcharge) * (permittivity - 1) / permittivity;
	double totalFakeASC = fake_asc.sum();
	std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
	EXPECT_NEAR(totalASC, totalFakeASC, 3e-3);
}
