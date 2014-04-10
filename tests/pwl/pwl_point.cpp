#include <iostream>
#include <vector>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "PWLSolver.hpp"
#include "UniformDielectric.hpp"
#include "Vacuum.hpp"
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


TEST(PWLSolver, pointCharge) 
{
	// Set up cavity
	Eigen::Vector3d N(0.0, 0.0, 0.0); 	
	std::vector<Sphere> spheres;      	
	Sphere sph1(N, 2.929075493);      	
	spheres.push_back(sph1);
	double probeRadius = 1.385;
	int patchLevel = 2;
	double coarsity = 0.5;
	WaveletCavity cavity(spheres, probeRadius, patchLevel, coarsity);
	cavity.readCavity("molec_dyadic.dat");
	// The point charge is located at the origin.
	// The potential at cavity point s_I is Q/|s_I|
	double permittivity = 78.39;
	Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(); 
	UniformDielectric<AD_directional> * gfOutside = new UniformDielectric<AD_directional>(permittivity);
	int firstKind = 0;
	PWLSolver solver(gfInside, gfOutside, firstKind);
	solver.buildSystemMatrix(cavity);
	cavity.uploadPoints(solver.getQuadratureLevel(), solver.getT_(), true);

	double charge = 8.0;
	int size = cavity.size();
	Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
	for (int i = 0; i < size; ++i)
	{
		Eigen::Vector3d center = cavity.elementCenter(i);
		double distance = center.norm();
		fake_mep(i) = charge / distance; 
	}
	// The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
	Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
	solver.compCharge(fake_mep, fake_asc);
	double totalASC = - charge * (permittivity - 1) / permittivity;
	double totalFakeASC = fake_asc.sum();
	std::cout << "totalASC -totalFakeASC = " << totalASC - totalFakeASC << std::endl;
	EXPECT_NEAR(totalASC, totalFakeASC, 3e-3);
}