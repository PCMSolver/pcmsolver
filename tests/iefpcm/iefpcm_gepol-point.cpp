#include <iostream>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include "DerivativeTypes.hpp"
#include "GePolCavity.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "IEFSolver.hpp"
#include "Symmetry.hpp"

// Disable obnoxious warnings from Google Test headers
#if defined (__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall" 
#pragma GCC diagnostic ignored "-Weffc++" 
#pragma GCC diagnostic ignored "-Wextra"
#include "gtest/gtest.h"
#pragma GCC diagnostic pop
#endif

TEST(IEFSolver, pointChargeGePol) 
{
	// Set up cavity
	Eigen::Vector3d N(0.0, 0.0, 0.0); 		
	std::vector<Sphere> spheres;      		
	Sphere sph1(N, 2.929075493);      		
	spheres.push_back(sph1);          		
	double area = 0.4;
	double probeRadius = 0.0;
	double minRadius = 100.0;
	// C1 
	Symmetry group = buildGroup(0, 0, 0, 0);
	GePolCavity cavity = GePolCavity(spheres, area, probeRadius, minRadius, group);
	cavity.saveCavity("point.npz");
	
	double permittivity = 78.39;
	Vacuum<AD_directional> * gfInside = new Vacuum<AD_directional>(); 
	UniformDielectric<AD_directional> * gfOutside = new UniformDielectric<AD_directional>(permittivity);
	bool symm = true;
	IEFSolver solver(gfInside, gfOutside, symm);
	solver.buildSystemMatrix(cavity);

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
	std::cout << "totalASC - totalFakeASC = " << totalASC - totalFakeASC << std::endl;
	EXPECT_NEAR(totalASC, totalFakeASC, 3e-3);
}
