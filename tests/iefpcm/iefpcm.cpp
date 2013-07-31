#include <iostream>

#include <Eigen/Dense>

#include "GePolCavity.hpp"
#include "TsLessCavity.hpp"
#include "Vacuum.hpp"
#include "UniformDielectric.hpp"
#include "IEFSolver.hpp"
#include "gtest/gtest.h"

TEST(IEFSolver, pointChargeGePol) 
{
	// Set up cavity
	Eigen::Vector3d N(0.0, 0.0, 0.0); 		
	std::vector<Sphere> spheres;      		
	Sphere sph1(N, 2.929075493);      		
	spheres.push_back(sph1);          		
	double area = 0.4;                		
	GePolCavity cavity(spheres, area);		
	// The point charge is located at the origin.
	// The potential at cavity point s_I is Q/|s_I|
	double permittivity = 78.39;
	Vacuum * gfInside = new Vacuum(2); // Automatic directional derivative
	UniformDielectric * gfOutside = new UniformDielectric(2, permittivity);
	IEFSolver solver(gfInside, gfOutside);
	solver.buildSystemMatrix(cavity);

	double charge = 8.0;
	int size = cavity.size();
	Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
	for (int i = 0; i < size; ++i)
	{
		Eigen::Vector3d center = cavity.getElementCenter(i);
		double distance = center.norm();
		fake_mep(i) = charge / distance; 
	}
	// The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
	Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
	solver.compCharge(fake_mep, fake_asc);
	double totalASC = - charge * (permittivity - 1) / permittivity;
	double totalFakeASC = fake_asc.sum();
	EXPECT_NEAR(totalASC, totalFakeASC, 3e-3);
}

TEST(IEFSolver, pointChargeTsLess) 
{
	// Set up cavity
	Eigen::Vector3d N(0.0, 0.0, 0.0); 		
	std::vector<Sphere> spheres;      		
	Sphere sph1(N, 2.929075493);      		
	spheres.push_back(sph1);          		
	double area = 0.4;                		
	TsLessCavity cavity(spheres, area);		
	// The point charge is located at the origin.
	// The potential at cavity point s_I is Q/|s_I|
	double permittivity = 78.39;
	Vacuum * gfInside = new Vacuum(2); // Automatic directional derivative
	UniformDielectric * gfOutside = new UniformDielectric(2, permittivity);
	IEFSolver solver(gfInside, gfOutside);
	solver.buildSystemMatrix(cavity);

	double charge = 8.0;
	int size = cavity.size();
	Eigen::VectorXd fake_mep = Eigen::VectorXd::Zero(size);
	for (int i = 0; i < size; ++i)
	{
		Eigen::Vector3d center = cavity.getElementCenter(i);
		double distance = center.norm();
		fake_mep(i) = charge / distance; 
	}
	// The total ASC for a dielectric is -Q*[(epsilon-1)/epsilon]
	Eigen::VectorXd fake_asc = Eigen::VectorXd::Zero(size);
	solver.compCharge(fake_mep, fake_asc);
	double totalASC = - charge * (permittivity - 1) / permittivity;
	double totalFakeASC = fake_asc.sum();
	EXPECT_NEAR(totalASC, totalFakeASC, 3e-3);
}
