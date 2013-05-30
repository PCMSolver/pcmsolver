#include <iostream>

#include <Eigen/Dense>

#include "Getkw.h"
#include "GePolCavity.h"
#include "PCMSolver.h"
#include "IEFSolver.h"
//#include "Interface.h"
#include "gtest/gtest.h"

TEST(iefpcm, pointCharge) {
// This vector contains the potential of a point charge q at the sphere surface.
	Eigen::Vector3d potential;
// The exact total surface charge is given by (eps - 1)/eps * q
	double totalCharge;
	Eigen::Vector3d origin;
	origin << 0.0, 0.0, 0.0;
	double radius = 1.0;
	Sphere point(origin, radius);
        std::vector<Sphere> cavSpheres;
	cavSpheres.push_back(point);
	const char *infile = "@pcmsolver.inp";
	Getkw Input = Getkw(infile, false, true);
////////string solvent = Input.getStr("Medium.Solvent");
	Section cavSect = Input.getSect("Cavity<gepol>");
//        GePolCavity cavity(cavSect, cavSpheres);		
// Read a dummy input that has already been parsed.
//	init_pcm_();
//  std::string name("foobar");
//  Eigen::MatrixX3d geometry;
//  geometry.resize(2, Eigen::NoChange);
//  geometry << 0.0, 0.0, 0.0,
//              1.0, 0.0, 0.0;
//  Eigen::VectorXd charges;
//  charges.resize(2);
//  charges << 1.0, 1.0;
//  Molecule foobar(charges, geometry, name);
//  foobar.moveToCOM();

//  Eigen::Matrix3d inertia;
//  foobar.moveToPAF();
//  inertia = foobar.getInertiaTensor();
//  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(inertia);
//  
//  Eigen::Vector3d moments;
//  moments << 0.0, 0.50397, 0.50397; 

//  for(int i = 0; i < 3; ++i) {
//      EXPECT_EQ(eigenSolver.eigenvalues()[i], moments[i]);
//  }
}
