#include <iostream>

#include <Eigen/Dense>

#include "GePolCavity.h"
#include "PCMSolver.h"
#include "CPCMSolver.h"
#include "gtest/gtest.h"

TEST(cpcm, point_charge) {
	Eigen::Vector3d fuffa;
	fuffa << 1.0, 1.0, 1.0;
	Eigen::Vector3d raboof;
	raboof << 1.0, 1.0, 1.0;
	for (int i = 0; i < 3; ++i) {
		EXPECT_EQ(fuffa(i), raboof(i));
	}
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
