/*! \file IEFSolver.cpp 
\brief PCM solver
*/

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "Getkw.h"
#include "taylor.hpp"
#include "GreensFunctionInterface.h"
#include "Cavity.h"
#include "GePolCavity.h"
#include "PCMSolver.h"
#include "IEFSolver.h"

IEFSolver::IEFSolver(GreensFunctionInterface & gfi, GreensFunctionInterface & gfo) : 
	PCMSolver(gfi, gfo) {
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
}

IEFSolver::IEFSolver(GreensFunctionInterface * gfi, GreensFunctionInterface * gfo) :
	PCMSolver(gfi, gfo) {
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
}

IEFSolver::IEFSolver(Section solver) : PCMSolver(solver) {
	builtIsotropicMatrix = false;
	builtAnisotropicMatrix = false;
}

IEFSolver::~IEFSolver(){
}

void IEFSolver::buildSystemMatrix(Cavity & cavity) {
    if (GePolCavity *gePolCavity = dynamic_cast<GePolCavity*> (&cavity)) {
		if (this->greenInside->isUniform() && this->greenOutside->isUniform()) {
			buildIsotropicMatrix(*gePolCavity);
		} else {
			buildAnisotropicMatrix(*gePolCavity);
		}
	} else {
		std::cout << "No other cavity than GePol for traditional PCM" << std::endl;
		exit(-1);
	}
}

void IEFSolver::buildAnisotropicMatrix(GePolCavity & cav){
    cavitySize = cav.size();
    MatrixXd SI(cavitySize, cavitySize);
    MatrixXd SE(cavitySize, cavitySize);
    MatrixXd DI(cavitySize, cavitySize);
    MatrixXd DE(cavitySize, cavitySize);
    for(int i = 0; i < cavitySize; i++){
		Vector3d p1 = cav.getTessCenter(i);
		Vector3d n1 = cav.getTessNormal(i);
		double area = cav.getTessArea(i);
		double radius = cav.getTessRadius(i);
		SI(i,i) =  greenInside->compDiagonalElementS(area); 
		SE(i,i) = greenOutside->compDiagonalElementS(area); 
		DI(i,i) =  greenInside->compDiagonalElementD(area, radius); 
		DE(i,i) = greenOutside->compDiagonalElementD(area, radius); 
		for (int j = 0; j < cavitySize; j++){
			Vector3d p2 = cav.getTessCenter(j);
			Vector3d n2 = cav.getTessNormal(j);
			if (i != j) {
				SI(i,j) = greenInside->evalf(p1, p2);
				SE(i,j) = greenOutside->evalf(p1, p2);
				DI(i,j) = greenInside->evald(n2, p1, p2);
				DE(i,j) = greenOutside->evald(n2, p1, p2);
			}
		}
    }
    MatrixXd a(cavitySize, cavitySize);
    MatrixXd aInv(cavitySize, cavitySize);
    a.setZero();
    aInv.setZero();
    for (int i = 0; i < cavitySize; i++) {
		a(i,i) = cav.getTessArea(i);
		aInv(i,i) = 2 * M_PI / cav.getTessArea(i);
    }
    PCMMatrix = ((aInv - DE) * a * SI + SE * a * (aInv + DI.transpose()));
    PCMMatrix = PCMMatrix.inverse();
    PCMMatrix *= ((aInv - DE) - SE * SI.inverse() * (aInv - DI));
    PCMMatrix = PCMMatrix * a;
	builtAnisotropicMatrix = true;
	builtIsotropicMatrix = false;
}


void IEFSolver::buildIsotropicMatrix(GePolCavity & cav){
	double epsilon = this->greenOutside->getDielectricConstant();
    cavitySize = cav.size();
    MatrixXd SI(cavitySize, cavitySize);
    MatrixXd DI(cavitySize, cavitySize);
    for(int i = 0; i < cavitySize; i++){
		Vector3d p1 = cav.getTessCenter(i);
		Vector3d n1 = cav.getTessNormal(i);
		double area = cav.getTessArea(i);
		double radius = cav.getTessRadius(i);
		SI(i,i) = greenInside->compDiagonalElementS(area); 
		DI(i,i) = greenInside->compDiagonalElementD(area, radius); 
		for (int j = 0; j < cavitySize; j++){
			Vector3d p2 = cav.getTessCenter(j);
			Vector3d n2 = cav.getTessNormal(j);
			if (i != j) {
				SI(i,j) = greenInside->evalf(p1, p2);
				DI(i,j) = greenInside->evald(n2, p1, p2);
			}
		}
    }
    MatrixXd a(cavitySize, cavitySize);
    MatrixXd aInv(cavitySize, cavitySize);
    a.setZero();
    aInv.setZero();
    for (int i = 0; i < cavitySize; i++) {
		a(i,i) = cav.getTessArea(i);
		aInv(i,i) = 2 * M_PI / cav.getTessArea(i);
    }
	double fact = (epsilon+1.0)/(epsilon-1.0);
    PCMMatrix = (fact * aInv - DI) * a * SI;
    PCMMatrix = PCMMatrix.inverse();
    PCMMatrix *= (aInv - DI);
    PCMMatrix = PCMMatrix * a;
	builtIsotropicMatrix = true;
	builtAnisotropicMatrix = false;
}

VectorXd IEFSolver::compCharge(const VectorXd &potential) {
	VectorXd charge;
	if (builtIsotropicMatrix or builtAnisotropicMatrix) {
		charge = - PCMMatrix * potential;
	} else {
		cout << "PCM matrix not initialized" << endl;
		exit(1);
	}
	return charge;
}

void IEFSolver::compCharge(const VectorXd & potential, VectorXd & charge) {
	if (builtIsotropicMatrix or builtAnisotropicMatrix) {
		charge = - PCMMatrix * potential;
	} else {
		cout << "PCM matrix not initialized" << endl;
		exit(1);
	}
}

ostream & IEFSolver::printObject(ostream & os) {
	string type = "Traditional";
	os << "~~~~~~~~~~ PCMSolver ~~~~~~~~~~\n" << endl;
	os << "========== Solver section" << endl;
	os << "Solver Type: " << type << endl;
	os << "Solvent: " << solvent << endl;
	return os;
}

